"""
Complete workflow using internal MongoDB (opendb.mp_2022) instead of MP API.

Step 1: Query internal DB for carbon-rich compounds (no internet needed)
        Include compounds that are experimental OR within e_hull threshold
Step 2: Select best per system by highest C content
Step 3: Annotate melting points (curated lookup + empirical, downstream only)
Step 4: Export results + funnel stats

Example usage:
    # Full search: binary + ternary + bimetal for rare earths
    python run_screening_internal.py --mode all --metal-group rare_earth --partner-group nonmetal

    # Binary rare earth carbides only
    python run_screening_internal.py --mode binary --metal-group rare_earth

    # Specific systems
    python run_screening_internal.py --systems B-C-La C-Hf-Ta C-La
"""

import argparse
import json

import pandas as pd

from screen_carbon_rich_compounds_internal import (
    METAL_GROUPS, PARTNER_GROUPS, SEARCH_MODES,
    MONGO_URI, MONGO_DB, MONGO_COLLECTION,
    screen_systems, find_highest_carbon_per_system, build_output_name,
)
from predict_melting_point import annotate_with_melting_points, filter_by_melting_point


def main():
    parser = argparse.ArgumentParser(
        description="Screening workflow using internal MongoDB (no internet needed)"
    )
    parser.add_argument(
        "--systems", nargs="+", default=None,
        help="Explicit chemical systems (e.g., B-C-La C-Hf-Ta). Bypasses --mode"
    )
    parser.add_argument(
        "--mode", default="ternary", choices=SEARCH_MODES,
        help="Search mode: binary/ternary/bimetal/all (default: ternary)"
    )
    parser.add_argument("--metal-group", default=None, choices=list(METAL_GROUPS.keys()),
                        help="Primary metal element group")
    parser.add_argument("--partner-group", default=None, choices=list(PARTNER_GROUPS.keys()),
                        help="Partner element group for ternary/bimetal modes")
    parser.add_argument("--partner-elements", nargs="+", default=None,
                        help="Explicit partner elements (overrides --partner-group)")
    parser.add_argument("--min-c-fraction", type=float, default=0.25,
                        help="Minimum carbon atomic fraction (default: 0.25)")
    parser.add_argument("--max-ehull", type=float, default=0.1,
                        help="Maximum energy above hull in eV/atom (default: 0.1)")
    parser.add_argument("--mp-min", type=float, default=2000,
                        help="Minimum melting point filter in degrees C (default: 2000)")
    parser.add_argument("--mp-max", type=float, default=2500,
                        help="Maximum melting point filter in degrees C (default: 2500)")
    parser.add_argument("--mongo-uri", default=MONGO_URI,
                        help="MongoDB URI")
    parser.add_argument("--db-name", default=MONGO_DB,
                        help=f"Database name (default: {MONGO_DB})")
    parser.add_argument("--collection", default=MONGO_COLLECTION,
                        help=f"Collection name (default: {MONGO_COLLECTION})")
    parser.add_argument("--output-prefix", default=None,
                        help="Prefix for output files (auto-generated if not specified)")
    args = parser.parse_args()

    # Auto-generate prefix from search parameters
    prefix = args.output_prefix or build_output_name(
        args.mode, args.metal_group, args.partner_group,
        args.partner_elements, args.systems).replace(".csv", "")

    print("=" * 70)
    print("Graphitization Catalyst Candidate Screening (Internal DB)")
    print(f"Mode: {args.mode}" + (" (bypassed by --systems)" if args.systems else ""))
    print(f"Database: {args.db_name}.{args.collection}")
    print(f"Target: Highest-carbon compounds with mp in {args.mp_min}-{args.mp_max} C")
    print("=" * 70)

    # Step 1: Query internal MongoDB
    print("\n[Step 1/4] Querying internal MongoDB for carbon-rich compounds...")
    print("  Filter: experimental OR e_hull <= threshold")
    df, stats = screen_systems(
        systems=args.systems,
        mode=args.mode,
        metal_group=args.metal_group,
        partner_group=args.partner_group,
        partner_elements=args.partner_elements,
        min_carbon_fraction=args.min_c_fraction,
        max_energy_above_hull=args.max_ehull,
        mongo_uri=args.mongo_uri,
        db_name=args.db_name,
        collection_name=args.collection,
    )

    if df.empty:
        print("No compounds found. Try relaxing filters (--min-c-fraction, --max-ehull).")
        stats.print_funnel(args.min_c_fraction, args.max_ehull,
                           args.mp_min, args.mp_max)
        return

    print(f"Found {len(df)} carbon-rich compounds across {df['chemsys'].nunique()} systems")

    # Step 2: Select best per system (by highest C content only)
    print("\n[Step 2/4] Selecting best compound per system (highest C content)...")
    best = find_highest_carbon_per_system(df)

    # Step 3: Annotate with melting points (downstream, does NOT affect selection)
    print("\n[Step 3/4] Annotating melting points (curated lookup + empirical)...")
    df = annotate_with_melting_points(df)
    best = annotate_with_melting_points(best)

    # Step 4: Save and display results
    print("\n[Step 4/4] Saving results...")

    full_file = f"{prefix}_all.csv"
    df.to_csv(full_file, index=False)
    print(f"  All compounds: {full_file} ({len(df)} rows)")

    best_file = f"{prefix}_best_per_system.csv"
    best.to_csv(best_file, index=False)
    print(f"  Best per system: {best_file} ({len(best)} rows)")

    in_range = filter_by_melting_point(df, args.mp_min, args.mp_max)
    if not in_range.empty:
        range_file = f"{prefix}_in_mp_range.csv"
        in_range.to_csv(range_file, index=False)
        print(f"  In mp range: {range_file} ({len(in_range)} rows)")

    # Update stats with melting point info
    stats.mp_in_range = len(in_range)
    if "mp_source" in df.columns:
        stats.mp_source_counts = df["mp_source"].value_counts().to_dict()

    # Save funnel statistics
    stats_file = f"{prefix}_funnel_stats.json"
    stats_dict = stats.to_dict()
    stats_dict["parameters"] = {
        "mode": args.mode,
        "metal_group": args.metal_group,
        "partner_group": args.partner_group,
        "partner_elements": args.partner_elements,
        "min_c_fraction": args.min_c_fraction,
        "max_ehull": args.max_ehull,
        "mp_min": args.mp_min,
        "mp_max": args.mp_max,
        "database": f"{args.db_name}.{args.collection}",
    }
    with open(stats_file, "w") as f:
        json.dump(stats_dict, f, indent=2)
    print(f"  Funnel stats: {stats_file}")

    # Print funnel
    stats.print_funnel(args.min_c_fraction, args.max_ehull,
                       args.mp_min, args.mp_max)

    # Display results
    print("\n" + "=" * 70)
    print("ALL SCREENED COMPOUNDS (ranked by C content)")
    print("=" * 70)
    display_cols = [
        "formula", "chemsys", "C_atomic_fraction", "C_weight_fraction",
        "energy_above_hull_eV", "experimental", "melting_point_C", "mp_source",
        "material_id",
    ]
    cols = [c for c in display_cols if c in df.columns]
    print(df[cols].head(20).to_string(index=False))

    print("\n" + "=" * 70)
    print("BEST COMPOUND PER SYSTEM (highest C content)")
    print("=" * 70)
    cols_best = ["formula", "chemsys", "C_atomic_fraction", "energy_above_hull_eV",
                 "experimental", "melting_point_C", "mp_source"]
    cols_best = [c for c in cols_best if c in best.columns]
    print(best[cols_best].to_string(index=False))


if __name__ == "__main__":
    main()
