"""
Complete workflow using internal MongoDB (opendb.mp_2022) instead of MP API.

Step 1: Query internal DB for carbon-rich compounds (no internet needed)
Step 2: Annotate melting points (curated lookup + empirical estimate)
Step 3: Rank candidates
Step 4: Export results + funnel stats

For MAPP GNN melting point predictions, export the CSV and run predict_melting_point.py
with --use-mapp on a machine with internet access.

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
from datetime import datetime

import pandas as pd

from screen_carbon_rich_compounds_internal import (
    METAL_GROUPS, PARTNER_GROUPS, SEARCH_MODES,
    MONGO_URI, MONGO_DB, MONGO_COLLECTION,
    screen_systems, find_highest_carbon_per_system,
)
from predict_melting_point import annotate_with_melting_points, filter_by_melting_point


def rank_candidates(
    df: pd.DataFrame,
    mp_min: float = 2000,
    mp_max: float = 2500,
    weight_c_fraction: float = 0.5,
    weight_stability: float = 0.3,
    weight_mp_proximity: float = 0.2,
) -> pd.DataFrame:
    """Rank compounds by composite score (C content + stability + mp proximity)."""
    df = df.copy()

    c_max = df["C_atomic_fraction"].max()
    c_min = df["C_atomic_fraction"].min()
    if c_max > c_min:
        df["score_c"] = (df["C_atomic_fraction"] - c_min) / (c_max - c_min)
    else:
        df["score_c"] = 1.0

    if "energy_above_hull_eV" in df.columns:
        e_max = df["energy_above_hull_eV"].max()
        if e_max > 0:
            df["score_stability"] = 1.0 - df["energy_above_hull_eV"] / e_max
        else:
            df["score_stability"] = 1.0
    else:
        df["score_stability"] = 0.5

    mp_center = (mp_min + mp_max) / 2
    mp_range = (mp_max - mp_min) / 2
    if "melting_point_C" in df.columns:
        df["score_mp"] = df["melting_point_C"].apply(
            lambda x: max(0, 1.0 - abs(x - mp_center) / mp_range) if pd.notna(x) else 0.0
        )
    else:
        df["score_mp"] = 0.0

    df["composite_score"] = (
        weight_c_fraction * df["score_c"]
        + weight_stability * df["score_stability"]
        + weight_mp_proximity * df["score_mp"]
    )

    return df.sort_values("composite_score", ascending=False).reset_index(drop=True)


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
    parser.add_argument("--experimental-only", action="store_true",
                        help="Only include experimentally synthesized compounds")
    parser.add_argument("--mp-min", type=float, default=2000,
                        help="Minimum melting point filter in degrees C (default: 2000)")
    parser.add_argument("--mp-max", type=float, default=2500,
                        help="Maximum melting point filter in degrees C (default: 2500)")
    parser.add_argument("--mongo-uri", default=MONGO_URI,
                        help=f"MongoDB URI (default: {MONGO_URI})")
    parser.add_argument("--db-name", default=MONGO_DB,
                        help=f"Database name (default: {MONGO_DB})")
    parser.add_argument("--collection", default=MONGO_COLLECTION,
                        help=f"Collection name (default: {MONGO_COLLECTION})")
    parser.add_argument("--output-prefix", default=None, help="Prefix for output files")
    args = parser.parse_args()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    prefix = args.output_prefix or f"screening_{timestamp}"

    print("=" * 70)
    print("Graphitization Catalyst Candidate Screening (Internal DB)")
    print(f"Mode: {args.mode}" + (" (bypassed by --systems)" if args.systems else ""))
    print(f"Database: {args.db_name}.{args.collection}")
    print(f"Target: Highest-carbon compounds with mp in {args.mp_min}-{args.mp_max} C")
    print("=" * 70)

    # Step 1: Query internal MongoDB
    print("\n[Step 1/4] Querying internal MongoDB for carbon-rich compounds...")
    df, stats = screen_systems(
        systems=args.systems,
        mode=args.mode,
        metal_group=args.metal_group,
        partner_group=args.partner_group,
        partner_elements=args.partner_elements,
        min_carbon_fraction=args.min_c_fraction,
        max_energy_above_hull=args.max_ehull,
        experimental_only=args.experimental_only,
        mongo_uri=args.mongo_uri,
        db_name=args.db_name,
        collection_name=args.collection,
    )

    if df.empty:
        print("No compounds found. Try relaxing filters (--min-c-fraction, --max-ehull).")
        stats.print_funnel(args.min_c_fraction, args.max_ehull,
                           args.mp_min, args.mp_max, args.experimental_only)
        return

    print(f"Found {len(df)} carbon-rich compounds across {df['chemsys'].nunique()} systems")

    # Step 2: Annotate with melting points (curated + empirical, no MAPP here)
    print("\n[Step 2/4] Annotating melting points (curated lookup + empirical)...")
    df = annotate_with_melting_points(df)

    # Step 3: Rank candidates
    print("\n[Step 3/4] Ranking candidates...")
    df = rank_candidates(df, args.mp_min, args.mp_max)

    # Step 4: Save and display results
    print("\n[Step 4/4] Saving results...")

    full_file = f"{prefix}_all.csv"
    df.to_csv(full_file, index=False)
    print(f"  All compounds: {full_file} ({len(df)} rows)")

    best = find_highest_carbon_per_system(df)
    best = annotate_with_melting_points(best)
    best = rank_candidates(best, args.mp_min, args.mp_max)
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
        "experimental_only": args.experimental_only,
        "mp_min": args.mp_min,
        "mp_max": args.mp_max,
        "database": f"{args.db_name}.{args.collection}",
    }
    with open(stats_file, "w") as f:
        json.dump(stats_dict, f, indent=2)
    print(f"  Funnel stats: {stats_file}")

    # Print funnel
    stats.print_funnel(args.min_c_fraction, args.max_ehull,
                       args.mp_min, args.mp_max, args.experimental_only)

    # Display top candidates
    print("\n" + "=" * 70)
    print("TOP CANDIDATES (ranked by composite score)")
    print("=" * 70)
    display_cols = [
        "formula", "chemsys", "C_atomic_fraction", "C_weight_fraction",
        "melting_point_C", "mp_source", "energy_above_hull_eV",
        "composite_score", "material_id",
    ]
    cols = [c for c in display_cols if c in df.columns]
    print(df[cols].head(20).to_string(index=False))

    # Best per system
    print("\n" + "=" * 70)
    print("HIGHEST CARBON COMPOUND PER SYSTEM")
    print("=" * 70)
    cols_best = ["formula", "chemsys", "C_atomic_fraction", "melting_point_C", "mp_source"]
    cols_best = [c for c in cols_best if c in best.columns]
    print(best[cols_best].to_string(index=False))

    print("\n" + "=" * 70)
    print("NEXT STEPS")
    print("=" * 70)
    print(f"""
To get MAPP GNN melting point predictions, take the CSV to a machine
with internet access and run:

    python predict_melting_point.py --input {full_file} --use-mapp

This will query the MAPP API for ML-based melting point predictions
for all compounds without curated experimental values.
""")


if __name__ == "__main__":
    main()
