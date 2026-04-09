"""
Complete workflow using Materials Project API (requires internet).

Step 1: Query MP for carbon-rich compounds
        Include compounds that are experimental OR within e_hull threshold
Step 2: Select best per system by highest C content
Step 3: Annotate melting points (curated lookup + empirical, downstream only)
Step 4: Export results + funnel stats

Example usage:
    python run_screening.py --api-key KEY --mode all --metal-group rare_earth --partner-group nonmetal
    python run_screening.py --api-key KEY --mode binary --metal-group rare_earth
    python run_screening.py --api-key KEY --systems La-B-C Hf-Ta-C La-C
"""

import argparse
import json

import pandas as pd

from screen_carbon_rich_compounds import (
    METAL_GROUPS, PARTNER_GROUPS, SEARCH_MODES,
    screen_systems, find_highest_carbon_per_system, build_output_name,
)
from predict_melting_point import annotate_with_melting_points, filter_by_melting_point


def main():
    parser = argparse.ArgumentParser(
        description="Screening workflow using Materials Project API"
    )
    parser.add_argument("--api-key", required=True, help="Materials Project API key")
    parser.add_argument("--systems", nargs="+", default=None,
                        help="Explicit chemical systems. Bypasses --mode")
    parser.add_argument("--mode", default="ternary", choices=SEARCH_MODES,
                        help="Search mode (default: ternary)")
    parser.add_argument("--metal-group", default=None, choices=list(METAL_GROUPS.keys()))
    parser.add_argument("--partner-group", default=None, choices=list(PARTNER_GROUPS.keys()))
    parser.add_argument("--partner-elements", nargs="+", default=None)
    parser.add_argument("--min-c-fraction", type=float, default=0.25)
    parser.add_argument("--max-ehull", type=float, default=0.1)
    parser.add_argument("--use-mapp", action="store_true",
                        help="Enable MAPP GNN melting point predictions")
    parser.add_argument("--mapp-url", default=None)
    parser.add_argument("--mp-min", type=float, default=2000)
    parser.add_argument("--mp-max", type=float, default=2500)
    parser.add_argument("--output-prefix", default=None,
                        help="Prefix for output files (auto-generated if not specified)")
    args = parser.parse_args()

    prefix = args.output_prefix or build_output_name(
        args.mode, args.metal_group, args.partner_group,
        args.partner_elements, args.systems).replace(".csv", "")

    print("=" * 70)
    print("Graphitization Catalyst Candidate Screening (MP API)")
    print(f"Mode: {args.mode}" + (" (bypassed by --systems)" if args.systems else ""))
    print(f"Target: Highest-carbon compounds with mp in {args.mp_min}-{args.mp_max} C")
    print("=" * 70)

    # Step 1: Query MP
    print("\n[Step 1/4] Querying Materials Project...")
    df, stats = screen_systems(
        api_key=args.api_key,
        systems=args.systems,
        mode=args.mode,
        metal_group=args.metal_group,
        partner_group=args.partner_group,
        partner_elements=args.partner_elements,
        min_carbon_fraction=args.min_c_fraction,
        max_energy_above_hull=args.max_ehull,
    )

    if df.empty:
        print("No compounds found.")
        stats.print_funnel(args.min_c_fraction, args.max_ehull, args.mp_min, args.mp_max)
        return

    print(f"Found {len(df)} carbon-rich compounds across {df['chemsys'].nunique()} systems")

    # Step 2: Best per system (highest C content)
    print("\n[Step 2/4] Selecting best compound per system (highest C content)...")
    best = find_highest_carbon_per_system(df)

    # Step 3: Annotate melting points (downstream)
    print("\n[Step 3/4] Annotating melting points...")
    mapp_kwargs = {}
    if args.use_mapp:
        mapp_kwargs["use_mapp"] = True
        if args.mapp_url:
            mapp_kwargs["mapp_url"] = args.mapp_url
    df = annotate_with_melting_points(df, **mapp_kwargs)
    best = annotate_with_melting_points(best, **mapp_kwargs)

    # Step 4: Save results
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

    stats.mp_in_range = len(in_range)
    if "mp_source" in df.columns:
        stats.mp_source_counts = df["mp_source"].value_counts().to_dict()

    stats_file = f"{prefix}_funnel_stats.json"
    with open(stats_file, "w") as f:
        json.dump(stats.to_dict(), f, indent=2)
    print(f"  Funnel stats: {stats_file}")

    stats.print_funnel(args.min_c_fraction, args.max_ehull, args.mp_min, args.mp_max)

    # Display
    print("\n" + "=" * 70)
    print("BEST COMPOUND PER SYSTEM (highest C content)")
    print("=" * 70)
    cols_best = ["formula", "chemsys", "C_atomic_fraction", "energy_above_hull_eV",
                 "experimental", "melting_point_C", "mp_source"]
    cols_best = [c for c in cols_best if c in best.columns]
    print(best[cols_best].to_string(index=False))


if __name__ == "__main__":
    main()
