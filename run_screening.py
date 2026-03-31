"""
Complete workflow: Screen carbon-rich compounds → Annotate melting points → Rank candidates.

This is the main entry point that combines:
1. Materials Project API querying for carbon-rich compounds
2. Melting point estimation/lookup
3. Ranking candidates for graphitization catalyst research

Example usage:
    # Screen all rare earth borocarbides (RE-B-C systems)
    python run_screening.py --api-key YOUR_KEY --metal-group rare_earth --nonmetal B

    # Screen specific systems
    python run_screening.py --api-key YOUR_KEY --systems La-B-C Ce-B-C Y-B-C Sc-B-C

    # Screen all rare earth systems with B, N, Si
    python run_screening.py --api-key YOUR_KEY --metal-group rare_earth --nonmetal all

    # Broader search: include metastable phases
    python run_screening.py --api-key YOUR_KEY --metal-group rare_earth --nonmetal B --max-ehull 0.3
"""

import argparse
from datetime import datetime

import pandas as pd

from screen_carbon_rich_compounds import METAL_GROUPS, NONMETALS, screen_systems, find_highest_carbon_per_system
from predict_melting_point import annotate_with_melting_points, filter_by_melting_point


def rank_candidates(
    df: pd.DataFrame,
    mp_min: float = 2000,
    mp_max: float = 2500,
    weight_c_fraction: float = 0.5,
    weight_stability: float = 0.3,
    weight_mp_proximity: float = 0.2,
) -> pd.DataFrame:
    """
    Rank compounds by a composite score considering:
    - Carbon atomic fraction (higher = better for graphitization)
    - Thermodynamic stability (lower e_above_hull = better)
    - Melting point proximity to target range center

    The ideal compound has:
    - High carbon content (more C available for precipitation)
    - Good thermodynamic stability (will actually form at temperature)
    - Melting point in the 2000-2500°C range (matches process conditions)
    """
    df = df.copy()

    # Normalize carbon fraction to [0, 1]
    c_max = df["C_atomic_fraction"].max()
    c_min = df["C_atomic_fraction"].min()
    if c_max > c_min:
        df["score_c"] = (df["C_atomic_fraction"] - c_min) / (c_max - c_min)
    else:
        df["score_c"] = 1.0

    # Normalize stability (invert: lower e_hull = higher score)
    if "energy_above_hull_eV" in df.columns:
        e_max = df["energy_above_hull_eV"].max()
        if e_max > 0:
            df["score_stability"] = 1.0 - df["energy_above_hull_eV"] / e_max
        else:
            df["score_stability"] = 1.0
    else:
        df["score_stability"] = 0.5

    # Melting point proximity score (1.0 at center of range, 0 at edges)
    mp_center = (mp_min + mp_max) / 2
    mp_range = (mp_max - mp_min) / 2
    if "melting_point_C" in df.columns:
        df["score_mp"] = df["melting_point_C"].apply(
            lambda x: max(0, 1.0 - abs(x - mp_center) / mp_range) if pd.notna(x) else 0.0
        )
    else:
        df["score_mp"] = 0.0

    # Composite score
    df["composite_score"] = (
        weight_c_fraction * df["score_c"]
        + weight_stability * df["score_stability"]
        + weight_mp_proximity * df["score_mp"]
    )

    return df.sort_values("composite_score", ascending=False).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(
        description="Complete screening workflow for graphitization catalyst candidates"
    )
    parser.add_argument("--api-key", required=True, help="Materials Project API key")
    parser.add_argument("--systems", nargs="+", default=None)
    parser.add_argument("--metal-group", default=None, choices=list(METAL_GROUPS.keys()))
    parser.add_argument("--nonmetal", default=None)
    parser.add_argument("--min-c-fraction", type=float, default=0.25)
    parser.add_argument("--max-ehull", type=float, default=0.1)
    parser.add_argument("--experimental-only", action="store_true",
                        help="Only include experimentally synthesized compounds")
    parser.add_argument("--mp-min", type=float, default=2000)
    parser.add_argument("--mp-max", type=float, default=2500)
    parser.add_argument("--output-prefix", default=None, help="Prefix for output files")
    args = parser.parse_args()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    prefix = args.output_prefix or f"screening_{timestamp}"

    print("=" * 70)
    print("Graphitization Catalyst Candidate Screening")
    print(f"Target: Highest-carbon compounds with mp in {args.mp_min}-{args.mp_max}°C")
    print("=" * 70)

    # Step 1: Query Materials Project
    print("\n[Step 1/4] Querying Materials Project for carbon-rich compounds...")
    df = screen_systems(
        api_key=args.api_key,
        systems=args.systems,
        metal_group=args.metal_group,
        nonmetal=args.nonmetal,
        min_carbon_fraction=args.min_c_fraction,
        max_energy_above_hull=args.max_ehull,
        experimental_only=args.experimental_only,
    )

    if df.empty:
        print("No compounds found. Try relaxing filters (--min-c-fraction, --max-ehull).")
        return

    print(f"Found {len(df)} carbon-rich compounds across {df['chemsys'].nunique()} systems")

    # Step 2: Annotate with melting points
    print("\n[Step 2/4] Annotating melting points...")
    df = annotate_with_melting_points(df)

    # Step 3: Rank candidates
    print("\n[Step 3/4] Ranking candidates...")
    df = rank_candidates(df, args.mp_min, args.mp_max)

    # Step 4: Save and display results
    print("\n[Step 4/4] Saving results...")

    # Full results
    full_file = f"{prefix}_all.csv"
    df.to_csv(full_file, index=False)
    print(f"  All compounds: {full_file} ({len(df)} rows)")

    # Best per system
    best = find_highest_carbon_per_system(df)
    best = annotate_with_melting_points(best)
    best = rank_candidates(best, args.mp_min, args.mp_max)
    best_file = f"{prefix}_best_per_system.csv"
    best.to_csv(best_file, index=False)
    print(f"  Best per system: {best_file} ({len(best)} rows)")

    # Filtered to mp range
    in_range = filter_by_melting_point(df, args.mp_min, args.mp_max)
    if not in_range.empty:
        range_file = f"{prefix}_in_mp_range.csv"
        in_range.to_csv(range_file, index=False)
        print(f"  In mp range: {range_file} ({len(in_range)} rows)")

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

    # Best per system summary
    print("\n" + "=" * 70)
    print("HIGHEST CARBON COMPOUND PER SYSTEM")
    print("=" * 70)
    cols_best = ["formula", "chemsys", "C_atomic_fraction", "melting_point_C", "mp_source"]
    cols_best = [c for c in cols_best if c in best.columns]
    print(best[cols_best].to_string(index=False))

    print("\n" + "=" * 70)
    print("NEXT STEPS")
    print("=" * 70)
    print("""
1. For compounds with mp_source='empirical_estimate', verify with:
   - MeLting GNN model: https://github.com/atomisticnet/MeLting
   - Literature search on specific compounds
   - CALPHAD databases (if available)

2. For promising candidates, check:
   - Phase diagram to confirm stability at target temperature
   - Carbon precipitation behavior (does C exsolve on cooling?)
   - Compatibility with your carbon matrix material

3. Consider the 2025 borocarbide screening paper for M-B-C systems:
   ML-guided search for energetically favorable metal borocarbide ternary compounds
   (Journal of Alloys and Compounds, 2025)
""")


if __name__ == "__main__":
    main()
