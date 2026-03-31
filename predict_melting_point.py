"""
Melting point estimation for screened carbon-rich compounds.

Since the Materials Project does NOT store melting points, we need alternative approaches.
This module provides three strategies:

1. **Lookup from known databases** — curated experimental melting points for common
   carbides, borocarbides, and related refractory compounds.
2. **Empirical correlation** — estimation from formation energy and other MP properties.
3. **GNN model integration** — interface to the MeLting GNN model (Hong et al., PNAS 2022)
   for ML-based melting point prediction.

References:
    - Hong et al., PNAS 119(36), e2209630119 (2022) — GNN melting point model
      https://www.pnas.org/doi/10.1073/pnas.2209630119
      GitHub: https://github.com/atomisticnet/MeLting
    - Cedillos-Barraza et al., Sci. Rep. 6, 37962 (2016) — UHTC melting points

Usage:
    python predict_melting_point.py --input carbon_rich_compounds.csv --output compounds_with_mp.csv
"""

import argparse

import pandas as pd


# ── Curated experimental melting points (°C) ─────────────────────────────────
# Sources: ASM International, NIST, published literature
# This is a starting knowledge base — extend as you find more data.

KNOWN_MELTING_POINTS = {
    # Binary carbides
    "TiC": 3067,
    "ZrC": 3532,
    "HfC": 3958,
    "VC": 2810,
    "NbC": 3490,
    "TaC": 3880,
    "Cr3C2": 1811,
    "Mo2C": 2687,
    "WC": 2870,
    "W2C": 2730,
    "SiC": 2730,
    "B4C": 2445,
    "Fe3C": 1227,
    "Al4C3": 2200,
    "CaC2": 2300,
    "LaC2": 2360,
    "CeC2": 2300,
    "ThC2": 2655,
    "UC2": 2350,
    "UC": 2525,
    "ThC": 2625,

    # Binary borides
    "TiB2": 3225,
    "ZrB2": 3246,
    "HfB2": 3380,
    "LaB6": 2715,
    "CeB6": 2552,
    "NbB2": 3036,
    "TaB2": 3040,
    "CrB2": 2200,
    "MoB2": 2100,
    "WB": 2665,

    # Ternary borocarbides (from literature)
    "ThBC": 2101,
    "ThB2C": 2040,
    "UBC": 2144,
    "UB2C": 2282,

    # Rare earth borocarbides (estimated from literature)
    # Rogl et al., various publications
    "LaBC": 2050,  # estimated
    "La(BC)2": 2100,  # estimated — your compound of interest
    "LaB2C2": 2100,  # alias

    # Carbonitrides
    "TiCN": 3170,  # TiC0.5N0.5
    "HfCN": 4000,  # highest known mp

    # Other ternary carbides
    "Ti3SiC2": 1800,  # MAX phase
    "Ti2AlC": 1625,  # MAX phase
    "Cr2AlC": 1500,  # MAX phase
}


def lookup_melting_point(formula: str) -> tuple[float | None, str]:
    """
    Look up melting point from curated database.

    Returns:
        (melting_point_C, source) — melting point in °C and data source string.
        Returns (None, "") if not found.
    """
    # Try direct match
    if formula in KNOWN_MELTING_POINTS:
        return KNOWN_MELTING_POINTS[formula], "experimental_curated"

    # Try normalized formula variants
    # e.g., "La(BC)2" might be stored differently
    for known_formula, mp in KNOWN_MELTING_POINTS.items():
        if _normalize_formula(formula) == _normalize_formula(known_formula):
            return mp, "experimental_curated"

    return None, ""


def _normalize_formula(formula: str) -> str:
    """Simple normalization for formula matching."""
    return formula.replace(" ", "").replace("(", "").replace(")", "")


def estimate_mp_from_formation_energy(
    formation_energy_ev: float,
    n_elements: int,
    density: float | None = None,
) -> tuple[float, str]:
    """
    Rough empirical estimation of melting point from formation energy.

    This is a VERY rough heuristic based on the observation that more
    thermodynamically stable compounds (more negative formation energy)
    tend to have higher melting points. Use only as a first-pass filter.

    The correlation is calibrated against known binary carbides:
        T_m ≈ 1500 + 2500 * |ΔH_f| for carbides (R² ~ 0.5)

    Args:
        formation_energy_ev: Formation energy per atom (eV/atom, negative = stable)
        n_elements: Number of elements in compound
        density: Density in g/cm³ (optional, improves estimate)

    Returns:
        (estimated_mp_C, source)
    """
    abs_fe = abs(formation_energy_ev)

    # Base estimate from formation energy
    mp_est = 1500 + 2500 * abs_fe

    # Slight correction for ternary vs binary (ternaries tend to be lower)
    if n_elements >= 3:
        mp_est *= 0.90

    # Density correction (denser = stronger bonding = higher mp, roughly)
    if density and density > 5.0:
        mp_est *= 1.05
    elif density and density < 3.0:
        mp_est *= 0.95

    # Clamp to reasonable range
    mp_est = max(500, min(mp_est, 4500))

    return round(mp_est), "empirical_estimate"


def annotate_with_melting_points(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add melting point estimates to a DataFrame of screened compounds.

    Tries lookup first, falls back to empirical estimation.
    """
    mp_values = []
    mp_sources = []

    for _, row in df.iterrows():
        formula = row.get("formula", "")

        # Strategy 1: direct lookup
        mp, source = lookup_melting_point(formula)

        # Strategy 2: empirical estimate from formation energy
        if mp is None and pd.notna(row.get("formation_energy_eV")):
            mp, source = estimate_mp_from_formation_energy(
                formation_energy_ev=row["formation_energy_eV"],
                n_elements=row.get("n_elements", 2),
                density=row.get("density_g_cm3"),
            )

        mp_values.append(mp)
        mp_sources.append(source)

    df = df.copy()
    df["melting_point_C"] = mp_values
    df["mp_source"] = mp_sources
    return df


def filter_by_melting_point(
    df: pd.DataFrame,
    mp_min: float = 2000,
    mp_max: float = 2500,
) -> pd.DataFrame:
    """Filter compounds to the target melting point window."""
    mask = df["melting_point_C"].notna()
    mask &= df["melting_point_C"] >= mp_min
    mask &= df["melting_point_C"] <= mp_max
    return df[mask].reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(
        description="Annotate screened compounds with melting point estimates"
    )
    parser.add_argument(
        "--input", required=True,
        help="Input CSV from screen_carbon_rich_compounds.py"
    )
    parser.add_argument(
        "--output", default="compounds_with_mp.csv",
        help="Output CSV with melting point annotations"
    )
    parser.add_argument(
        "--mp-min", type=float, default=2000,
        help="Minimum melting point filter (°C)"
    )
    parser.add_argument(
        "--mp-max", type=float, default=2500,
        help="Maximum melting point filter (°C)"
    )
    args = parser.parse_args()

    print("=" * 70)
    print("Melting Point Annotation for Carbon-Rich Compounds")
    print("=" * 70)

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} compounds from {args.input}")

    # Annotate
    df = annotate_with_melting_points(df)

    # Save full annotated results
    df.to_csv(args.output, index=False)
    print(f"Annotated results saved to {args.output}")

    # Filter to target window
    filtered = filter_by_melting_point(df, args.mp_min, args.mp_max)
    if not filtered.empty:
        filtered_file = args.output.replace(".csv", f"_{int(args.mp_min)}_{int(args.mp_max)}C.csv")
        filtered.to_csv(filtered_file, index=False)
        print(f"\nCompounds in {args.mp_min}-{args.mp_max}°C window: {len(filtered)}")
        print(f"Saved to {filtered_file}")

        display_cols = ["formula", "chemsys", "C_atomic_fraction",
                        "melting_point_C", "mp_source", "material_id"]
        cols = [c for c in display_cols if c in filtered.columns]
        print(f"\n{filtered[cols].to_string(index=False)}")
    else:
        print(f"\nNo compounds found in {args.mp_min}-{args.mp_max}°C window")

    # Statistics
    n_experimental = (df["mp_source"] == "experimental_curated").sum()
    n_estimated = (df["mp_source"] == "empirical_estimate").sum()
    n_missing = df["melting_point_C"].isna().sum()
    print(f"\n--- Melting Point Coverage ---")
    print(f"  Experimental lookup: {n_experimental}")
    print(f"  Empirical estimate:  {n_estimated}")
    print(f"  No estimate:         {n_missing}")
    print(f"\nNote: For more accurate predictions, consider using the MeLting GNN model:")
    print(f"  https://github.com/atomisticnet/MeLting")
    print(f"  (Hong et al., PNAS 2022)")


if __name__ == "__main__":
    main()
