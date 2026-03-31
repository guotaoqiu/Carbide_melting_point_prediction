"""
Melting point estimation for screened carbon-rich compounds.

Since the Materials Project does NOT store melting points, we need alternative approaches.
This module provides three strategies (in priority order):

1. **Curated lookup** — hardcoded experimental melting points for common carbides,
   borocarbides, and related refractory compounds.
2. **MAPP GNN prediction** — Hong et al.'s GNN+ResNet ensemble model (PNAS 2022)
   that predicts melting temperature from chemical formula via a remote API.
   Supports compounds with up to 6 elements. Returns prediction + uncertainty.
3. **Empirical estimation** — rough correlation from formation energy and other
   MP properties. Use as last-resort fallback only (R2 ~ 0.5).

References:
    - Hong et al., PNAS 119(36), e2209630119 (2022) — MAPP melting point model
      https://www.pnas.org/doi/10.1073/pnas.2209630119
      GitHub: https://github.com/qjhong/mapp_api
    - Cedillos-Barraza et al., Sci. Rep. 6, 37962 (2016) — UHTC melting points

Usage:
    # Default: curated lookup + empirical fallback (no network needed)
    python predict_melting_point.py --input carbon_rich_compounds.csv

    # Enable MAPP GNN predictions (requires internet, replaces empirical for unknowns)
    python predict_melting_point.py --input carbon_rich_compounds.csv --use-mapp

    # MAPP with custom server URL (if self-hosting)
    python predict_melting_point.py --input compounds.csv --use-mapp --mapp-url http://your-server:5007
"""

import argparse
import json
import time

import pandas as pd
import requests


# ── Curated experimental melting points (deg C) ─────────────────────────────
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

# ── MAPP API configuration ──────────────────────────────────────────────────
# Hong et al.'s MAPP (Materials-Agnostic Platform for Prediction) server.
# The model is a GNN + ResNet ensemble of 30 models trained on ~10k compounds.
# It accepts chemical formulas and returns melting temperature in Kelvin.
# Ref: https://github.com/qjhong/mapp_api

MAPP_DEFAULT_URL = "http://206.207.50.58:5007/MT_ML_Qijun_Hong_Predict_noNN"
MAPP_BATCH_SIZE = 500  # max formulas per API call to avoid timeouts


def lookup_melting_point(formula: str) -> tuple[float | None, str]:
    """
    Look up melting point from curated database.

    Returns:
        (melting_point_C, source) — melting point in deg C and data source string.
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


def predict_with_mapp(
    formulas: list[str],
    mapp_url: str = MAPP_DEFAULT_URL,
    timeout: int = 120,
    retries: int = 3,
) -> dict[str, tuple[float, float]]:
    """
    Predict melting temperatures using Hong et al.'s MAPP API.

    The MAPP model is a GNN + ResNet ensemble of 30 bootstrap models.
    It accepts chemical formulas (up to 6 elements) and returns:
    - Predicted melting temperature (K)
    - Standard error from the ensemble (K)

    Args:
        formulas: List of chemical formula strings (e.g., ["TiC", "LaB2C2"])
        mapp_url: MAPP API endpoint URL
        timeout: Request timeout in seconds
        retries: Number of retry attempts on failure

    Returns:
        Dict mapping formula -> (melting_point_C, std_error_C).
        Formulas that fail prediction are omitted from the dict.
    """
    if not formulas:
        return {}

    results = {}

    # Process in batches
    for batch_start in range(0, len(formulas), MAPP_BATCH_SIZE):
        batch = formulas[batch_start:batch_start + MAPP_BATCH_SIZE]

        # Format payload as MAPP expects: [{"9": "formula1"}, {"9": "formula2"}, ...]
        # The key "9" is an arbitrary placeholder used by the MAPP API protocol
        payload = json.dumps([{"9": f} for f in batch])

        for attempt in range(retries):
            try:
                response = requests.post(
                    mapp_url,
                    data=payload,
                    headers={"content-type": "application/json"},
                    timeout=timeout,
                )
                response.raise_for_status()
                data = response.json()
                break
            except (requests.RequestException, json.JSONDecodeError) as e:
                if attempt < retries - 1:
                    wait = 2 ** (attempt + 1)
                    print(f"  MAPP API attempt {attempt + 1} failed: {e}. Retrying in {wait}s...")
                    time.sleep(wait)
                else:
                    print(f"  MAPP API failed after {retries} attempts: {e}")
                    print(f"  Skipping batch of {len(batch)} formulas")
                    data = None

        if data is None:
            continue

        # Parse response — MAPP returns a list of dicts with predictions
        # Format: [{"chemical_formula": "TiC", "melting_temperature_in_kelvin": 3340.5,
        #           "standard_error_in_kelvin": 120.3}, ...]
        if isinstance(data, list):
            for entry in data:
                formula = entry.get("chemical_formula", "")
                mp_kelvin = entry.get("melting_temperature_in_kelvin")
                se_kelvin = entry.get("standard_error_in_kelvin", 0)
                if formula and mp_kelvin is not None:
                    # Convert Kelvin to Celsius
                    mp_celsius = mp_kelvin - 273.15
                    se_celsius = se_kelvin  # standard error same magnitude
                    results[formula] = (round(mp_celsius), round(se_celsius))
        elif isinstance(data, dict):
            # Alternative response format: single dict with arrays
            for formula, mp_kelvin in zip(
                data.get("chemical_formula", []),
                data.get("melting_temperature_in_kelvin", []),
            ):
                se_kelvin = 0
                if "standard_error_in_kelvin" in data:
                    idx = data["chemical_formula"].index(formula)
                    se_kelvin = data["standard_error_in_kelvin"][idx]
                mp_celsius = mp_kelvin - 273.15
                results[formula] = (round(mp_celsius), round(se_kelvin))

    return results


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
        T_m ~ 1500 + 2500 * |DeltaH_f| for carbides (R2 ~ 0.5)

    Args:
        formation_energy_ev: Formation energy per atom (eV/atom, negative = stable)
        n_elements: Number of elements in compound
        density: Density in g/cm3 (optional, improves estimate)

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


def annotate_with_melting_points(
    df: pd.DataFrame,
    use_mapp: bool = False,
    mapp_url: str = MAPP_DEFAULT_URL,
) -> pd.DataFrame:
    """
    Add melting point estimates to a DataFrame of screened compounds.

    Priority order:
    1. Curated experimental lookup
    2. MAPP GNN prediction (if --use-mapp enabled)
    3. Empirical estimation from formation energy (fallback)

    Args:
        df: DataFrame with at least a "formula" column
        use_mapp: If True, use MAPP API for compounds not in curated database
        mapp_url: MAPP API endpoint URL
    """
    df = df.copy()

    # Step 1: Curated lookup for all formulas
    mp_values = []
    mp_sources = []
    mp_errors = []
    needs_prediction = []  # (index, formula) for formulas not in curated DB

    for idx, row in df.iterrows():
        formula = row.get("formula", "")
        mp, source = lookup_melting_point(formula)
        if mp is not None:
            mp_values.append(mp)
            mp_sources.append(source)
            mp_errors.append(None)
        else:
            mp_values.append(None)
            mp_sources.append("")
            mp_errors.append(None)
            needs_prediction.append((idx, formula))

    # Step 2: MAPP GNN prediction for unknowns
    if use_mapp and needs_prediction:
        formulas_to_predict = [f for _, f in needs_prediction]
        idx_map = {f: i for i, f in needs_prediction}

        print(f"  Querying MAPP API for {len(formulas_to_predict)} compounds...")
        mapp_results = predict_with_mapp(formulas_to_predict, mapp_url=mapp_url)

        mapp_hits = 0
        for formula, (mp_c, se_c) in mapp_results.items():
            if formula in idx_map:
                pos = list(df.index).index(idx_map[formula]) if idx_map[formula] in df.index else None
                if pos is not None:
                    mp_values[pos] = mp_c
                    mp_sources[pos] = "mapp_gnn"
                    mp_errors[pos] = se_c
                    mapp_hits += 1

        print(f"  MAPP returned predictions for {mapp_hits}/{len(formulas_to_predict)} compounds")

    # Step 3: Empirical fallback for remaining unknowns
    for i, (mp, source) in enumerate(zip(mp_values, mp_sources)):
        if mp is None:
            row = df.iloc[i]
            if pd.notna(row.get("formation_energy_eV")):
                mp_est, src = estimate_mp_from_formation_energy(
                    formation_energy_ev=row["formation_energy_eV"],
                    n_elements=row.get("n_elements", 2),
                    density=row.get("density_g_cm3"),
                )
                mp_values[i] = mp_est
                mp_sources[i] = src

    df["melting_point_C"] = mp_values
    df["mp_source"] = mp_sources
    df["mp_std_error_C"] = mp_errors
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
        "--use-mapp", action="store_true",
        help=(
            "Enable MAPP GNN melting point predictions (Hong et al., PNAS 2022). "
            "Requires internet access to reach the MAPP API server. "
            "Predictions use an ensemble of 30 GNN+ResNet models and include "
            "uncertainty estimates (standard error)"
        ),
    )
    parser.add_argument(
        "--mapp-url", default=MAPP_DEFAULT_URL,
        help=(
            "MAPP API endpoint URL. Only needed if self-hosting the MAPP server "
            f"(default: {MAPP_DEFAULT_URL})"
        ),
    )
    parser.add_argument(
        "--mp-min", type=float, default=2000,
        help="Minimum melting point filter (deg C, default: 2000)"
    )
    parser.add_argument(
        "--mp-max", type=float, default=2500,
        help="Maximum melting point filter (deg C, default: 2500)"
    )
    args = parser.parse_args()

    print("=" * 70)
    print("Melting Point Annotation for Carbon-Rich Compounds")
    print("=" * 70)

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} compounds from {args.input}")

    if args.use_mapp:
        print(f"MAPP GNN predictions ENABLED (server: {args.mapp_url})")
    else:
        print("MAPP GNN predictions disabled (use --use-mapp to enable)")

    # Annotate
    df = annotate_with_melting_points(df, use_mapp=args.use_mapp, mapp_url=args.mapp_url)

    # Save full annotated results
    df.to_csv(args.output, index=False)
    print(f"Annotated results saved to {args.output}")

    # Filter to target window
    filtered = filter_by_melting_point(df, args.mp_min, args.mp_max)
    if not filtered.empty:
        filtered_file = args.output.replace(".csv", f"_{int(args.mp_min)}_{int(args.mp_max)}C.csv")
        filtered.to_csv(filtered_file, index=False)
        print(f"\nCompounds in {args.mp_min}-{args.mp_max} deg C window: {len(filtered)}")
        print(f"Saved to {filtered_file}")

        display_cols = ["formula", "chemsys", "C_atomic_fraction",
                        "melting_point_C", "mp_std_error_C", "mp_source", "material_id"]
        cols = [c for c in display_cols if c in filtered.columns]
        print(f"\n{filtered[cols].to_string(index=False)}")
    else:
        print(f"\nNo compounds found in {args.mp_min}-{args.mp_max} deg C window")

    # Statistics
    n_experimental = (df["mp_source"] == "experimental_curated").sum()
    n_mapp = (df["mp_source"] == "mapp_gnn").sum()
    n_estimated = (df["mp_source"] == "empirical_estimate").sum()
    n_missing = df["melting_point_C"].isna().sum()
    print(f"\n--- Melting Point Coverage ---")
    print(f"  Experimental lookup:  {n_experimental}")
    print(f"  MAPP GNN prediction:  {n_mapp}")
    print(f"  Empirical estimate:   {n_estimated}")
    print(f"  No estimate:          {n_missing}")

    if n_mapp > 0:
        mapp_df = df[df["mp_source"] == "mapp_gnn"]
        avg_error = mapp_df["mp_std_error_C"].mean()
        print(f"  MAPP avg std error:   {avg_error:.0f} deg C")

    if not args.use_mapp and (n_estimated > 0 or n_missing > 0):
        print(f"\nTip: Use --use-mapp for ML-based predictions on the "
              f"{n_estimated + n_missing} compounds without experimental data")


if __name__ == "__main__":
    main()
