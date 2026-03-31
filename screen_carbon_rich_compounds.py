"""
Screen carbon-rich compounds in ternary M-X-C systems using the Materials Project API.

Purpose: For graphitization catalyst research, find compounds with the highest carbon
content in specific ternary chemical systems (e.g., La-B-C → La(BC)2).
These are candidate phases that may form during high-temperature graphitization
and act as carbon-precipitation intermediates.

Usage:
    python screen_carbon_rich_compounds.py --api-key YOUR_MP_API_KEY
    python screen_carbon_rich_compounds.py --api-key YOUR_MP_API_KEY --systems La-B-C Fe-B-C
    python screen_carbon_rich_compounds.py --api-key YOUR_MP_API_KEY --metal-group rare_earth --nonmetal B

Requirements:
    pip install mp-api pymatgen pandas
"""

import argparse
import itertools
from typing import Optional

import pandas as pd
from mp_api.client import MPRester
from pymatgen.core import Composition


# ── Element groups for systematic screening ──────────────────────────────────

METAL_GROUPS = {
    "rare_earth": [
        "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd",
        "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Y", "Sc",
    ],
    "transition_3d": ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"],
    "transition_4d": ["Zr", "Nb", "Mo", "Ru", "Rh", "Pd"],
    "transition_5d": ["Hf", "Ta", "W", "Re", "Os", "Ir", "Pt"],
    "actinide": ["Th", "U"],
    "alkaline_earth": ["Ca", "Sr", "Ba"],
    "alkali": ["Li", "Na", "K"],
}

# Non-metal partners commonly forming ternary compounds with C
NONMETALS = ["B", "N", "Si", "P", "S"]


def query_chemsys(
    mpr: MPRester,
    chemsys: str,
    min_carbon_fraction: float = 0.2,
    max_energy_above_hull: float = 0.1,  # eV/atom, for thermodynamic stability
    experimental_only: bool = False,
) -> list[dict]:
    """
    Query all compounds in a chemical system, filter for C-containing phases.

    Args:
        mpr: MPRester client instance
        chemsys: Chemical system string, e.g. "La-B-C"
        min_carbon_fraction: Minimum atomic fraction of carbon (0-1)
        max_energy_above_hull: Maximum energy above hull in eV/atom (stability filter)
        experimental_only: If True, only return experimentally synthesized compounds
                           (excludes theoretical/predicted structures)

    Returns:
        List of dicts with compound info, sorted by carbon fraction descending
    """
    elements = chemsys.split("-")
    if "C" not in elements:
        raise ValueError(f"System {chemsys} does not contain carbon")

    search_kwargs = dict(
        chemsys=chemsys,
        energy_above_hull=(0, max_energy_above_hull),
        fields=[
            "material_id",
            "formula_pretty",
            "composition_reduced",
            "energy_above_hull",
            "formation_energy_per_atom",
            "symmetry",
            "volume",
            "density",
            "nelements",
            "theoretical",
        ],
    )
    if experimental_only:
        search_kwargs["theoretical"] = False

    try:
        docs = mpr.materials.summary.search(**search_kwargs)
    except Exception as e:
        print(f"  Warning: failed to query {chemsys}: {e}")
        return []

    results = []
    for doc in docs:
        comp = doc.composition_reduced
        c_frac = comp.get_atomic_fraction("C")

        # Must contain carbon and meet minimum fraction
        if c_frac < min_carbon_fraction:
            continue

        # Must be a true ternary (or at least contain elements beyond just C)
        n_elements = len(comp.elements)
        if n_elements < 2:
            continue

        # Compute carbon weight fraction as well
        c_wt_frac = comp.get_wt_fraction("C")

        spacegroup = ""
        crystal_system = ""
        if hasattr(doc, "symmetry") and doc.symmetry:
            spacegroup = getattr(doc.symmetry, "symbol", "")
            crystal_system = getattr(doc.symmetry, "crystal_system", "")

        results.append({
            "material_id": str(doc.material_id),
            "formula": doc.formula_pretty,
            "chemsys": chemsys,
            "C_atomic_fraction": round(c_frac, 4),
            "C_weight_fraction": round(c_wt_frac, 4),
            "n_elements": n_elements,
            "energy_above_hull_eV": round(doc.energy_above_hull, 4) if doc.energy_above_hull else None,
            "formation_energy_eV": round(doc.formation_energy_per_atom, 4) if doc.formation_energy_per_atom else None,
            "spacegroup": spacegroup,
            "crystal_system": crystal_system,
            "density_g_cm3": round(doc.density, 2) if doc.density else None,
            "theoretical": getattr(doc, "theoretical", None),
        })

    results.sort(key=lambda x: x["C_atomic_fraction"], reverse=True)
    return results


def screen_systems(
    api_key: str,
    systems: Optional[list[str]] = None,
    metal_group: Optional[str] = None,
    nonmetal: Optional[str] = None,
    min_carbon_fraction: float = 0.2,
    max_energy_above_hull: float = 0.1,
    experimental_only: bool = False,
) -> pd.DataFrame:
    """
    Screen multiple ternary M-X-C systems for carbon-rich compounds.

    If `systems` is provided, query those specific systems.
    Otherwise, generate M-X-C combinations from metal_group and nonmetal.
    """
    if systems is None:
        if metal_group is None or nonmetal is None:
            raise ValueError("Provide either --systems or both --metal-group and --nonmetal")

        metals = METAL_GROUPS.get(metal_group)
        if metals is None:
            raise ValueError(f"Unknown metal group: {metal_group}. Options: {list(METAL_GROUPS.keys())}")

        nonmetals = [nonmetal] if nonmetal != "all" else NONMETALS
        systems = [f"{m}-{x}-C" for m, x in itertools.product(metals, nonmetals)]

    all_results = []
    with MPRester(api_key) as mpr:
        for i, sys in enumerate(systems):
            print(f"[{i+1}/{len(systems)}] Querying {sys} ...")
            results = query_chemsys(mpr, sys, min_carbon_fraction, max_energy_above_hull, experimental_only)
            if results:
                print(f"  Found {len(results)} carbon-rich compounds")
                all_results.extend(results)
            else:
                print(f"  No carbon-rich compounds found")

    df = pd.DataFrame(all_results)
    if not df.empty:
        df = df.sort_values("C_atomic_fraction", ascending=False).reset_index(drop=True)
    return df


def find_highest_carbon_per_system(df: pd.DataFrame) -> pd.DataFrame:
    """For each chemsys, return only the compound with highest C atomic fraction."""
    if df.empty:
        return df
    idx = df.groupby("chemsys")["C_atomic_fraction"].idxmax()
    return df.loc[idx].sort_values("C_atomic_fraction", ascending=False).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(
        description="Screen ternary M-X-C systems for carbon-rich compounds via Materials Project API"
    )
    parser.add_argument("--api-key", required=True, help="Materials Project API key")
    parser.add_argument(
        "--systems", nargs="+", default=None,
        help="Specific chemical systems to query, e.g. La-B-C Fe-B-C"
    )
    parser.add_argument(
        "--metal-group", default=None,
        choices=list(METAL_GROUPS.keys()),
        help="Metal group for systematic screening"
    )
    parser.add_argument(
        "--nonmetal", default=None,
        help=f"Non-metal partner element (or 'all' for {NONMETALS})"
    )
    parser.add_argument(
        "--min-c-fraction", type=float, default=0.2,
        help="Minimum carbon atomic fraction (default: 0.2)"
    )
    parser.add_argument(
        "--max-ehull", type=float, default=0.1,
        help="Maximum energy above hull in eV/atom (default: 0.1, set higher for metastable phases)"
    )
    parser.add_argument(
        "--experimental-only", action="store_true",
        help="Only include experimentally synthesized compounds (exclude theoretical/predicted structures)"
    )
    parser.add_argument(
        "--output", default="carbon_rich_compounds.csv",
        help="Output CSV filename"
    )
    args = parser.parse_args()

    print("=" * 70)
    print("Carbon-Rich Compound Screening via Materials Project")
    print("=" * 70)

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
        print("\nNo compounds found matching criteria.")
        return

    # Save full results
    df.to_csv(args.output, index=False)
    print(f"\nFull results saved to {args.output} ({len(df)} compounds)")

    # Show the highest-carbon compound per system
    best = find_highest_carbon_per_system(df)
    best_file = args.output.replace(".csv", "_best_per_system.csv")
    best.to_csv(best_file, index=False)
    print(f"Best per system saved to {best_file} ({len(best)} systems)")

    print("\n" + "=" * 70)
    print("Top carbon-rich compounds (highest C fraction per system):")
    print("=" * 70)
    display_cols = ["formula", "chemsys", "C_atomic_fraction", "C_weight_fraction",
                    "energy_above_hull_eV", "material_id"]
    print(best[display_cols].to_string(index=False))

    # Summary statistics
    print(f"\n--- Summary ---")
    print(f"Total compounds found: {len(df)}")
    print(f"Unique systems with hits: {df['chemsys'].nunique()}")
    print(f"Highest C atomic fraction: {df['C_atomic_fraction'].max():.4f} "
          f"({df.loc[df['C_atomic_fraction'].idxmax(), 'formula']})")


if __name__ == "__main__":
    main()
