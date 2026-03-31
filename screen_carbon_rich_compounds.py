"""
Screen carbon-rich compounds in various chemical systems using the Materials Project API.

Purpose: For graphitization catalyst research, find compounds with the highest carbon
content in specific chemical systems. Supports multiple search modes:
  - binary:    M-C       (metal carbides, e.g., LaC2, TiC)
  - ternary:   M-X-C     (metal + any partner + carbon, e.g., La-B-C, Ti-Si-C)
  - bimetal:   M1-M2-C   (two metals + carbon, e.g., Ta-Hf-C)
  - all:       all of the above combined

Usage:
    python screen_carbon_rich_compounds.py --api-key KEY --mode all --metal-group rare_earth
    python screen_carbon_rich_compounds.py --api-key KEY --mode binary --metal-group transition_3d
    python screen_carbon_rich_compounds.py --api-key KEY --mode ternary --metal-group rare_earth --partner-group nonmetal
    python screen_carbon_rich_compounds.py --api-key KEY --mode bimetal --metal-group rare_earth --partner-group transition_5d
    python screen_carbon_rich_compounds.py --api-key KEY --systems La-B-C Hf-Ta-C La-C

Requirements:
    pip install mp-api pymatgen pandas
"""

import argparse
import itertools
from dataclasses import dataclass, field
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

PARTNER_GROUPS = {
    # Non-metals
    "nonmetal": ["B", "N", "Si", "P", "S"],
    "nonmetal_extended": ["B", "N", "Si", "P", "S", "Se", "Te", "O", "F", "Cl"],
    # Metal groups (same as METAL_GROUPS, available as partner choices)
    "rare_earth": METAL_GROUPS["rare_earth"],
    "transition_3d": METAL_GROUPS["transition_3d"],
    "transition_4d": METAL_GROUPS["transition_4d"],
    "transition_5d": METAL_GROUPS["transition_5d"],
    "actinide": METAL_GROUPS["actinide"],
    "alkaline_earth": METAL_GROUPS["alkaline_earth"],
    "alkali": METAL_GROUPS["alkali"],
}

SEARCH_MODES = ["binary", "ternary", "bimetal", "all"]


# ── Screening funnel statistics ──────────────────────────────────────────────

@dataclass
class ScreeningStats:
    """Track compound counts at each stage of the screening funnel.

    All counts are from the SAME query (with API-level filters applied).
    No extra API calls are made — we just count as we apply each client-side filter.
    """
    systems_queried: int = 0
    total_returned: int = 0           # C-containing compounds returned from MP query
    experimental: int = 0             # subset that are experimentally synthesized
    theoretical: int = 0              # subset that are theoretical/predicted
    stable_on_hull: int = 0           # e_above_hull == 0
    pass_carbon_filter: int = 0       # meet min carbon fraction
    highest_c_per_system: int = 0     # best compound per system
    # These are filled later by run_screening.py after melting point annotation
    mp_in_range: int = 0              # melting point in target window
    mp_source_counts: dict = field(default_factory=dict)

    def print_funnel(self, min_c_frac: float, max_ehull: float,
                     mp_min: float = 0, mp_max: float = 0,
                     experimental_only: bool = False):
        """Print a visual funnel summary."""
        print("\n" + "=" * 70)
        print("SCREENING FUNNEL STATISTICS")
        print("=" * 70)
        print(f"  Chemical systems queried:             {self.systems_queried}")
        print(f"  C-containing compounds from MP:       {self.total_returned}")
        print(f"    - Experimentally synthesized:        {self.experimental}")
        print(f"    - Theoretical/predicted:             {self.theoretical}")
        print(f"  Stable on convex hull (e_hull = 0):   {self.stable_on_hull}")
        print(f"  Pass C fraction filter (>= {min_c_frac}):    {self.pass_carbon_filter}")
        print(f"  Highest-C compound per system:        {self.highest_c_per_system}")
        if mp_min > 0 or mp_max > 0:
            print(f"  Melting point in {mp_min}-{mp_max} C:       {self.mp_in_range}")
        if self.mp_source_counts:
            print(f"  Melting point sources:")
            for src, count in sorted(self.mp_source_counts.items()):
                print(f"    - {src}: {count}")
        if experimental_only:
            print(f"\n  Note: --experimental-only was ON, theoretical compounds excluded at API level")
        print(f"  Note: stability filter e_hull <= {max_ehull} eV/atom applied at API level")

    def to_dict(self) -> dict:
        """Convert to dict for CSV/JSON export."""
        d = {
            "systems_queried": self.systems_queried,
            "total_C_compounds_from_MP": self.total_returned,
            "experimental_synthesized": self.experimental,
            "theoretical_predicted": self.theoretical,
            "stable_on_hull": self.stable_on_hull,
            "pass_carbon_filter": self.pass_carbon_filter,
            "highest_C_per_system": self.highest_c_per_system,
            "mp_in_range": self.mp_in_range,
        }
        for src, count in self.mp_source_counts.items():
            d[f"mp_source_{src}"] = count
        return d


def generate_systems(
    mode: str,
    metal_group: Optional[str] = None,
    partner_group: Optional[str] = None,
    partner_elements: Optional[list[str]] = None,
) -> list[str]:
    """
    Generate chemical system strings based on search mode.

    Args:
        mode: One of "binary", "ternary", "bimetal", "all"
        metal_group: Key into METAL_GROUPS for the primary metals
        partner_group: Key into PARTNER_GROUPS for the partner elements (ternary/bimetal)
        partner_elements: Explicit list of partner elements (overrides partner_group)

    Returns:
        Deduplicated list of chemical system strings (e.g., ["La-C", "La-B-C", ...])
    """
    if metal_group is None:
        raise ValueError("--metal-group is required when not using --systems")

    metals = METAL_GROUPS.get(metal_group)
    if metals is None:
        raise ValueError(f"Unknown metal group: {metal_group}. Options: {list(METAL_GROUPS.keys())}")

    # Resolve partner elements
    partners = None
    if partner_elements:
        partners = partner_elements
    elif partner_group:
        partners = PARTNER_GROUPS.get(partner_group)
        if partners is None:
            raise ValueError(f"Unknown partner group: {partner_group}. Options: {list(PARTNER_GROUPS.keys())}")

    systems = set()

    # Binary: M-C
    if mode in ("binary", "all"):
        for m in metals:
            systems.add(f"{m}-C")

    # Ternary: M-partner-C (partner can be anything)
    if mode in ("ternary", "all"):
        if partners is None:
            raise ValueError("--partner-group or --partner-elements required for ternary/all mode")
        for m, x in itertools.product(metals, partners):
            if m != x:
                elements = sorted([m, x, "C"])
                systems.add("-".join(elements))

    # Bimetal: M1-M2-C (two metals from potentially different groups)
    if mode in ("bimetal", "all"):
        if partners is None:
            for m1, m2 in itertools.combinations(metals, 2):
                elements = sorted([m1, m2, "C"])
                systems.add("-".join(elements))
        else:
            for m, p in itertools.product(metals, partners):
                if m != p:
                    elements = sorted([m, p, "C"])
                    systems.add("-".join(elements))

    return sorted(systems)


def query_chemsys(
    mpr: MPRester,
    chemsys: str,
    max_energy_above_hull: float = 0.1,
    experimental_only: bool = False,
) -> list[dict]:
    """
    Query compounds in a chemical system from Materials Project.

    Applies stability and experimental filters at the API level to minimize data transfer.
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

        # Must contain carbon
        if c_frac == 0:
            continue

        # Must contain at least 2 elements (not pure C)
        n_elements = len(comp.elements)
        if n_elements < 2:
            continue

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
            "energy_above_hull_eV": round(doc.energy_above_hull, 4) if doc.energy_above_hull is not None else None,
            "formation_energy_eV": round(doc.formation_energy_per_atom, 4) if doc.formation_energy_per_atom is not None else None,
            "spacegroup": spacegroup,
            "crystal_system": crystal_system,
            "density_g_cm3": round(doc.density, 2) if doc.density else None,
            "theoretical": getattr(doc, "theoretical", None),
        })

    return results


def screen_systems(
    api_key: str,
    systems: Optional[list[str]] = None,
    mode: str = "ternary",
    metal_group: Optional[str] = None,
    partner_group: Optional[str] = None,
    partner_elements: Optional[list[str]] = None,
    min_carbon_fraction: float = 0.2,
    max_energy_above_hull: float = 0.1,
    experimental_only: bool = False,
) -> tuple[pd.DataFrame, ScreeningStats]:
    """
    Screen chemical systems for carbon-rich compounds.

    Uses the same filtered API query as before (no extra calls).
    Counts are tracked at each client-side filtering step.

    Returns:
        (filtered_df, stats) tuple
    """
    if systems is None:
        systems = generate_systems(mode, metal_group, partner_group, partner_elements)

    stats = ScreeningStats()
    stats.systems_queried = len(systems)
    print(f"Will query {len(systems)} chemical systems")

    # Collect results with API-level filters, deduplicate by material_id
    seen_ids = set()
    all_results = []
    with MPRester(api_key) as mpr:
        for i, sys in enumerate(systems):
            print(f"[{i+1}/{len(systems)}] Querying {sys} ...")
            results = query_chemsys(mpr, sys, max_energy_above_hull, experimental_only)
            new_count = 0
            for r in results:
                if r["material_id"] not in seen_ids:
                    seen_ids.add(r["material_id"])
                    all_results.append(r)
                    new_count += 1
            if new_count > 0:
                print(f"  Found {len(results)} compounds ({new_count} new)")
            else:
                print(f"  No new compounds found")

    if not all_results:
        return pd.DataFrame(), stats

    df = pd.DataFrame(all_results)

    # ── Count at each step (no extra queries) ────────────────────────────
    # Step 1: total returned from MP (after API-level stability + experimental filters)
    stats.total_returned = len(df)
    stats.experimental = int((df["theoretical"] == False).sum())  # noqa: E712
    stats.theoretical = int((df["theoretical"] == True).sum())  # noqa: E712
    stats.stable_on_hull = int((df["energy_above_hull_eV"] == 0).sum())

    # Step 2: apply carbon fraction filter (client-side)
    df = df[df["C_atomic_fraction"] >= min_carbon_fraction].reset_index(drop=True)
    stats.pass_carbon_filter = len(df)

    # Step 3: count unique systems (= highest-C per system count)
    if not df.empty:
        stats.highest_c_per_system = df["chemsys"].nunique()
        df = df.sort_values("C_atomic_fraction", ascending=False).reset_index(drop=True)

    return df, stats


def find_highest_carbon_per_system(df: pd.DataFrame) -> pd.DataFrame:
    """For each chemsys, return only the compound with highest C atomic fraction."""
    if df.empty:
        return df
    idx = df.groupby("chemsys")["C_atomic_fraction"].idxmax()
    return df.loc[idx].sort_values("C_atomic_fraction", ascending=False).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(
        description="Screen chemical systems for carbon-rich compounds via Materials Project API"
    )
    parser.add_argument("--api-key", required=True, help="Materials Project API key")
    parser.add_argument(
        "--systems", nargs="+", default=None,
        help="Explicit chemical systems to query (e.g., La-B-C Hf-Ta-C La-C). Bypasses --mode"
    )
    parser.add_argument(
        "--mode", default="ternary", choices=SEARCH_MODES,
        help=(
            "Search mode: "
            "'binary' = M-C only; "
            "'ternary' = M-partner-C; "
            "'bimetal' = M1-M2-C; "
            "'all' = binary + ternary + bimetal combined. "
            "(default: ternary)"
        ),
    )
    parser.add_argument(
        "--metal-group", default=None,
        choices=list(METAL_GROUPS.keys()),
        help="Primary metal element group"
    )
    parser.add_argument(
        "--partner-group", default=None,
        choices=list(PARTNER_GROUPS.keys()),
        help=(
            "Partner element group for ternary/bimetal modes. "
            "Can be a non-metal group (e.g., 'nonmetal', 'nonmetal_extended') "
            "or a metal group (e.g., 'transition_5d') for M1-M2-C screening"
        ),
    )
    parser.add_argument(
        "--partner-elements", nargs="+", default=None,
        help="Explicit partner elements (overrides --partner-group). E.g., --partner-elements B N Si Hf Ta"
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
    print(f"Mode: {args.mode}" + (" (bypassed by --systems)" if args.systems else ""))
    print("=" * 70)

    df, stats = screen_systems(
        api_key=args.api_key,
        systems=args.systems,
        mode=args.mode,
        metal_group=args.metal_group,
        partner_group=args.partner_group,
        partner_elements=args.partner_elements,
        min_carbon_fraction=args.min_c_fraction,
        max_energy_above_hull=args.max_ehull,
        experimental_only=args.experimental_only,
    )

    if df.empty:
        print("\nNo compounds found matching criteria.")
        stats.print_funnel(args.min_c_fraction, args.max_ehull,
                           experimental_only=args.experimental_only)
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

    # Print funnel
    stats.print_funnel(args.min_c_fraction, args.max_ehull,
                       experimental_only=args.experimental_only)

    # Summary
    print(f"\n--- Summary ---")
    print(f"Total compounds found: {len(df)}")
    print(f"Unique systems with hits: {df['chemsys'].nunique()}")
    print(f"Highest C atomic fraction: {df['C_atomic_fraction'].max():.4f} "
          f"({df.loc[df['C_atomic_fraction'].idxmax(), 'formula']})")


if __name__ == "__main__":
    main()
