"""
Screen carbon-rich compounds from internal MongoDB mirror of Materials Project.

This module queries the internal opendb.mp_2022 database (154718 entries) instead
of the Materials Project API, avoiding network restrictions on company machines.

The database schema matches MP 2022 with fields:
    material_id, formula_pretty, chemsys, elements, nelements, composition,
    composition_reduced, energy_above_hull, formation_energy_per_atom,
    symmetry, structure, volume, density, nsites, etc.

Screening logic:
    - Include compounds that are EITHER experimentally synthesized OR within
      the e_above_hull threshold (or both)
    - Rank purely by carbon atomic fraction (highest first)
    - Experimental compounds are trusted to exist regardless of e_above_hull

Usage:
    python screen_carbon_rich_compounds_internal.py --mode all --metal-group rare_earth --partner-group nonmetal
    python screen_carbon_rich_compounds_internal.py --mode binary --metal-group transition_3d
    python screen_carbon_rich_compounds_internal.py --systems La-B-C Hf-Ta-C La-C

Requirements:
    pip install pymongo pymatgen pandas
"""

import argparse
import itertools
from dataclasses import dataclass, field
from typing import Optional

import pandas as pd
from pymongo import MongoClient
from pymatgen.core import Composition


# ── MongoDB connection ───────────────────────────────────────────────────────

MONGO_URI = "mongodb://QiuGT:woshiQiuGT6%40@10.156.204.60:27017/?readPreference=primary&appname=MongoDB%20Compass&ssl=false"
MONGO_DB = "opendb"
MONGO_COLLECTION = "mp_2022"


def get_collection(mongo_uri: str = MONGO_URI, db_name: str = MONGO_DB,
                   collection_name: str = MONGO_COLLECTION):
    """Connect to the internal MongoDB and return the collection handle."""
    client = MongoClient(mongo_uri)
    db = client[db_name]
    return db[collection_name]


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
    "nonmetal": ["B", "N", "Si", "P", "S"],
    "nonmetal_extended": ["B", "N", "Si", "P", "S", "Se", "Te", "O", "F", "Cl"],
    "rare_earth": METAL_GROUPS["rare_earth"],
    "transition_3d": METAL_GROUPS["transition_3d"],
    "transition_4d": METAL_GROUPS["transition_4d"],
    "transition_5d": METAL_GROUPS["transition_5d"],
    "actinide": METAL_GROUPS["actinide"],
    "alkaline_earth": METAL_GROUPS["alkaline_earth"],
    "alkali": METAL_GROUPS["alkali"],
}

SEARCH_MODES = ["binary", "ternary", "bimetal", "all", "combinations", "comprehensive"]


# ── Screening funnel statistics ──────────────────────────────────────────────

@dataclass
class ScreeningStats:
    """Track compound counts at each stage of the screening funnel."""
    systems_queried: int = 0
    total_returned: int = 0
    experimental: int = 0
    not_experimental: int = 0
    stable_on_hull: int = 0
    pass_stability_or_experimental: int = 0
    pass_carbon_filter: int = 0
    highest_c_per_system: int = 0
    mp_in_range: int = 0
    mp_source_counts: dict = field(default_factory=dict)

    def print_funnel(self, min_c_frac: float, max_ehull: float,
                     mp_min: float = 0, mp_max: float = 0):
        """Print a visual funnel summary."""
        print("\n" + "=" * 70)
        print("SCREENING FUNNEL STATISTICS")
        print("=" * 70)
        if self.systems_queried > 0:
            print(f"  Chemical systems queried:                 {self.systems_queried}")
        else:
            print(f"  Mode: comprehensive (all carbides in DB)")
        print(f"  C-containing compounds from DB:           {self.total_returned}")
        print(f"    - Experimental:                          {self.experimental}")
        print(f"    - Not experimental:                      {self.not_experimental}")
        print(f"  Stable on convex hull (e_hull = 0):       {self.stable_on_hull}")
        print(f"  Pass filter (experimental OR e_hull<={max_ehull}): {self.pass_stability_or_experimental}")
        print(f"  Pass C fraction filter (>= {min_c_frac}):        {self.pass_carbon_filter}")
        print(f"  Highest-C compound per system:            {self.highest_c_per_system}")
        if mp_min > 0 or mp_max > 0:
            print(f"  Melting point in {mp_min}-{mp_max} C:           {self.mp_in_range}")
        if self.mp_source_counts:
            print(f"  Melting point sources:")
            for src, count in sorted(self.mp_source_counts.items()):
                print(f"    - {src}: {count}")

    def to_dict(self) -> dict:
        d = {
            "systems_queried": self.systems_queried,
            "total_C_compounds_from_DB": self.total_returned,
            "experimental": self.experimental,
            "not_experimental": self.not_experimental,
            "stable_on_hull": self.stable_on_hull,
            "pass_stability_or_experimental": self.pass_stability_or_experimental,
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
    """Generate chemical system strings based on search mode.

    Modes:
        binary:        M-C for one metal group
        ternary:       M-partner-C for one metal group + one partner group
        bimetal:       M1-M2-C for one metal group + one partner group
        all:           binary + ternary + bimetal for one metal group + partner
        combinations:  ALL metal groups x ALL partner groups, including:
                       - binary M-C for every metal group
                       - ternary M-X-C for every metal group x every nonmetal group
                       - bimetal M1-M2-C for every metal group x every metal group
                       No --metal-group or --partner-group needed.
        comprehensive: queries entire DB directly, no system generation
    """
    if mode == "comprehensive":
        return []

    if mode == "combinations":
        systems = set()
        all_metals = []
        for group_name, elements in METAL_GROUPS.items():
            all_metals.extend(elements)

        # Deduplicate metals
        all_metals = list(dict.fromkeys(all_metals))

        # Binary: every metal + C
        for m in all_metals:
            systems.add("-".join(sorted([m, "C"])))

        # Ternary: every metal x every nonmetal (including extended) + C
        all_nonmetals = list(dict.fromkeys(
            PARTNER_GROUPS["nonmetal_extended"]
        ))
        for m, x in itertools.product(all_metals, all_nonmetals):
            if m != x:
                systems.add("-".join(sorted([m, x, "C"])))

        # Bimetal: every metal x every metal + C (unique pairs)
        for m1, m2 in itertools.combinations(all_metals, 2):
            systems.add("-".join(sorted([m1, m2, "C"])))

        return sorted(systems)

    if metal_group is None:
        raise ValueError("--metal-group is required when not using --systems")

    metals = METAL_GROUPS.get(metal_group)
    if metals is None:
        raise ValueError(f"Unknown metal group: {metal_group}. Options: {list(METAL_GROUPS.keys())}")

    partners = None
    if partner_elements:
        partners = partner_elements
    elif partner_group:
        partners = PARTNER_GROUPS.get(partner_group)
        if partners is None:
            raise ValueError(f"Unknown partner group: {partner_group}. Options: {list(PARTNER_GROUPS.keys())}")

    systems = set()

    if mode in ("binary", "all"):
        for m in metals:
            elements = sorted([m, "C"])
            systems.add("-".join(elements))

    if mode in ("ternary", "all"):
        if partners is None:
            raise ValueError("--partner-group or --partner-elements required for ternary/all mode")
        for m, x in itertools.product(metals, partners):
            if m != x:
                elements = sorted([m, x, "C"])
                systems.add("-".join(elements))

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


def build_output_name(mode: str, metal_group: Optional[str],
                      partner_group: Optional[str],
                      partner_elements: Optional[list[str]],
                      systems: Optional[list[str]]) -> str:
    """Build a descriptive output filename from the search parameters."""
    if mode == "comprehensive":
        return "carbon_rich_comprehensive.csv"

    if mode == "combinations":
        return "carbon_rich_combinations.csv"

    if systems:
        sys_str = "_".join(systems[:3])
        if len(systems) > 3:
            sys_str += f"_etc{len(systems)}"
        return f"carbon_rich_{sys_str}.csv"

    parts = [mode]
    if metal_group:
        parts.append(metal_group)
    if partner_group:
        parts.append(partner_group)
    elif partner_elements:
        parts.append("-".join(partner_elements[:3]))

    return f"carbon_rich_{'_'.join(parts)}.csv"


def _parse_doc(doc, chemsys_override: str = "") -> Optional[dict]:
    """Parse a MongoDB document into a result dict. Returns None if invalid."""
    comp_dict = doc.get("composition_reduced", {})
    if not comp_dict or "C" not in comp_dict:
        return None

    comp = Composition(comp_dict)
    c_frac = comp.get_atomic_fraction("C")
    if c_frac == 0:
        return None

    n_elements = doc.get("nelements", len(comp.elements))
    if n_elements < 2:
        return None

    c_wt_frac = comp.get_wt_fraction("C")

    spacegroup = ""
    crystal_system = ""
    sym = doc.get("symmetry", {})
    if sym:
        spacegroup = sym.get("symbol", "")
        crystal_system = sym.get("crystal_system", "")

    ehull = doc.get("energy_above_hull")
    fe = doc.get("formation_energy_per_atom")

    is_theoretical = doc.get("theoretical", None)
    if is_theoretical is not None:
        is_experimental = not is_theoretical
    else:
        is_experimental = None

    return {
        "material_id": doc.get("material_id", str(doc.get("_id", ""))),
        "formula": doc.get("formula_pretty", str(comp.reduced_formula)),
        "chemsys": doc.get("chemsys", chemsys_override),
        "C_atomic_fraction": round(c_frac, 4),
        "C_weight_fraction": round(c_wt_frac, 4),
        "n_elements": n_elements,
        "energy_above_hull_eV": round(ehull, 4) if ehull is not None else None,
        "formation_energy_eV": round(fe, 4) if fe is not None else None,
        "spacegroup": spacegroup,
        "crystal_system": crystal_system,
        "density_g_cm3": round(doc.get("density", 0), 2) if doc.get("density") else None,
        "experimental": is_experimental,
    }


def query_all_carbides(col) -> list[dict]:
    """
    Query ALL carbon-containing compounds from the entire database in one go.

    Uses elements array to find any compound containing C with at least 2 elements.
    This is the comprehensive mode — no system-by-system enumeration needed.
    """
    query = {
        "elements": {"$all": ["C"]},
        "nelements": {"$gte": 2},
    }

    projection = {
        "material_id": 1, "formula_pretty": 1, "chemsys": 1,
        "composition_reduced": 1, "nelements": 1, "nsites": 1,
        "energy_above_hull": 1, "formation_energy_per_atom": 1,
        "symmetry": 1, "volume": 1, "density": 1, "theoretical": 1,
    }

    try:
        cursor = col.find(query, projection)
        docs = list(cursor)
    except Exception as e:
        print(f"  Warning: failed to query all carbides: {e}")
        return []

    results = []
    for doc in docs:
        r = _parse_doc(doc)
        if r is not None:
            results.append(r)

    return results


def query_chemsys_mongo(
    col,
    chemsys: str,
) -> list[dict]:
    """
    Query ALL C-containing compounds in a chemical system from the internal MongoDB.

    No stability filter at the DB level — we fetch everything and filter client-side
    so we can count the funnel properly. The chemsys field in mp_2022 uses
    alphabetically sorted elements joined by "-".
    """
    target_elements = chemsys.split("-")
    if "C" not in target_elements:
        raise ValueError(f"System {chemsys} does not contain carbon")

    # Build all sub-chemsys that contain C
    non_c_elements = [e for e in target_elements if e != "C"]
    sub_systems = set()
    for r in range(len(non_c_elements) + 1):
        for combo in itertools.combinations(non_c_elements, r):
            sub_elements = sorted(list(combo) + ["C"])
            sub_systems.add("-".join(sub_elements))
    sub_systems.discard("C")

    query = {"chemsys": {"$in": list(sub_systems)}}

    projection = {
        "material_id": 1, "formula_pretty": 1, "chemsys": 1,
        "composition_reduced": 1, "nelements": 1, "nsites": 1,
        "energy_above_hull": 1, "formation_energy_per_atom": 1,
        "symmetry": 1, "volume": 1, "density": 1, "theoretical": 1,
    }

    try:
        cursor = col.find(query, projection)
        docs = list(cursor)
    except Exception as e:
        print(f"  Warning: failed to query {chemsys}: {e}")
        return []

    results = []
    for doc in docs:
        r = _parse_doc(doc, chemsys)
        if r is not None:
            results.append(r)

    return results


def screen_systems(
    systems: Optional[list[str]] = None,
    mode: str = "ternary",
    metal_group: Optional[str] = None,
    partner_group: Optional[str] = None,
    partner_elements: Optional[list[str]] = None,
    min_carbon_fraction: float = 0.2,
    max_energy_above_hull: float = 0.1,
    mongo_uri: str = MONGO_URI,
    db_name: str = MONGO_DB,
    collection_name: str = MONGO_COLLECTION,
) -> tuple[pd.DataFrame, ScreeningStats]:
    """
    Screen chemical systems for carbon-rich compounds from internal MongoDB.

    A compound passes the filter if it is:
    - Experimentally synthesized (regardless of e_above_hull), OR
    - Within the e_above_hull threshold

    Then ranked purely by carbon atomic fraction (highest first).

    Returns:
        (filtered_df, stats) tuple
    """
    col = get_collection(mongo_uri, db_name, collection_name)
    stats = ScreeningStats()

    # ── Comprehensive mode: single query for ALL carbides ────────────────
    if mode == "comprehensive":
        print("Comprehensive mode: querying ALL carbon-containing compounds in database...")
        all_results = query_all_carbides(col)
        stats.systems_queried = 0  # not system-based
        print(f"  Found {len(all_results)} C-containing compounds total")

    # ── System-by-system mode ────────────────────────────────────────────
    else:
        if systems is None:
            systems = generate_systems(mode, metal_group, partner_group, partner_elements)

        stats.systems_queried = len(systems)
        print(f"Will query {len(systems)} chemical systems from internal DB")

        seen_ids = set()
        all_results = []
        for i, sys in enumerate(systems):
            print(f"[{i+1}/{len(systems)}] Querying {sys} ...")
            results = query_chemsys_mongo(col, sys)
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

    # ── Count funnel stages ──────────────────────────────────────────────
    stats.total_returned = len(df)
    if "experimental" in df.columns:
        stats.experimental = int((df["experimental"] == True).sum())  # noqa: E712
        stats.not_experimental = int((df["experimental"] == False).sum())  # noqa: E712
    stats.stable_on_hull = int((df["energy_above_hull_eV"] == 0).sum())

    # ── Filter: experimental OR within e_hull threshold ──────────────────
    mask_experimental = df["experimental"] == True  # noqa: E712
    mask_stable = df["energy_above_hull_eV"] <= max_energy_above_hull
    df = df[mask_experimental | mask_stable].reset_index(drop=True)
    stats.pass_stability_or_experimental = len(df)

    # ── Filter: carbon fraction ──────────────────────────────────────────
    df = df[df["C_atomic_fraction"] >= min_carbon_fraction].reset_index(drop=True)
    stats.pass_carbon_filter = len(df)

    # ── Rank by C content (highest first) ────────────────────────────────
    if not df.empty:
        df = df.sort_values("C_atomic_fraction", ascending=False).reset_index(drop=True)
        stats.highest_c_per_system = df["chemsys"].nunique()

    return df, stats


def find_highest_carbon_per_system(df: pd.DataFrame) -> pd.DataFrame:
    """For each chemsys, return the compound with the highest C atomic fraction."""
    if df.empty:
        return df
    idx = df.groupby("chemsys")["C_atomic_fraction"].idxmax()
    return df.loc[idx].sort_values("C_atomic_fraction", ascending=False).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(
        description="Screen chemical systems for carbon-rich compounds from internal MongoDB (opendb.mp_2022)"
    )
    parser.add_argument(
        "--systems", nargs="+", default=None,
        help="Explicit chemical systems to query (e.g., B-C-La C-Hf-Ta). Bypasses --mode"
    )
    parser.add_argument(
        "--mode", default="ternary", choices=SEARCH_MODES,
        help=("Search mode: 'binary'=M-C, 'ternary'=M-partner-C, 'bimetal'=M1-M2-C, "
              "'all'=combined, 'combinations'=all metal/partner group combos, "
              "'comprehensive'=ALL carbides in DB (default: ternary)")
    )
    parser.add_argument(
        "--metal-group", default=None, choices=list(METAL_GROUPS.keys()),
        help="Primary metal element group"
    )
    parser.add_argument(
        "--partner-group", default=None, choices=list(PARTNER_GROUPS.keys()),
        help="Partner element group for ternary/bimetal modes"
    )
    parser.add_argument(
        "--partner-elements", nargs="+", default=None,
        help="Explicit partner elements (overrides --partner-group)"
    )
    parser.add_argument("--min-c-fraction", type=float, default=0.2,
                        help="Minimum carbon atomic fraction (default: 0.2)")
    parser.add_argument("--max-ehull", type=float, default=0.1,
                        help="Maximum energy above hull in eV/atom (default: 0.1)")
    parser.add_argument("--mongo-uri", default=MONGO_URI,
                        help="MongoDB connection URI")
    parser.add_argument("--db-name", default=MONGO_DB,
                        help=f"Database name (default: {MONGO_DB})")
    parser.add_argument("--collection", default=MONGO_COLLECTION,
                        help=f"Collection name (default: {MONGO_COLLECTION})")
    parser.add_argument("--output", default=None,
                        help="Output CSV filename (auto-generated if not specified)")
    args = parser.parse_args()

    # Auto-generate output filename from search parameters
    output = args.output or build_output_name(
        args.mode, args.metal_group, args.partner_group,
        args.partner_elements, args.systems)

    print("=" * 70)
    print("Carbon-Rich Compound Screening (Internal MongoDB)")
    print(f"Mode: {args.mode}" + (" (bypassed by --systems)" if args.systems else ""))
    print(f"Database: {args.db_name}.{args.collection}")
    print("=" * 70)

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
        print("\nNo compounds found matching criteria.")
        stats.print_funnel(args.min_c_fraction, args.max_ehull)
        return

    # Save full results
    df.to_csv(output, index=False)
    print(f"\nFull results saved to {output} ({len(df)} compounds)")

    # Highest-carbon per system
    best = find_highest_carbon_per_system(df)
    best_file = output.replace(".csv", "_best_per_system.csv")
    best.to_csv(best_file, index=False)
    print(f"Best per system saved to {best_file} ({len(best)} systems)")

    print("\n" + "=" * 70)
    print("Top carbon-rich compounds (highest C fraction per system):")
    print("=" * 70)
    display_cols = ["formula", "chemsys", "C_atomic_fraction", "C_weight_fraction",
                    "energy_above_hull_eV", "experimental", "material_id"]
    print(best[display_cols].to_string(index=False))

    stats.print_funnel(args.min_c_fraction, args.max_ehull)


if __name__ == "__main__":
    main()
