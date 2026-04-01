"""
Screen carbon-rich compounds from internal MongoDB mirror of Materials Project.

This module queries the internal opendb.mp_2022 database (154718 entries) instead
of the Materials Project API, avoiding network restrictions on company machines.

The database schema matches MP 2022 with fields:
    material_id, formula_pretty, chemsys, elements, nelements, composition,
    composition_reduced, energy_above_hull, formation_energy_per_atom,
    symmetry, structure, volume, density, nsites, etc.

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

SEARCH_MODES = ["binary", "ternary", "bimetal", "all"]


# ── Screening funnel statistics ──────────────────────────────────────────────

@dataclass
class ScreeningStats:
    """Track compound counts at each stage of the screening funnel."""
    systems_queried: int = 0
    total_returned: int = 0
    experimental: int = 0
    theoretical: int = 0
    stable_on_hull: int = 0
    pass_carbon_filter: int = 0
    highest_c_per_system: int = 0
    mp_in_range: int = 0
    mp_source_counts: dict = field(default_factory=dict)

    def print_funnel(self, min_c_frac: float, max_ehull: float,
                     mp_min: float = 0, mp_max: float = 0,
                     experimental_only: bool = False):
        """Print a visual funnel summary."""
        print("\n" + "=" * 70)
        print("SCREENING FUNNEL STATISTICS")
        print("=" * 70)
        print(f"  Chemical systems queried:             {self.systems_queried}")
        print(f"  C-containing compounds from DB:       {self.total_returned}")
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
            print(f"\n  Note: --experimental-only was ON, theoretical compounds excluded")
        print(f"  Note: stability filter e_hull <= {max_ehull} eV/atom applied")

    def to_dict(self) -> dict:
        d = {
            "systems_queried": self.systems_queried,
            "total_C_compounds_from_DB": self.total_returned,
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
    """Generate chemical system strings based on search mode."""
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
            # chemsys in DB is alphabetically sorted, e.g. "C-La" not "La-C"
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


def query_chemsys_mongo(
    col,
    chemsys: str,
    max_energy_above_hull: float = 0.1,
    experimental_only: bool = False,
) -> list[dict]:
    """
    Query compounds in a chemical system from the internal MongoDB.

    The chemsys field in mp_2022 uses alphabetically sorted elements joined by "-",
    and querying a chemsys returns exact matches AND all sub-systems.
    We query for the exact chemsys, then also query sub-systems containing C.
    """
    target_elements = chemsys.split("-")
    if "C" not in target_elements:
        raise ValueError(f"System {chemsys} does not contain carbon")

    # Build all sub-chemsys that contain C
    # e.g., for "B-C-La", sub-systems with C are: "C", "B-C", "C-La", "B-C-La"
    non_c_elements = [e for e in target_elements if e != "C"]
    sub_systems = set()
    for r in range(len(non_c_elements) + 1):
        for combo in itertools.combinations(non_c_elements, r):
            sub_elements = sorted(list(combo) + ["C"])
            sub_systems.add("-".join(sub_elements))
    # Remove pure "C" — we want at least one other element
    sub_systems.discard("C")

    # MongoDB query
    query = {
        "chemsys": {"$in": list(sub_systems)},
        "energy_above_hull": {"$lte": max_energy_above_hull},
    }

    # The mp_2022 collection may not have a 'theoretical' field.
    # MP marks entries as theoretical if they lack ICSD IDs.
    # We check if the field exists before filtering.
    if experimental_only:
        query["theoretical"] = False

    projection = {
        "material_id": 1,
        "formula_pretty": 1,
        "chemsys": 1,
        "composition_reduced": 1,
        "nelements": 1,
        "nsites": 1,
        "energy_above_hull": 1,
        "formation_energy_per_atom": 1,
        "symmetry": 1,
        "volume": 1,
        "density": 1,
        "theoretical": 1,
    }

    try:
        cursor = col.find(query, projection)
        docs = list(cursor)
    except Exception as e:
        print(f"  Warning: failed to query {chemsys}: {e}")
        return []

    results = []
    for doc in docs:
        # Compute carbon fraction from composition_reduced
        comp_dict = doc.get("composition_reduced", {})
        if not comp_dict or "C" not in comp_dict:
            continue

        comp = Composition(comp_dict)
        c_frac = comp.get_atomic_fraction("C")
        if c_frac == 0:
            continue

        n_elements = doc.get("nelements", len(comp.elements))
        if n_elements < 2:
            continue

        c_wt_frac = comp.get_wt_fraction("C")

        spacegroup = ""
        crystal_system = ""
        sym = doc.get("symmetry", {})
        if sym:
            spacegroup = sym.get("symbol", "")
            crystal_system = sym.get("crystal_system", "")

        ehull = doc.get("energy_above_hull")
        fe = doc.get("formation_energy_per_atom")

        results.append({
            "material_id": doc.get("material_id", str(doc.get("_id", ""))),
            "formula": doc.get("formula_pretty", str(comp.reduced_formula)),
            "chemsys": doc.get("chemsys", chemsys),
            "C_atomic_fraction": round(c_frac, 4),
            "C_weight_fraction": round(c_wt_frac, 4),
            "n_elements": n_elements,
            "energy_above_hull_eV": round(ehull, 4) if ehull is not None else None,
            "formation_energy_eV": round(fe, 4) if fe is not None else None,
            "spacegroup": spacegroup,
            "crystal_system": crystal_system,
            "density_g_cm3": round(doc.get("density", 0), 2) if doc.get("density") else None,
            "theoretical": doc.get("theoretical", None),
        })

    return results


def screen_systems(
    systems: Optional[list[str]] = None,
    mode: str = "ternary",
    metal_group: Optional[str] = None,
    partner_group: Optional[str] = None,
    partner_elements: Optional[list[str]] = None,
    min_carbon_fraction: float = 0.2,
    max_energy_above_hull: float = 0.1,
    experimental_only: bool = False,
    mongo_uri: str = MONGO_URI,
    db_name: str = MONGO_DB,
    collection_name: str = MONGO_COLLECTION,
) -> tuple[pd.DataFrame, ScreeningStats]:
    """
    Screen chemical systems for carbon-rich compounds from internal MongoDB.

    Returns:
        (filtered_df, stats) tuple
    """
    if systems is None:
        systems = generate_systems(mode, metal_group, partner_group, partner_elements)

    stats = ScreeningStats()
    stats.systems_queried = len(systems)
    print(f"Will query {len(systems)} chemical systems from internal DB")

    col = get_collection(mongo_uri, db_name, collection_name)

    # Collect results, deduplicate by material_id
    seen_ids = set()
    all_results = []
    for i, sys in enumerate(systems):
        print(f"[{i+1}/{len(systems)}] Querying {sys} ...")
        results = query_chemsys_mongo(col, sys, max_energy_above_hull, experimental_only)
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

    # ── Count at each step ───────────────────────────────────────────────
    stats.total_returned = len(df)
    if "theoretical" in df.columns:
        stats.experimental = int((df["theoretical"] == False).sum())  # noqa: E712
        stats.theoretical = int((df["theoretical"] == True).sum())  # noqa: E712
    else:
        # If theoretical field doesn't exist, all are from MP so count as unknown
        stats.experimental = len(df)
        stats.theoretical = 0
    stats.stable_on_hull = int((df["energy_above_hull_eV"] == 0).sum())

    # Apply carbon fraction filter
    df = df[df["C_atomic_fraction"] >= min_carbon_fraction].reset_index(drop=True)
    stats.pass_carbon_filter = len(df)

    if not df.empty:
        stats.highest_c_per_system = df["chemsys"].nunique()
        df = df.sort_values("C_atomic_fraction", ascending=False).reset_index(drop=True)

    return df, stats


def find_highest_carbon_per_system(
    df: pd.DataFrame,
    weight_c: float = 0.7,
    weight_stability: float = 0.3,
) -> pd.DataFrame:
    """For each chemsys, return the best compound by C fraction + stability.

    Selection criteria (no melting point involved — that's a downstream step):
    - Carbon atomic fraction (higher = better), weighted by weight_c
    - Thermodynamic stability (lower e_above_hull = better), weighted by weight_stability

    Args:
        df: DataFrame with C_atomic_fraction and energy_above_hull_eV columns
        weight_c: Weight for carbon fraction score (default 0.7)
        weight_stability: Weight for stability score (default 0.3)
    """
    if df.empty:
        return df

    df = df.copy()

    # Normalize C fraction to [0, 1]
    c_max = df["C_atomic_fraction"].max()
    c_min = df["C_atomic_fraction"].min()
    if c_max > c_min:
        df["_score_c"] = (df["C_atomic_fraction"] - c_min) / (c_max - c_min)
    else:
        df["_score_c"] = 1.0

    # Normalize stability (invert: lower e_hull = higher score)
    if "energy_above_hull_eV" in df.columns:
        e_max = df["energy_above_hull_eV"].max()
        if e_max > 0:
            df["_score_stab"] = 1.0 - df["energy_above_hull_eV"] / e_max
        else:
            df["_score_stab"] = 1.0
    else:
        df["_score_stab"] = 0.5

    df["_selection_score"] = weight_c * df["_score_c"] + weight_stability * df["_score_stab"]

    # Pick the best per system
    idx = df.groupby("chemsys")["_selection_score"].idxmax()
    result = df.loc[idx].sort_values("_selection_score", ascending=False).reset_index(drop=True)

    # Clean up temp columns
    result = result.drop(columns=["_score_c", "_score_stab", "_selection_score"])
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Screen chemical systems for carbon-rich compounds from internal MongoDB (opendb.mp_2022)"
    )
    parser.add_argument(
        "--systems", nargs="+", default=None,
        help="Explicit chemical systems to query (e.g., La-B-C Hf-Ta-C La-C). Bypasses --mode"
    )
    parser.add_argument(
        "--mode", default="ternary", choices=SEARCH_MODES,
        help="Search mode: 'binary'=M-C, 'ternary'=M-partner-C, 'bimetal'=M1-M2-C, 'all'=combined (default: ternary)"
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
    parser.add_argument("--experimental-only", action="store_true",
                        help="Only include experimentally synthesized compounds")
    parser.add_argument("--mongo-uri", default=MONGO_URI,
                        help=f"MongoDB connection URI (default: {MONGO_URI})")
    parser.add_argument("--db-name", default=MONGO_DB,
                        help=f"Database name (default: {MONGO_DB})")
    parser.add_argument("--collection", default=MONGO_COLLECTION,
                        help=f"Collection name (default: {MONGO_COLLECTION})")
    parser.add_argument("--output", default="carbon_rich_compounds.csv",
                        help="Output CSV filename")
    args = parser.parse_args()

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
        experimental_only=args.experimental_only,
        mongo_uri=args.mongo_uri,
        db_name=args.db_name,
        collection_name=args.collection,
    )

    if df.empty:
        print("\nNo compounds found matching criteria.")
        stats.print_funnel(args.min_c_fraction, args.max_ehull,
                           experimental_only=args.experimental_only)
        return

    # Save full results
    df.to_csv(args.output, index=False)
    print(f"\nFull results saved to {args.output} ({len(df)} compounds)")

    # Highest-carbon per system
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

    stats.print_funnel(args.min_c_fraction, args.max_ehull,
                       experimental_only=args.experimental_only)


if __name__ == "__main__":
    main()
