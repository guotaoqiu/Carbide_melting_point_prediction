# Carbide Melting Point Prediction & Carbon-Rich Compound Screening

High-throughput screening toolkit for identifying carbon-rich compounds in ternary M-X-C chemical systems as candidates for graphitization catalysis research.

## Installation

```bash
pip install -r requirements.txt
```

**Requirements:**
- Python >= 3.10
- `mp-api` >= 0.41.0 (Materials Project API client)
- `pymatgen` >= 2024.1.1 (Python Materials Genomics)
- `pandas` >= 2.0

You also need a **Materials Project API key**. Get one free at: https://materialsproject.org/api

## Files Overview

| File | Purpose |
|------|---------|
| `run_screening.py` | **Main entry point.** Complete workflow: query MP -> annotate melting points -> rank candidates |
| `screen_carbon_rich_compounds.py` | Step 1: Query Materials Project API for carbon-rich compounds in ternary systems |
| `predict_melting_point.py` | Step 2: Annotate compounds with melting point estimates (lookup + empirical) |
| `requirements.txt` | Python dependencies |

---

## run_screening.py

**Complete screening workflow.** Combines MP API querying, melting point annotation, and candidate ranking into a single pipeline.

### Usage

```bash
python run_screening.py --api-key YOUR_KEY [OPTIONS]
```

### Flags

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--api-key` | Yes | — | Your Materials Project API key |
| `--systems` | No | — | Specific chemical systems to query (e.g., `La-B-C Fe-B-C`). Mutually exclusive with `--metal-group`/`--nonmetal` |
| `--metal-group` | No | — | Metal element group for systematic screening. See [Metal Groups](#metal-groups) below |
| `--nonmetal` | No | — | Non-metal partner element (e.g., `B`, `N`, `Si`), or `all` to screen all of `[B, N, Si, P, S]` |
| `--min-c-fraction` | No | `0.25` | Minimum carbon atomic fraction (0 to 1). Compounds below this threshold are excluded |
| `--max-ehull` | No | `0.1` | Maximum energy above convex hull in eV/atom. Controls thermodynamic stability filter. Increase (e.g., `0.3`) to include metastable phases |
| `--experimental-only` | No | `False` | Only include experimentally synthesized compounds (excludes theoretical/predicted structures from MP) |
| `--mp-min` | No | `2000` | Minimum melting point filter in degrees C |
| `--mp-max` | No | `2500` | Maximum melting point filter in degrees C |
| `--output-prefix` | No | `screening_<timestamp>` | Prefix for output CSV filenames |

### Examples

```bash
# Screen all rare earth borocarbides (RE-B-C systems)
python run_screening.py --api-key YOUR_KEY --metal-group rare_earth --nonmetal B

# Screen specific systems only
python run_screening.py --api-key YOUR_KEY --systems La-B-C Ce-B-C Y-B-C Sc-B-C

# Screen all rare earth systems with B, N, and Si as non-metal partner
python run_screening.py --api-key YOUR_KEY --metal-group rare_earth --nonmetal all

# Broader search: include metastable phases, lower C threshold
python run_screening.py --api-key YOUR_KEY --metal-group rare_earth --nonmetal B \
    --max-ehull 0.3 --min-c-fraction 0.15

# Only experimentally known compounds in transition metal systems
python run_screening.py --api-key YOUR_KEY --metal-group transition_3d --nonmetal B \
    --experimental-only

# Custom melting point window
python run_screening.py --api-key YOUR_KEY --systems La-B-C --mp-min 1800 --mp-max 2200
```

### Output Files

The script generates three CSV files (all prefixed with `--output-prefix` or a timestamp):

| File | Content |
|------|---------|
| `<prefix>_all.csv` | All compounds found, ranked by composite score |
| `<prefix>_best_per_system.csv` | Only the highest-carbon compound per chemical system |
| `<prefix>_in_mp_range.csv` | Compounds with melting points inside the target window |

### Composite Scoring

Candidates are ranked by a weighted composite score:

| Component | Weight | Meaning |
|-----------|--------|---------|
| Carbon atomic fraction | 0.5 | Higher C content = more C available for graphite precipitation |
| Thermodynamic stability | 0.3 | Lower energy above hull = more likely to form |
| Melting point proximity | 0.2 | Closer to center of target mp range = better match to process conditions |

---

## screen_carbon_rich_compounds.py

**Materials Project API querying module.** Queries ternary M-X-C chemical systems and filters for carbon-rich phases.

### Usage (standalone)

```bash
python screen_carbon_rich_compounds.py --api-key YOUR_KEY [OPTIONS]
```

### Flags

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--api-key` | Yes | — | Your Materials Project API key |
| `--systems` | No | — | Specific chemical systems to query (e.g., `La-B-C Fe-B-C`). Space-separated list |
| `--metal-group` | No | — | Metal element group for systematic screening. See [Metal Groups](#metal-groups) below |
| `--nonmetal` | No | — | Non-metal partner element, or `all` for `[B, N, Si, P, S]` |
| `--min-c-fraction` | No | `0.2` | Minimum carbon atomic fraction (0 to 1) |
| `--max-ehull` | No | `0.1` | Maximum energy above hull in eV/atom. Set higher (e.g., `0.3`) for metastable phases |
| `--experimental-only` | No | `False` | Only return experimentally synthesized compounds. Uses MP API `theoretical=False` filter |
| `--output` | No | `carbon_rich_compounds.csv` | Output CSV filename |

### Examples

```bash
# Query La-B-C system with default settings
python screen_carbon_rich_compounds.py --api-key YOUR_KEY --systems La-B-C

# All rare earth + boron systems, only experimental compounds
python screen_carbon_rich_compounds.py --api-key YOUR_KEY \
    --metal-group rare_earth --nonmetal B --experimental-only

# Include metastable phases with very low carbon threshold
python screen_carbon_rich_compounds.py --api-key YOUR_KEY \
    --systems Hf-Ta-C Zr-Nb-C --max-ehull 0.5 --min-c-fraction 0.1
```

### Output Columns

The output CSV contains the following columns:

| Column | Description |
|--------|-------------|
| `material_id` | Materials Project ID (e.g., `mp-1234`) |
| `formula` | Pretty chemical formula |
| `chemsys` | Chemical system (e.g., `B-C-La`) |
| `C_atomic_fraction` | Atomic fraction of carbon (0 to 1) |
| `C_weight_fraction` | Weight fraction of carbon (0 to 1) |
| `n_elements` | Number of distinct elements |
| `energy_above_hull_eV` | Energy above convex hull (eV/atom); 0 = thermodynamically stable |
| `formation_energy_eV` | Formation energy per atom (eV/atom) |
| `spacegroup` | Space group symbol |
| `crystal_system` | Crystal system (cubic, hexagonal, etc.) |
| `density_g_cm3` | Density in g/cm3 |
| `theoretical` | `True` if the structure is theoretical/predicted, `False` if experimentally observed |

---

## predict_melting_point.py

**Melting point annotation module.** Since Materials Project does not store melting points, this module provides estimation through two strategies:

1. **Curated lookup** — hardcoded experimental melting points for ~40 common carbides, borides, borocarbides, and related compounds (from ASM International, NIST, published literature).
2. **Empirical estimation** — rough correlation from formation energy: `T_m ~ 1500 + 2500 * |DeltaH_f|`, with corrections for element count and density. Use as a first-pass filter only (R2 ~ 0.5).

### Usage (standalone)

```bash
python predict_melting_point.py --input INPUT_CSV [OPTIONS]
```

### Flags

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--input` | Yes | — | Input CSV file (output from `screen_carbon_rich_compounds.py`) |
| `--output` | No | `compounds_with_mp.csv` | Output CSV with melting point annotations added |
| `--mp-min` | No | `2000` | Minimum melting point for filtering (degrees C) |
| `--mp-max` | No | `2500` | Maximum melting point for filtering (degrees C) |

### Examples

```bash
# Annotate and filter
python predict_melting_point.py --input carbon_rich_compounds.csv

# Custom melting point window
python predict_melting_point.py --input carbon_rich_compounds.csv \
    --mp-min 1800 --mp-max 3000 --output high_mp_compounds.csv
```

### Output Columns (added)

Two columns are appended to the input CSV:

| Column | Description |
|--------|-------------|
| `melting_point_C` | Estimated or looked-up melting point in degrees C |
| `mp_source` | Source of the value: `experimental_curated` or `empirical_estimate` |

### Curated Melting Point Database

The built-in lookup table currently covers:

| Category | Examples | Count |
|----------|----------|-------|
| Binary carbides | TiC, ZrC, HfC, TaC, SiC, B4C, WC, ... | 20 |
| Binary borides | TiB2, ZrB2, HfB2, LaB6, ... | 10 |
| Ternary borocarbides | ThBC, UBC, UB2C, LaBC, La(BC)2, ... | 6 |
| Carbonitrides | TiCN, HfCN | 2 |
| MAX phases | Ti3SiC2, Ti2AlC, Cr2AlC | 3 |

To add new entries, edit the `KNOWN_MELTING_POINTS` dictionary in `predict_melting_point.py`.

### For More Accurate Predictions

The empirical estimate is a rough heuristic. For production-quality predictions, consider:

- **MeLting GNN model** (Hong et al., PNAS 2022): https://github.com/atomisticnet/MeLting
- **CALPHAD databases** (e.g., SGTE, COST507) for phase diagram-based melting points
- **Literature search** for specific compound classes

---

## Metal Groups

The `--metal-group` flag accepts the following predefined groups:

| Group Name | Elements |
|------------|----------|
| `rare_earth` | La, Ce, Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Y, Sc |
| `transition_3d` | Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn |
| `transition_4d` | Zr, Nb, Mo, Ru, Rh, Pd |
| `transition_5d` | Hf, Ta, W, Re, Os, Ir, Pt |
| `actinide` | Th, U |
| `alkaline_earth` | Ca, Sr, Ba |
| `alkali` | Li, Na, K |

When combined with `--nonmetal`, the script generates all M-X-C ternary combinations. For example, `--metal-group rare_earth --nonmetal B` generates 16 systems: La-B-C, Ce-B-C, ..., Sc-B-C.

Using `--nonmetal all` expands to all five non-metals `[B, N, Si, P, S]`, giving 16 x 5 = 80 systems for rare earths.

---

## Typical Workflow

```
1. Run the full pipeline:
   python run_screening.py --api-key KEY --metal-group rare_earth --nonmetal B --experimental-only

2. Review the output CSVs:
   - *_best_per_system.csv  ->  one top compound per RE-B-C system
   - *_in_mp_range.csv      ->  compounds in your 2000-2500 C window

3. For compounds with mp_source='empirical_estimate':
   - Cross-check with literature
   - Run through MeLting GNN model for better predictions
   - Consult CALPHAD phase diagrams

4. For promising candidates, verify:
   - Phase stability at your target processing temperature
   - Carbon precipitation behavior on cooling
   - Compatibility with your carbon matrix
```

---

## References

- Hong et al., "Melting temperature prediction using a graph neural network model", *PNAS* 119(36), e2209630119 (2022). https://www.pnas.org/doi/10.1073/pnas.2209630119
- ML-guided search for energetically favorable metal borocarbide ternary compounds, *Journal of Alloys and Compounds* (2025). https://www.sciencedirect.com/science/article/abs/pii/S0925838825062401
- Rogl, P., "Phase Equilibria and Structural Chemistry within Ternary Systems: Actinide Metal-Boron-Carbon", Springer (1983).
- Materials Project API: https://materialsproject.org/api
- Pymatgen: https://pymatgen.org
