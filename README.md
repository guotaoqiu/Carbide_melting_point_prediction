# Carbide Melting Point Prediction & Carbon-Rich Compound Screening

High-throughput screening toolkit for identifying carbon-rich compounds in various chemical systems as candidates for graphitization catalysis research.

## Background

The catalytic effect of LaB6 on graphitization is believed to stem from the formation of La(BC)2 — a high-carbon borocarbide — from LaB6 at high temperature in a carbon-rich environment. As carbon precipitates, it adheres to La(BC)2 surfaces, gradually forming graphite. This toolkit systematically finds the **highest-carbon-content compound** in each chemical system with melting points in a target range (default: 2000-2500 C), replacing the inefficient manual browsing of Materials Project phase diagrams.

The screening supports multiple system types:
- **Binary M-C** — simple metal carbides (e.g., LaC2, TiC)
- **Ternary M-X-C** — metal + non-metal + carbon (e.g., La-B-C, Ti-Si-C)
- **Bimetal M1-M2-C** — two metals + carbon (e.g., Ta-Hf-C)
- **All combined** — comprehensive search across binary, ternary, and bimetal systems

## Installation

```bash
pip install -r requirements.txt
```

**Requirements:**
- Python >= 3.10
- `mp-api` >= 0.41.0 (Materials Project API client)
- `pymatgen` >= 2024.1.1 (Python Materials Genomics)
- `pandas` >= 2.0
- `requests` >= 2.28 (for MAPP API calls)

You also need a **Materials Project API key**. Get one free at: https://materialsproject.org/api

## Files Overview

| File | Purpose |
|------|---------|
| `run_screening.py` | **Main entry point.** Complete workflow: query MP -> annotate melting points -> rank candidates |
| `screen_carbon_rich_compounds.py` | Step 1: Query Materials Project API for carbon-rich compounds across multiple search modes |
| `predict_melting_point.py` | Step 2: Annotate compounds with melting point estimates (curated lookup, MAPP GNN, or empirical) |
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
| `--systems` | No | — | Explicit chemical systems to query (e.g., `La-B-C Hf-Ta-C La-C`). Bypasses `--mode` when provided |
| `--mode` | No | `ternary` | Search mode. Controls what types of systems to generate. See [Search Modes](#search-modes) |
| `--metal-group` | No | — | Primary metal element group. See [Metal Groups](#metal-groups) |
| `--partner-group` | No | — | Partner element group for ternary/bimetal modes. See [Partner Groups](#partner-groups) |
| `--partner-elements` | No | — | Explicit partner elements, space-separated (overrides `--partner-group`). E.g., `B N Si Hf Ta` |
| `--min-c-fraction` | No | `0.25` | Minimum carbon atomic fraction (0 to 1). Compounds below this threshold are excluded |
| `--max-ehull` | No | `0.1` | Maximum energy above convex hull in eV/atom. Controls thermodynamic stability filter. Increase (e.g., `0.3`) to include metastable phases |
| `--experimental-only` | No | `False` | Only include experimentally synthesized compounds (excludes theoretical/predicted structures from MP) |
| `--use-mapp` | No | `False` | Enable MAPP GNN melting point predictions (Hong et al., PNAS 2022). Requires internet. See [MAPP GNN Predictions](#mapp-gnn-predictions) |
| `--mapp-url` | No | `http://206.207.50.58:5007/...` | Custom MAPP API endpoint URL (only needed if self-hosting the MAPP server) |
| `--mp-min` | No | `2000` | Minimum melting point filter in degrees C |
| `--mp-max` | No | `2500` | Maximum melting point filter in degrees C |
| `--output-prefix` | No | `screening_<timestamp>` | Prefix for output CSV filenames |

### Examples

```bash
# Full search: binary + ternary + bimetal for rare earths with non-metal partners
python run_screening.py --api-key KEY --mode all --metal-group rare_earth --partner-group nonmetal

# Same but with MAPP GNN melting point predictions
python run_screening.py --api-key KEY --mode all --metal-group rare_earth \
    --partner-group nonmetal --use-mapp

# Just binary rare earth carbides (M-C only, no partner needed)
python run_screening.py --api-key KEY --mode binary --metal-group rare_earth

# Ternary: rare earth + boron + carbon
python run_screening.py --api-key KEY --mode ternary --metal-group rare_earth --partner-elements B

# Bimetal: rare earth + refractory metals + carbon (e.g., La-Hf-C, Ce-Ta-C)
python run_screening.py --api-key KEY --mode bimetal --metal-group rare_earth --partner-group transition_5d

# Ternary with extended non-metals (includes O, Se, Te, F, Cl)
python run_screening.py --api-key KEY --mode ternary --metal-group rare_earth --partner-group nonmetal_extended

# Custom partner elements: mix metals and non-metals freely
python run_screening.py --api-key KEY --mode ternary --metal-group rare_earth --partner-elements B N Hf Ta W

# Specific systems (bypasses --mode entirely)
python run_screening.py --api-key KEY --systems La-B-C Hf-Ta-C La-C Ce-B-C

# Only experimentally known, broader stability window, with ML melting points
python run_screening.py --api-key KEY --mode all --metal-group transition_3d \
    --partner-group nonmetal --experimental-only --max-ehull 0.3 --use-mapp

# Custom melting point window
python run_screening.py --api-key KEY --systems La-B-C --mp-min 1800 --mp-max 2200
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

**Materials Project API querying module.** Queries chemical systems across multiple search modes and filters for carbon-rich phases.

### Usage (standalone)

```bash
python screen_carbon_rich_compounds.py --api-key YOUR_KEY [OPTIONS]
```

### Flags

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--api-key` | Yes | — | Your Materials Project API key |
| `--systems` | No | — | Explicit chemical systems to query (e.g., `La-B-C Hf-Ta-C La-C`). Bypasses `--mode` |
| `--mode` | No | `ternary` | Search mode: `binary`, `ternary`, `bimetal`, or `all`. See [Search Modes](#search-modes) |
| `--metal-group` | No | — | Primary metal element group. See [Metal Groups](#metal-groups) |
| `--partner-group` | No | — | Partner element group. See [Partner Groups](#partner-groups) |
| `--partner-elements` | No | — | Explicit partner elements (overrides `--partner-group`). E.g., `--partner-elements B N Si Hf Ta` |
| `--min-c-fraction` | No | `0.2` | Minimum carbon atomic fraction (0 to 1) |
| `--max-ehull` | No | `0.1` | Maximum energy above hull in eV/atom. Set higher (e.g., `0.3`) for metastable phases |
| `--experimental-only` | No | `False` | Only return experimentally synthesized compounds. Uses MP API `theoretical=False` filter |
| `--output` | No | `carbon_rich_compounds.csv` | Output CSV filename |

### Examples

```bash
# Binary: all rare earth carbides
python screen_carbon_rich_compounds.py --api-key KEY --mode binary --metal-group rare_earth

# Ternary: rare earth + boron + carbon
python screen_carbon_rich_compounds.py --api-key KEY --mode ternary --metal-group rare_earth --partner-elements B

# Bimetal: transition metals paired with each other + carbon
python screen_carbon_rich_compounds.py --api-key KEY --mode bimetal --metal-group transition_3d --partner-group transition_5d

# Everything: binary + ternary + bimetal for rare earths
python screen_carbon_rich_compounds.py --api-key KEY --mode all --metal-group rare_earth --partner-group nonmetal

# Explicit systems
python screen_carbon_rich_compounds.py --api-key KEY --systems La-C La-B-C Hf-Ta-C --experimental-only
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
| `theoretical` | `True` if theoretical/predicted, `False` if experimentally observed |

**Note:** Results are deduplicated by `material_id` — if a compound appears in overlapping subsystems (e.g., LaC2 appears in both La-C and La-B-C queries), it is only listed once.

---

## predict_melting_point.py

**Melting point annotation module.** Since Materials Project does not store melting points, this module provides three estimation strategies (in priority order):

1. **Curated lookup** — hardcoded experimental melting points for ~40 common carbides, borides, borocarbides, and related compounds (from ASM International, NIST, published literature).
2. **MAPP GNN prediction** (opt-in via `--use-mapp`) — Hong et al.'s GNN+ResNet ensemble model that predicts melting temperature from chemical formula alone. Uses a remote API server. Returns prediction + uncertainty (standard error from 30-model ensemble). Supports compounds with up to 6 elements.
3. **Empirical estimation** (fallback) — rough correlation from formation energy: `T_m ~ 1500 + 2500 * |DeltaH_f|`, with corrections for element count and density. Use as last-resort only (R2 ~ 0.5).

### Usage (standalone)

```bash
python predict_melting_point.py --input INPUT_CSV [OPTIONS]
```

### Flags

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--input` | Yes | — | Input CSV file (output from `screen_carbon_rich_compounds.py`) |
| `--output` | No | `compounds_with_mp.csv` | Output CSV with melting point annotations added |
| `--use-mapp` | No | `False` | Enable MAPP GNN melting point predictions. Requires internet access. See [MAPP GNN Predictions](#mapp-gnn-predictions) |
| `--mapp-url` | No | `http://206.207.50.58:5007/...` | Custom MAPP API endpoint URL (only needed if self-hosting) |
| `--mp-min` | No | `2000` | Minimum melting point for filtering (degrees C) |
| `--mp-max` | No | `2500` | Maximum melting point for filtering (degrees C) |

### Examples

```bash
# Default: curated lookup + empirical fallback (no internet needed)
python predict_melting_point.py --input carbon_rich_compounds.csv

# Enable MAPP GNN predictions for compounds not in curated database
python predict_melting_point.py --input carbon_rich_compounds.csv --use-mapp

# Custom melting point window + MAPP
python predict_melting_point.py --input carbon_rich_compounds.csv \
    --use-mapp --mp-min 1800 --mp-max 3000 --output high_mp_compounds.csv
```

### Output Columns (added)

Three columns are appended to the input CSV:

| Column | Description |
|--------|-------------|
| `melting_point_C` | Estimated or looked-up melting point in degrees C |
| `mp_source` | Source: `experimental_curated`, `mapp_gnn`, or `empirical_estimate` |
| `mp_std_error_C` | Standard error of prediction in degrees C (only for `mapp_gnn` source, null otherwise) |

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

---

## MAPP GNN Predictions

The MAPP model (Materials-Agnostic Platform for Prediction) is from:

> Hong et al., "Melting temperature prediction using a graph neural network model: From ancient minerals to new materials", *PNAS* 119(36), e2209630119 (2022).

### How it works

- **Model**: Ensemble of 30 GNN + ResNet models trained on ~10,000 compounds with known melting points
- **Input**: Chemical formula only (no crystal structure needed)
- **Output**: Predicted melting temperature (K) + standard error from ensemble disagreement
- **Scope**: Compounds with up to 6 elements
- **Accuracy**: The model's top predictions for highest-melting compounds are overwhelmingly carbides and nitrides, consistent with experimental knowledge. The ensemble standard error provides built-in uncertainty quantification

### Architecture

The MAPP model runs as a **remote API** hosted by Hong's group at ASU. The trained model weights are not publicly distributed — predictions are obtained via HTTP POST to their server.

```
Your formulas (CSV) --> HTTP POST --> MAPP server (30 GNN models) --> Predictions (K) + std error
```

### Usage

```bash
# Enable in the standalone melting point script
python predict_melting_point.py --input compounds.csv --use-mapp

# Enable in the full screening workflow
python run_screening.py --api-key KEY --mode all --metal-group rare_earth \
    --partner-group nonmetal --use-mapp

# Use a self-hosted MAPP server
python predict_melting_point.py --input compounds.csv --use-mapp \
    --mapp-url http://your-server:5007/MT_ML_Qijun_Hong_Predict_noNN
```

### Priority order

When `--use-mapp` is enabled, melting points are determined in this order:
1. **Curated lookup** — if the formula is in our hardcoded database, use the experimental value (most reliable)
2. **MAPP GNN** — for all other formulas, query the MAPP API (ML prediction with uncertainty)
3. **Empirical estimate** — only if MAPP fails (server down, formula too complex, etc.)

Without `--use-mapp`, step 2 is skipped and the empirical fallback is used directly.

### Caveats

- Requires internet access to reach the MAPP API server (`http://206.207.50.58:5007`)
- The server may be temporarily unavailable — the script retries 3 times with exponential backoff
- Predictions are batched (500 formulas per request) to avoid timeouts
- The MAPP model was trained primarily on binary/ternary compounds; predictions for quaternary+ systems may have larger errors
- Always check the `mp_std_error_C` column — large standard errors (>200 deg C) indicate low model confidence

### MAPP API Reference

- **Paper**: https://www.pnas.org/doi/10.1073/pnas.2209630119
- **GitHub**: https://github.com/qjhong/mapp_api
- **Group page**: https://faculty.engineering.asu.edu/hong/melting-temperature-predictor/

---

## Search Modes

The `--mode` flag controls what types of chemical systems are generated:

| Mode | Systems Generated | When to Use |
|------|-------------------|-------------|
| `binary` | M-C | Find simple metal carbides. No `--partner-group` needed |
| `ternary` | M-partner-C | Metal + any partner element + carbon. Requires `--partner-group` or `--partner-elements` |
| `bimetal` | M1-M2-C | Two metals + carbon. If `--partner-group` is a metal group, generates cross-group pairs. If omitted, generates pairs within `--metal-group` |
| `all` | M-C + M-partner-C + M1-M2-C | Comprehensive search. Generates all binary, ternary, and bimetal systems. Requires `--partner-group` or `--partner-elements` for the ternary/bimetal components |

**Note:** `--systems` always takes priority over `--mode`. If you provide explicit systems, mode is ignored.

### How systems are generated

- **binary**: For each metal M in `--metal-group`, generates `M-C`
- **ternary**: For each (M, partner) pair, generates the canonical sorted system (e.g., `B-C-La`)
- **bimetal**: For each (M1, M2) pair across groups, generates `M1-M2-C`. If no partner is specified, generates all pairs within `--metal-group` using combinations
- **all**: Union of all three above

Duplicate systems are automatically removed. Results are deduplicated by `material_id`.

---

## Metal Groups

The `--metal-group` flag selects the primary metal elements:

| Group Name | Elements |
|------------|----------|
| `rare_earth` | La, Ce, Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Y, Sc |
| `transition_3d` | Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn |
| `transition_4d` | Zr, Nb, Mo, Ru, Rh, Pd |
| `transition_5d` | Hf, Ta, W, Re, Os, Ir, Pt |
| `actinide` | Th, U |
| `alkaline_earth` | Ca, Sr, Ba |
| `alkali` | Li, Na, K |

---

## Partner Groups

The `--partner-group` flag selects elements to pair with the primary metals. It can be a non-metal group or another metal group (for bimetal screening):

| Group Name | Elements | Typical Use |
|------------|----------|-------------|
| `nonmetal` | B, N, Si, P, S | Standard ternary M-X-C screening |
| `nonmetal_extended` | B, N, Si, P, S, Se, Te, O, F, Cl | Broader non-metal screening including chalcogenides/halides |
| `rare_earth` | La, Ce, Pr, ... (same as metal group) | Bimetal: pair metals with rare earths |
| `transition_3d` | Ti, V, Cr, ... (same as metal group) | Bimetal: pair with 3d transition metals |
| `transition_4d` | Zr, Nb, Mo, ... (same as metal group) | Bimetal: pair with 4d transition metals |
| `transition_5d` | Hf, Ta, W, ... (same as metal group) | Bimetal: pair with refractory 5d metals |
| `actinide` | Th, U | Bimetal: pair with actinides |
| `alkaline_earth` | Ca, Sr, Ba | Bimetal: pair with alkaline earths |
| `alkali` | Li, Na, K | Bimetal: pair with alkali metals |

**Tip:** Use `--partner-elements` to freely mix metals and non-metals as partners, e.g., `--partner-elements B N Hf Ta`.

---

## Typical Workflow

```
1. Start broad — screen everything for your metal group with ML melting points:
   python run_screening.py --api-key KEY --mode all --metal-group rare_earth \
       --partner-group nonmetal --experimental-only --use-mapp

2. Review the output CSVs:
   - *_best_per_system.csv  ->  one top compound per system
   - *_in_mp_range.csv      ->  compounds in your 2000-2500 C window
   - Check mp_source column:  experimental_curated > mapp_gnn > empirical_estimate
   - Check mp_std_error_C:    large values (>200) = low confidence

3. Drill into interesting bimetal systems:
   python run_screening.py --api-key KEY --mode bimetal --metal-group rare_earth \
       --partner-group transition_5d --max-ehull 0.3 --use-mapp

4. For high-confidence candidates, verify:
   - Phase stability at your target processing temperature
   - Carbon precipitation behavior on cooling
   - Compatibility with your carbon matrix
   - Cross-check MAPP predictions with literature where possible

5. For low-confidence predictions (large std error or empirical_estimate):
   - Search literature for specific compound melting points
   - Consult CALPHAD databases (e.g., SGTE, COST507)
```

---

## References

- Hong et al., "Melting temperature prediction using a graph neural network model", *PNAS* 119(36), e2209630119 (2022). https://www.pnas.org/doi/10.1073/pnas.2209630119
- MAPP API GitHub: https://github.com/qjhong/mapp_api
- ML-guided search for energetically favorable metal borocarbide ternary compounds, *Journal of Alloys and Compounds* (2025). https://www.sciencedirect.com/science/article/abs/pii/S0925838825062401
- Rogl, P., "Phase Equilibria and Structural Chemistry within Ternary Systems: Actinide Metal-Boron-Carbon", Springer (1983).
- Materials Project API: https://materialsproject.org/api
- Pymatgen: https://pymatgen.org
