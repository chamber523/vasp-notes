# BayesianOpt4dftu Tutorial

A tutorial for automatically determining Hubbard U parameters in DFT+U using Bayesian Optimization (BO).

---

## What is BayesianOpt4dftu?

Standard DFT+U requires manually choosing U values, which is system-dependent and often tuned empirically. BayesianOpt4dftu automates this by:

1. Running a high-accuracy **baseline** calculation (HSE06 or GW) to obtain a reference band structure
2. Iteratively running **DFT+U** calculations with different U values
3. Using **Bayesian Optimization** to minimize a loss function that measures the difference between DFT+U and the baseline band structure
4. Returning the optimal U values that best reproduce the reference

The loss function combines:
- **Δgap**: difference in band gap between DFT+U and baseline
- **Δband**: difference in band structure shape (weighted by k-point density)
- **Δmag** *(optional)*: difference in magnetic moments

---

## Installation

```bash
pip install git+https://github.com/caizefeng/BayesianOpt4dftu.git
```

Dependencies installed automatically: `numpy`, `pandas`, `ase`, `pymatgen`, `bayesian-optimization`, `vaspvis`.

VASP must be installed and accessible separately.

---

## Workflow Overview

```
HSE/GW baseline calculation (reference band structure)
        |
        v
Bayesian Optimization loop:
    DFT+U SCF  →  DFT+U Band  →  evaluate loss (Δgap + Δband + Δmag)
        ↑_______________|
        (BO proposes next U values)
        |
        v
Output: optimal U values + BO convergence plot
```

### Directory Structure

```
working_dir/
├── input.json          # Main configuration file
├── hse/                # Baseline HSE calculation (pre-computed)
│   ├── scf/
│   └── band/
├── dftu/               # DFT+U calculations (auto-generated each iteration)
│   ├── scf/
│   └── band/
├── u_xxx.txt           # U values and objective scores per iteration
├── formatted_u_xxx.txt # Formatted version of above
└── 1D_xxx.png / 2D_xxx.png  # BO convergence visualization
```

---

## Configuration: `input.json`

### `vasp_env` — VASP Environment

| Key | Description | Example |
|-----|-------------|---------|
| `vasp_run_command` | Command to run VASP | `"mpirun -np 64 /path/to/vasp_std"` |
| `out_file_name` | VASP stdout filename | `"slurm-vasp.out"` |
| `vasp_pp_path` | Path to pseudopotential directory (containing `potpaw_PBE/`) | `"/path/to/Pseudopotentials/"` |
| `dry_run` | Generate input files only, without running VASP | `false` |
| `dftu_only` | Skip baseline calculation; use pre-computed HSE/GW from `<working_dir>/hse/` | `false` |
| `get_optimal_band` | Run one final DFT+U calculation at the optimal U after BO finishes | `true` |

> **Tip**: Set `dftu_only: true` when the HSE baseline is already computed (e.g., submitted as a separate job). Place the completed HSE results in `<working_dir>/hse/scf/` and `<working_dir>/hse/band/`.

---

### `bo` — Bayesian Optimization Settings

| Key | Default | Description |
|-----|---------|-------------|
| `resume_checkpoint` | `false` | Resume from a previous run saved in `u_tmp.txt` and `input_tmp.json` |
| `baseline` | `"hse"` | Reference method: `"hse"` or `"gw"` (`"gw"` requires `dftu_only: true`) |
| `which_u` | — | Which elements to optimize; `1` = optimize, `0` = fix at 0. Format: one integer per species. For a unary substance: `[1,]` |
| `br` | `[5, 5]` | Band range: `[n_valence, n_conduction]` bands from the Fermi level included in Δband |
| `kappa` | `5` | Exploration–exploitation trade-off. Low (~0) = exploit known optima; high (~10) = explore new regions |
| `alpha_gap` | `0.25` | Weight for Δgap in the loss function |
| `alpha_band` | `0.75` | Weight for Δband in the loss function |
| `alpha_mag` | `0.0` | Weight for Δmagnetization; requires `LORBIT` set in all INCAR files. Set to `0` to exclude |
| `mag_axis` | `"all"` | Magnetic moment component for Δmag: `"x"`, `"y"`, `"z"`, or `"all"`. Only relevant for non-collinear calculations |
| `threshold` | `0.0001` | Stop BO when objective function changes by less than this between consecutive iterations. Set to `0.0` to disable |
| `threshold_opt_u` | `0.0` | Stop BO when optimal U values change by less than this over `report_optimum_interval` iterations. Set to `0.0` to disable |
| `urange` | `[-10, 10]` | Search range for U values (eV). Same range applied to all optimized elements |
| `elements` | — | Element symbols in your system, used for plot labels. E.g., `["In", "As"]` |
| `iteration` | `50` | Maximum number of BO iterations |
| `report_optimum_interval` | `10` | Interval (iterations) at which optimal U is extracted from the GP mean and logged |
| `print_magmom` | `false` | Print magnetic moment at each iteration |

---

### `structure_info` — Crystal Structure

| Key | Description |
|-----|-------------|
| `lattice_param` | Lattice constant (Å); corresponds to the scale factor in POSCAR |
| `cell` | 3×3 matrix of lattice vectors (in units of `lattice_param`) |
| `atoms` | List of `[element, position, magmom]` for each atom. `magmom` is a scalar for collinear, or `[mx, my, mz]` for non-collinear. Use a small non-zero value (e.g., `1e-6`) instead of exact `0` to avoid ASE errors |
| `kgrid_hse` | SCF k-grid for HSE baseline, e.g., `[7, 7, 7]` |
| `kgrid_pbe` | SCF k-grid for DFT+U, e.g., `[7, 7, 7]` |
| `num_kpts` | Number of k-points per segment along the band path. Use `"auto"` to match the baseline |
| `kpath` | High-symmetry k-path string, e.g., `"G X W L G K"` |
| `custom_kpoints` | Override default BZ coordinates for specific k-points, e.g., `{"H": [0.5, -0.5, 0.5]}`. Set to `null` to use defaults |
| `custom_POTCAR_path` | Path to a custom POTCAR file. Set to `null` to auto-generate from `vasp_pp_path` |

#### Example: InAs

```json
"structure_info": {
    "lattice_param": 6.0584,
    "cell": [
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0]
    ],
    "atoms": [
        ["In", [0, 0, 0],       [0, 0, 1e-6]],
        ["As", [0.75, 0.75, 0.75], [0, 0, 1e-6]]
    ],
    "kgrid_hse": [7, 7, 7],
    "kgrid_pbe": [7, 7, 7],
    "num_kpts": 50,
    "kpath": "G X W L G K",
    "custom_kpoints": null,
    "custom_POTCAR_path": null
}
```

---

### INCAR Flags Sections

These sections use lowercase ASE VASP calculator keys (mostly the same as INCAR tags but lowercase). See the [ASE VASP calculator docs](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html) for the full list.

| Section | Applied to |
|---------|-----------|
| `general_flags` | All VASP calculations |
| `scf` | SCF calculations only |
| `band` | Band structure calculations only |
| `pbe` | DFT+U calculations only |
| `hse` | HSE baseline calculations only |

#### Example: `general_flags`

```json
"general_flags": {
    "encut": 400,
    "sigma": 0.05,
    "ediff": 1e-5,
    "prec": "N",
    "algo": "D",
    "lsorbit": true,
    "saxis": [0, 0, 1],
    "nbands": 21,
    "kpar": 6,
    "ncore": 3,
    "lmaxmix": 4,
    "lorbit": 11
}
```

#### Example: `pbe` (DFT+U specific)

The `ldau_luj` key sets per-element U and J values. `L` is the angular momentum channel (1 = p, 2 = d, 3 = f).

```json
"pbe": {
    "xc": "pbe",
    "ldau": true,
    "ldau_luj": {
        "In": {"L": 1, "U": 0.0, "J": 0.0},
        "As": {"L": 1, "U": 0.0, "J": 0.0}
    }
}
```

> **Note**: The U values here are just placeholders. BO will override them at each iteration.

---

## Running

```bash
cd /path/to/working_dir
bo_dftu
```

To resume from a checkpoint:

```bash
# Set "resume_checkpoint": true in input.json, then:
bo_dftu
```

---

## Output Files

| File | Description |
|------|-------------|
| `u_xxx.txt` | Tab-separated log: U values, band gap, Δgap, Δband, Δmag per iteration |
| `formatted_u_xxx.txt` | Human-readable version of above |
| `1D_xxx.png` | BO visualization for 1 optimized element: GP predicted mean + acquisition function |
| `2D_xxx.png` | BO visualization for 2 optimized elements |
| `u_tmp.txt` / `input_tmp.json` | Checkpoint files for resuming |

> **Note**: Visualization is only generated for 1 or 2 optimized U parameters. Three or more will skip the plot.

The filename suffix encodes BO settings, e.g.:
```
2D_kappa_5.0_ag_0.25_ab_0.75_am_0.0.png
       ↑          ↑       ↑       ↑
     kappa    alpha_gap  alpha_band  alpha_mag
```

---

## Tips

- **Start with `dry_run: true`** to verify that input files are generated correctly before launching real VASP jobs.
- **Alpha weights**: `alpha_band` typically carries more information than `alpha_gap`. The default `0.25/0.75` split is a reasonable starting point.
- **`kappa`**: If BO converges too quickly to a local optimum, increase `kappa` (more exploration). If convergence is slow, decrease it.
- **`urange`**: Narrow the range if you have prior knowledge (e.g., `[0, 8]` for transition metal d-electrons) to speed up convergence.
- **`LMAXMIX`**: Must be set correctly — `4` for d-electron systems, `6` for f-electron systems — for proper charge density mixing with DFT+U.

---

## References

Yu et al. (2020) — "Machine learning the Hubbard U parameter in DFT+U using Bayesian optimization", *npj Computational Materials*, 6(1):1–6.
