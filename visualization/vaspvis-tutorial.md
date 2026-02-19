# vaspvis Tutorial

vaspvis is a Python library for visualizing electronic structure data (band structure and DOS) from VASP calculations.

---

## Installation

```bash
pip install vaspvis
```

---

## Required VASP Files

Each calculation folder must contain:

| Plot type | Required files |
|-----------|---------------|
| Band structure | `EIGENVAL`, `KPOINTS`, `POSCAR`, `INCAR`, `OUTCAR` |
| Projected band | above + `PROCAR` (need `LORBIT = 11` in INCAR) |
| DOS | `DOSCAR`, `POSCAR`, `INCAR`, `OUTCAR` |
| Projected DOS | above, `LORBIT >= 11` in INCAR |

**Note:** vaspvis caches parsed data as `.npy` files in the same folder. Subsequent loads are much faster.

---

## Two Usage Modes

### Mode 1 — `standard` module (quick, no matplotlib needed)

```python
from vaspvis import standard

standard.band_plain(folder='./band')
```

One line → figure saved as PNG.

### Mode 2 — `Band` / `Dos` classes (flexible, pass your own axes)

```python
from vaspvis import Band, Dos
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(5, 4))
bs = Band(folder='./band', projected=True)
bs.plot_spd(ax)
plt.tight_layout()
plt.savefig('band_spd.png', dpi=300)
```

---

## Band Structure

### Plain band structure

```python
from vaspvis import standard

standard.band_plain(
    folder='./band',
    output='band_plain.png',
    erange=[-4, 4],          # energy window (eV), relative to E_F
    color='black',
    linewidth=1.25,
    figsize=(4, 3),
    fontsize=12,
)
```

### s/p/d projected

```python
standard.band_spd(
    folder='./band',
    output='band_spd.png',
    erange=[-4, 4],
    scale_factor=5,           # size of scatter points
    orbitals='spd',           # 's', 'p', 'd', 'f', or combinations like 'sp'
)
```

### Orbital projected

Orbital indices:

| Index | Orbital | Index | Orbital |
|-------|---------|-------|---------|
| 0 | s | 4 | d_xy |
| 1 | p_y | 5 | d_yz |
| 2 | p_z | 6 | d_z² |
| 3 | p_x | 7 | d_xz |
| — | — | 8 | d_x²-y² |

```python
standard.band_orbitals(
    folder='./band',
    orbitals=[0, 1, 2, 3, 4, 5, 6, 7, 8],  # which orbitals to show
    erange=[-4, 4],
)
```

### Atom projected

Atom indices follow the order in POSCAR (0-based).

```python
standard.band_atoms(
    folder='./band',
    atoms=[0, 1],             # atom indices from POSCAR
    erange=[-4, 4],
)
```

### Atom-orbital projected

```python
standard.band_atom_orbitals(
    folder='./band',
    atom_orbital_dict={
        0: [1, 3],            # atom 0 → p_y and p_x
        1: [4, 8],            # atom 1 → d_xy and d_x²-y²
    },
    erange=[-4, 4],
)
```

### Element projected

```python
standard.band_elements(
    folder='./band',
    elements=['In', 'As'],
    erange=[-4, 4],
)
```

### Element s/p/d projected

```python
standard.band_element_spd(
    folder='./band',
    element_spd_dict={'As': 'spd'},   # project spd onto As only
    erange=[-4, 4],
)
```

### Element-orbital projected

```python
standard.band_element_orbitals(
    folder='./band',
    element_orbital_dict={
        'As': [2],            # As → p_z
        'In': [3],            # In → p_x
    },
    erange=[-4, 4],
)
```

---

## Density of States

`energyaxis` controls whether energy is on the x or y axis.

### Plain DOS

```python
standard.dos_plain(
    folder='./dos',
    energyaxis='x',           # 'x' (horizontal energy) or 'y' (vertical)
    erange=[-4, 4],
)
```

### s/p/d projected DOS

```python
standard.dos_spd(
    folder='./dos',
    energyaxis='x',
    erange=[-4, 4],
)
```

### Orbital projected DOS

```python
standard.dos_orbitals(
    folder='./dos',
    orbitals=[0, 1, 2, 3, 4, 5, 6, 7, 8],
    energyaxis='x',
    erange=[-4, 4],
)
```

### Atom projected DOS

```python
standard.dos_atoms(
    folder='./dos',
    atoms=[0, 1],
    energyaxis='x',
)
```

### Atom-orbital projected DOS

```python
standard.dos_atom_orbitals(
    folder='./dos',
    atom_orbital_dict={0: [1, 3], 1: [4, 8]},
    energyaxis='x',
)
```

### Element projected DOS

```python
standard.dos_elements(
    folder='./dos',
    elements=['In', 'As'],
    energyaxis='x',
)
```

### Element s/p/d projected DOS

```python
standard.dos_element_spd(
    folder='./dos',
    element_spd_dict={'As': 'spd'},
    energyaxis='x',
)
```

### Layer-resolved DOS (for slabs)

```python
standard.dos_layers(
    folder='./dos',
    energyaxis='x',
    erange=[-4, 4],
)
```

---

## Combined Band + DOS

All band functions have a `band_dos_*` equivalent. Pass two separate folders.

```python
# Plain
standard.band_dos_plain(
    band_folder='./band',
    dos_folder='./dos',
    erange=[-4, 4],
)

# s/p/d projected
standard.band_dos_spd(
    band_folder='./band',
    dos_folder='./dos',
    erange=[-4, 4],
)

# Element projected
standard.band_dos_elements(
    band_folder='./band',
    dos_folder='./dos',
    elements=['In', 'As'],
    erange=[-4, 4],
)

# Element s/p/d
standard.band_dos_element_spd(
    band_folder='./band',
    dos_folder='./dos',
    element_spd_dict={'As': 'spd'},
    erange=[-4, 4],
)
```

---

## Custom Plots (Band/Dos classes)

Use this when you need full control over the figure layout.

```python
from vaspvis import Band, Dos
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(7, 4),
                         gridspec_kw={'width_ratios': [2, 1]})
ax_band, ax_dos = axes

# Band
bs = Band(folder='./band', projected=True)
bs.plot_element_spd(ax_band, element_spd_dict={'As': 'spd'})
ax_band.set_ylim(-4, 4)
ax_band.set_ylabel('$E - E_F$ (eV)')

# DOS (energyaxis='y' to match band structure)
dos = Dos(folder='./dos')
dos.plot_spd(ax_dos, energyaxis='y')
ax_dos.set_ylim(-4, 4)
ax_dos.set_yticklabels([])

plt.tight_layout()
plt.savefig('custom_band_dos.png', dpi=300)
```

---

## Spin-Polarized Calculations

For `ISPIN = 2` calculations, load spin channels separately:

```python
# Standard module — spin-polarized variants
standard.band_plain_spin_polarized(folder='./band', erange=[-4, 4])
standard.dos_plain_spin_polarized(folder='./dos', energyaxis='x')
standard.band_dos_plain_spin_polarized(band_folder='./band', dos_folder='./dos')

# Or load manually
bs_up   = Band(folder='./band', spin='up')
bs_down = Band(folder='./band', spin='down')

dos_up   = Dos(folder='./dos', spin='up')
dos_down = Dos(folder='./dos', spin='down')
dos_both = Dos(folder='./dos', spin='both',
               combination_method='sub')   # magnetization = up - down
```

---

## Spin-Orbit Coupling (SOC)

For `LSORBIT = .TRUE.` calculations, select a spin component to track:

```python
bs  = Band(folder='./band', soc_axis='z')  # track S_z component
dos = Dos(folder='./dos',   soc_axis='z')
```

`spin='up'` → positive S_z, `spin='down'` → negative S_z.

---

## Band Unfolding (Supercell / Slab)

Useful for visualizing band structures of supercells back in the primitive BZ.

### Step 1 — Get transformation matrix

```python
from vaspvis.utils import convert_slab

M = convert_slab(
    bulk_path='POSCAR_bulk',   # primitive bulk POSCAR
    slab_path='POSCAR_slab',   # slab POSCAR
    index=[1, 1, 1],           # Miller index
)
# M is printed and returned as a 3×3 integer matrix
```

### Step 2 — Generate KPOINTS

```python
from vaspvis.utils import generate_kpoints

high_symm_points = [
    [0.5, 0.0, 0.5],   # X
    [0.0, 0.0, 0.0],   # Gamma
    [0.5, 0.0, 0.5],   # X
]

generate_kpoints(
    M=M,
    high_symmetry_points=high_symm_points,
    n=50,               # k-points per segment
)
```

### Step 3 — Run VASP, then plot

```python
standard.band_plain(
    folder='./band_unfold',
    unfold=True,
    M=M,
    kpath='XGX',
    high_symm_points=high_symm_points,
    n=50,
    erange=[-4, 4],
    heatmap=True,       # density heatmap (recommended for unfolded)
)
```

---

## Common Parameters Reference

| Parameter | Type | Description |
|-----------|------|-------------|
| `folder` | str | Path to VASP output folder |
| `erange` | [min, max] | Energy window relative to E_F (eV) |
| `figsize` | (w, h) | Figure size in inches |
| `fontsize` | float | Font size |
| `output` | str | Output filename |
| `save` | bool | If False, returns `(fig, ax)` instead of saving |
| `spin` | str | `'up'`, `'down'`, or `'both'` |
| `scale_factor` | float | Size of projected scatter points |
| `energyaxis` | str | `'x'` or `'y'` for DOS plots |
| `shift_efermi` | float | Manually shift the Fermi energy (eV) |
| `efermi_folder` | str | Use OUTCAR from a separate SCF folder for E_F |

---

## Clearing Cache

If you rerun a VASP calculation and want vaspvis to re-parse the files, delete the cached `.npy` files:

```bash
rm ./band/eigenvalues.npy ./band/projected_eigenvalues.npy
rm ./dos/dos.npy ./dos/projected_dos.npy
```
