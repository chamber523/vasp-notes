# Magnetic Calculation

Non-collinear magnetic calculations with SOC for PdCrO₂-based systems.

## Example: PdCrO2 Non-collinear Magnetism

```
magnetic-calculation/
└── example-PdCrO2/
    ├── model/                        # Structure preparation
    │   ├── shift_z.py                 # Reorient cell (Cr to bottom)
    │   ├── make_supercell.py          # Generate supercell
    │   ├── generate_magmom.py         # Calculate magnetic moments (Takatsu model)
    │   ├── generate_mcif.py           # VESTA visualization
    │   ├── 01_POSCAR_unitcell
    │   ├── 02_POSCAR_Cr-bottom
    │   ├── 03_POSCAR_supercell_6x
    │   ├── MAGMOM.txt
    │   └── magnetic_structure.mcif
    │
    └── calculation/                  # DFT calculations
        ├── scf_noncollinear/          # Self-consistent field
        ├── dos_nonllinear/            # Density of states
        └── band_unfold/               # Band structure + unfolding
```

## Workflow

### 1. Structure Preparation

```bash
cd example-PdCrO2/model/

# Reorient to Cr-bottom
python shift_z.py 01_POSCAR_unitcell 02_POSCAR_Cr-bottom

# Generate supercell (transformation: [[2,1,0],[1,2,0],[0,0,2]])
python make_supercell.py 02_POSCAR_Cr-bottom 03_POSCAR_supercell_6x

# Calculate magnetic moments
python generate_magmom.py

# Generate VESTA file
python generate_mcif.py 03_POSCAR_supercell_6x
```

### 2. DFT Calculations

```bash
# SCF
cd ../calculation/scf_noncollinear/
sbatch submit_vasp6.4.3_cpu.slurm

# DOS (after SCF converges)
cd ../dos_nonllinear/
cp ../scf_noncollinear/CHGCAR .
sbatch submit_vasp6.4.3_cpu.slurm

# Band structure
cd ../band_unfold/
cp ../scf_noncollinear/CHGCAR .
sbatch submit_vasp6.4.3_cpu.slurm
```

## Key Parameters

### Non-collinear Magnetism

In VASP, collinear spin, non-collinear magnetism, and SOC are related but distinct concepts, controlled by two INCAR tags: `LNONCOLLINEAR` and `LSORBIT`.

---

#### Collinear Spin

```
ISPIN         = 2
LNONCOLLINEAR = .FALSE.
LSORBIT       = .FALSE.
```

- Spin is restricted to ↑ and ↓ along a fixed axis (typically z)
- The two spin channels are solved independently
- `MAGMOM` is a scalar per atom (e.g., `MAGMOM = 3 -3 0`)
- Computational cost: lowest
- Use for: ferromagnetism, simple antiferromagnetism, most routine magnetic calculations

---

#### Non-collinear Magnetism (without SOC)

```
LNONCOLLINEAR = .TRUE.
LSORBIT       = .FALSE.
```

- Spin can point in any direction; magnetic moment is a 3D vector **m** = (m_x, m_y, m_z)
- `MAGMOM` requires 3 components per atom (e.g., `MAGMOM = 0 0 3`)
- `ISPIN` is ignored — spin is no longer a simple ↑/↓ flag but a full 3D vector, so `ISPIN` becomes meaningless. Do not set `ISPIN = 2` alongside `LNONCOLLINEAR = .TRUE.` (will cause an error in VASP ≥ 6.5)
- Computational cost: ~2× collinear
- Use for: spin spirals, skyrmions, frustrated magnetism, canted antiferromagnets

---

#### SOC (Spin–Orbit Coupling)

```
LSORBIT       = .TRUE.
LNONCOLLINEAR = .TRUE.    # automatically enabled when LSORBIT = .TRUE.
SAXIS         = 0 0 1     # spin quantization axis
```

- SOC couples the orbital (**L**) and spin (**S**) degrees of freedom via the **L·S** interaction
- Spin is no longer a good quantum number; ↑ and ↓ are mixed
- `LNONCOLLINEAR` is automatically set to `.TRUE.` when `LSORBIT = .TRUE.` — no need to set it manually
- `MAGMOM` requires 3 components per atom, same as non-collinear
- `SAXIS` defines the spin quantization axis; default is `0 0 1` (z-axis)
- Computational cost: ~2–4× collinear
- Use for: band splitting, magnetic anisotropy energy (MAE), topological properties

---

#### Summary

| | Collinear | Non-collinear | SOC |
|---|---|---|---|
| `LNONCOLLINEAR` | `.FALSE.` | `.TRUE.` | `.TRUE.` (auto) |
| `LSORBIT` | `.FALSE.` | `.FALSE.` | `.TRUE.` |
| MAGMOM format | scalar per atom | 3 components per atom | 3 components per atom |
| Spin channels | independent | coupled | coupled + orbital mixing |
| Relative cost | 1× | ~2× | ~2–4× |

---

#### Recommended Additional Tags for Non-collinear / SOC

```
GGA_COMPAT = .FALSE.    # Restores full rotational invariance of GGA (recommended)
LASPH      = .TRUE.     # Non-spherical gradient corrections (recommended)
LMAXMIX    = 4          # For d-electron systems (6 for f-electron systems)
```

> **Note**: `GGA_COMPAT = .FALSE.` and `LASPH = .TRUE.` are both recommended for non-collinear calculations to improve numerical precision of GGA.

---

#### PdCrO₂ Example

```
LSORBIT = True
MAGMOM = 36*0 \        # O: non-magnetic (3 components × 12 atoms)
    <18 Cr 3D moments> \
    18*0               # Pd: non-magnetic (3 components × 6 atoms)
```


## References

- Takatsu et al., PRB 89, 104408 (2014): Magnetic structure
