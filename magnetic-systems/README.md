# Magnetic Systems Workflow

Non-collinear magnetic calculations with SOC for PdCrO₂-based systems.

## Directory Structure

```
magnetic-systems/
├── model/                        # Structure preparation
│   ├── shift_z.py                 # Reorient cell (Cr to bottom)
│   ├── make_supercell.py          # Generate supercell with transformation matrix
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
cd model/

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
cd calculation/scf_noncollinear/
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

```bash
LSORBIT = True         # SOC + non-collinear
MAGMOM = 36*0 \        # O: non-magnetic
    <18 Cr 3D moments> \
    18*0               # Pd: non-magnetic
```

### DFT+U

```bash
LDAU = True
LDAUTYPE = 2
LDAUU = 0 4 0          # U = 4.0 eV for Cr
```

### Band Unfolding

- Transformation matrix: `[[2,1,0],[1,2,0],[0,0,2]]`
- Visualization: vaspvis
- Requirement: `LWAVE = True` in band_unfold/INCAR

## References

- Takatsu et al., PRB 89, 104408 (2014): Magnetic structure
