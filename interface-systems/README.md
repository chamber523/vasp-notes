# Interface Systems: Schottky Barrier Height Calculations

Metal-semiconductor interface calculations with SOC and van der Waals corrections.

## Directory Structure

```
interface-systems/
└── PbSe_Al2Se3_Al_111_001_111_20_8_8/    # Example: PbSe/Al2Se3/Al interface
    ├── INCAR                              # VASP input parameters
    ├── KPOINTS                            # K-point mesh
    ├── POSCAR                             # Interface structure
    ├── POTCAR                             # Pseudopotentials
    ├── submit.slurm                       # Job submission script
    ├── LOCPOT                             # Local potential for analysis
    ├── MACROSCOPIC_AVERAGE_*.dat          # Macroscopic averaged potentials
    ├── SBH.ipynb                          # Schottky barrier height analysis
    └── pot_combined.png                   # Potential plot
```

## Workflow

### 1. Structure Preparation

- Build metal-semiconductor interface slab
- Add sufficient vacuum layer (~20 Å)
- Ensure symmetric termination to avoid dipole issues

### 2. VASP Calculation

```bash
cd PbSe_Al2Se3_Al_111_001_111_20_8_8/
sbatch submit.slurm
```

Key INCAR parameters:
```bash
# SOC
LSORBIT = True
GGA_COMPAT = False

# van der Waals correction
IVDW = 20              # Tkatchenko-Scheffler
LVDW_EWALD = True

# Slab dipole correction
IDIPOL = 3             # z-direction
LDIPOL = True
DIPOL = <x y z>        # Center of mass

# Potential output
LVHAR = True           # Output LOCPOT
```

### 3. Analysis

Extract macroscopic averaged potential and calculate Schottky barrier height using `SBH.ipynb`.

**Analysis steps:**
1. Load MACROSCOPIC_AVERAGE_*.dat files
2. Identify potential plateaus in metal and semiconductor regions
3. Calculate work function offset
4. Determine barrier height

## Key Parameters

### Dipole Correction

For asymmetric slabs, dipole correction is essential:
```bash
IDIPOL = 3             # Dipole along z
LDIPOL = True          # Enable correction
DIPOL = x y z          # Dipole center (typically center of mass)
```

### Vacuum Layer

- Minimum: 15 Å
- Recommended: 20+ Å
- Check LOCPOT to ensure potential plateaus are reached

### K-point Mesh

```
KPOINTS
Automatic
0
Gamma
8 8 1    # Dense in xy, sparse along z (slab direction)
```

## References

- Van de Walle & Martin, PRB 35, 8154 (1987): Band alignment methodology
