# Band Alignment Calculation

Work function and Schottky barrier height calculations for metal-semiconductor interfaces.

## Example: PbSe/Al2Se3/Al Interface

```
band-alignment-calculation/
└── example-PbSe_Al2Se3_Al/           # PbSe(111)||Al2Se3(001)||Al(111)
    ├── INCAR                          # VASP parameters
    ├── KPOINTS                        # 8×8×1 mesh
    ├── POSCAR                         # Interface structure (OgreInterface)
    ├── submit.slurm                   # Job script
    ├── MACROSCOPIC_AVERAGE_*.dat      # Potential data (VASPKIT 427)
    ├── SBH.ipynb                      # Analysis notebook
    └── pot_combined.png               # Potential visualization
```

## Workflow

### 1. Structure Preparation

Build interface using **OgreInterface**:
- Match lattice: PbSe(111) || Al2Se3(001) || Al(111)
- Add vacuum layer (~20 Å)
- Interface naming: `Material1_Material2_Material3_hkl1_hkl2_hkl3_vacuumThickness`

### 2. VASP Calculation

```bash
cd example-PbSe_Al2Se3_Al/
sbatch submit.slurm
```

**Key INCAR settings:**
```bash
# Basic
PREC = Normal
ENCUT = 400
EDIFF = 1e-5
SIGMA = 0.05
ISMEAR = 0

# SOC
LSORBIT = True
GGA_COMPAT = False
MAGMOM = 1944*0.0          # Total number of atoms × 3

# van der Waals
IVDW = 20                  # Tkatchenko-Scheffler
LVDW_EWALD = True

# Slab dipole correction (CRITICAL!)
IDIPOL = 3                 # z-direction
LDIPOL = True
DIPOL = 0.51925 0.49223 0.42848    # Center of mass from POSCAR

# Output potential
LVHAR = True               # Write LOCPOT for potential analysis

# Single-point calculation
NSW = 1
IBRION = -1
LCHARG = False             # Save memory
LWAVE = False
```

### 3. Extract Macroscopic Averaged Potential

Use **VASPKIT 427** after VASP completes:

```bash
vaspkit
# Select: 427 (Planar/Macroscopic averaged potential)

# Prompts:
# 1. Select direction: 3 (z-direction for slab)
# 2. Period length: ~layer distance (e.g., 3.5 Å)
# 3. Iteration number: 3-5 (for smoother curve)
```

**Output files:**
- `PLANAR_AVERAGE.dat`: Raw planar averaged potential
- `MACROSCOPIC_AVERAGE_<period>.dat`: Smoothed macroscopic average

### 4. Analysis

Open `SBH.ipynb` to:
1. Load multiple `MACROSCOPIC_AVERAGE_*.dat` files
2. Plot potential along z-direction
3. Identify plateaus in metal and semiconductor regions
4. Calculate work function differences
5. Determine Schottky barrier height

## Key Concepts

### Work Function

Energy required to move an electron from inside the solid to vacuum.

**Work function = Vacuum level - Fermi level**

### Dipole Correction

For asymmetric slabs (different terminations on two sides), the potential in vacuum becomes a slope instead of a flat plateau. Dipole correction adds a step function to make the vacuum potential flat on both sides.

```bash
LDIPOL = True           # Enable dipole correction
IDIPOL = 3              # Direction (3 = z-axis for slab)
DIPOL = x y z           # Dipole center (typically center of mass)
```

### Macroscopic Averaged Potential

VASPKIT calculates macroscopic-averaged potential by converging the planar averaged potential with a moving average:

- **Planar average**: z = average(x, y) for every z
- **Macroscopic average**: z = average(x, y, z ± Δz') for every z, where 2Δz' is the layer distance

This removes atomic-scale oscillations to produce a smooth potential curve.

## References

- https://vaspkit.com/tutorials.html
