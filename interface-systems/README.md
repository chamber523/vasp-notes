# Interface Systems: Work Function and Schottky Barrier

Metal-semiconductor interface slab calculations for work function and Schottky barrier height.

## Example: PbSe/Al2Se3/Al Interface

```
interface-systems/
└── PbSe_Al2Se3_Al_111_001_111_20_8_8/    # PbSe(111)||Al2Se3(001)||Al(111)
    ├── INCAR                              # VASP parameters
    ├── KPOINTS                            # 8×8×1 mesh
    ├── POSCAR                             # Interface structure (OgreInterface)
    ├── submit.slurm                       # Job script
    ├── MACROSCOPIC_AVERAGE_*.dat          # Potential data (VASPKIT 427)
    ├── SBH.ipynb                          # Analysis notebook
    └── pot_combined.png                   # Potential visualization
```

## Key Concepts

### Work Function (功函数)

Energy required to move an electron from inside the solid to vacuum (analogous to photoelectric effect).

**Work function = Vacuum level - Fermi level**

### Vacuum Level (真空能级)

Energy of a free electron in vacuum far from the surface.

### Dipole Correction (偶极校正)

**Problem:** For asymmetric slabs (e.g., different surfaces on two sides), the potential in vacuum becomes a **slope** instead of a flat plateau.

**Solution:** Dipole correction adds a step function to make the vacuum potential flat on both sides:
```
Without correction:    With correction:
    /                      ___
   /  (slope)             |   | (step)
  /                       |___|
```

INCAR parameters:
```bash
LDIPOL = True           # Enable dipole correction
IDIPOL = 3              # Direction (3 = z-axis for slab)
DIPOL = x y z           # Dipole center (typically center of mass)
```

## Workflow

### 1. Structure Preparation

Build interface using **OgreInterface**:
- Match lattice: PbSe(111) || Al2Se3(001) || Al(111)
- Add vacuum layer (~20 Å)
- Interface naming: `PbSe_Al2Se3_Al_111_001_111_20_8_8`
  - Format: `Material1_Material2_Material3_hkl1_hkl2_hkl3_vacuumThickness_latticeParam1_latticeParam2`

### 2. VASP Calculation

```bash
cd PbSe_Al2Se3_Al_111_001_111_20_8_8/
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
#    (Should be close to but smaller than actual layer spacing)
# 3. Iteration number: 3-5 (for smoother curve)
```

**Output files:**
- `PLANAR_AVERAGE.dat`: Raw planar averaged potential
- `MACROSCOPIC_AVERAGE_<period>.dat`: Smoothed macroscopic average

**What is macroscopic averaging?**

Removes atomic-scale oscillations in potential:
```
Planar average:      Macroscopic average:
  /\/\/\               _____
 /      \/\           /     \___
         \/\         /          \___
(atomic oscillation)    (smooth curve)
```

Formula: For each z, average over z ± Δz' where 2Δz' ≈ layer distance

### 4. Analysis

Open `SBH.ipynb` to:
1. Load multiple `MACROSCOPIC_AVERAGE_*.dat` files
2. Plot potential along z-direction
3. Identify plateaus in metal and semiconductor regions
4. Calculate work function differences
5. Determine Schottky barrier height

**Typical plot:**
```
Potential (eV)
   |     Metal    |  Semiconductor  |  Metal
   |    ______    |     _______     |  ______
   |   /      \___|____/       \____|_/
   |
   +-----------------------------------------> z (Å)
        Barrier   Interface   Barrier
```

## Important Notes

### Why Dipole Correction?

For **asymmetric** slabs (different terminations on two sides):
- Without correction: vacuum potential slopes → cannot define vacuum level
- With correction: vacuum potential is flat on both sides → well-defined vacuum level

For **symmetric** slabs: dipole correction may not be necessary but doesn't hurt.

### Vacuum Layer Size

- Minimum: 15 Å
- Recommended: 20+ Å
- Check: Potential should plateau in vacuum region (flat line)

If potential doesn't plateau → increase vacuum layer.

### K-point Convergence

```
KPOINTS
Automatic
0
Gamma
8 8 1
```
- Dense in xy (parallel to interface): 8×8
- Sparse in z (perpendicular): 1
- Larger systems may need denser mesh

## References

- Van de Walle & Martin, PRB 35, 8154 (1987): Macroscopic averaging method
- VASPKIT manual: Tool 427 for potential averaging
