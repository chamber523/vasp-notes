# Non-collinear Magnetism in VASP

Comprehensive guide for setting up and troubleshooting non-collinear magnetic calculations with spin-orbit coupling (SOC).

## Table of Contents

1. [Theory Background](#theory-background)
2. [INCAR Setup](#incar-setup)
3. [MAGMOM Format](#magmom-format)
4. [Practical Workflow](#practical-workflow)
5. [Convergence Tips](#convergence-tips)
6. [Analysis](#analysis)
7. [Troubleshooting](#troubleshooting)

---

## Theory Background

### Collinear vs Non-collinear Magnetism

**Collinear magnetism** (`ISPIN=2`):
- Spin-up and spin-down states are antiparallel
- Magnetic moments align along a single axis (usually z)
- 2-component spinor: |↑⟩ and |↓⟩
- MAGMOM: single value per atom

**Non-collinear magnetism** (`LSORBIT=True`):
- Magnetic moments can point in any 3D direction
- Required for complex magnetic structures (spirals, canted moments)
- 4-component spinor (includes spin-orbit coupling)
- MAGMOM: 3 components (mx, my, mz) per atom

### When to Use Non-collinear Calculations

**Required for:**
- Systems with spin-orbit coupling effects
- Complex magnetic structures (120° triangular lattice, spin spirals)
- Magnetic anisotropy calculations
- Materials with heavy elements (Pd, Pt, Ir)

**Example systems:**
- PdCrO₂: 120° non-collinear ordering in triangular lattice
- Pyrochlore iridates: canted antiferromagnetism
- Topological magnets: Dzyaloshinskii-Moriya interaction

---

## INCAR Setup

### Minimal Setup

```bash
# ========== Enable Non-collinear Magnetism ==========
LSORBIT = True           # Automatically enables LNONCOLLINEAR
GGA_COMPAT = False       # Must be False for SOC

# ========== Initial Magnetic Moments (3D vectors) ==========
MAGMOM = mx1 my1 mz1  mx2 my2 mz2  mx3 my3 mz3  ...
```

**Note:** VASP automatically sets `LNONCOLLINEAR = True` when `LSORBIT = True`

### Complete Production Setup

```bash
# ========== Basic Parameters ==========
PREC = Normal
ENCUT = 400              # Adjust based on POTCAR ENMAX
EDIFF = 1e-7             # Tight convergence for magnetic systems
ALGO = Normal            # or All for difficult cases

# ========== Non-collinear Magnetism + SOC ==========
LSORBIT = True           # Enable SOC (auto-enables non-collinear)
GGA_COMPAT = False       # REQUIRED for SOC
SAXIS = 0 0 1            # Spin quantization axis (default: z)

# ========== Initial Magnetic Moments ==========
# Format: mx1 my1 mz1  mx2 my2 mz2  ...
MAGMOM = 36*0 \          # 36 O atoms (non-magnetic)
    0.250611  0.150583  0.956305 \  # Cr1
   -0.479195 -0.621832  0.621832 \  # Cr2
    ...                  # (18 Cr magnetic atoms)
    18*0                 # 18 Pd atoms (non-magnetic)

# ========== Symmetry ==========
ISYM = 0                 # Often needed to break symmetry
# or ISYM = -1           # No symmetry, but use time-reversal

# ========== Electronic Structure ==========
ISMEAR = 0               # Gaussian smearing
SIGMA = 0.05             # Small smearing width
NELM = 200               # Max electronic steps

# ========== Charge Density ==========
LCHARG = True            # Write CHGCAR for DOS/band
LWAVE = True             # Write WAVECAR (large file!)

# ========== Output ==========
LORBIT = 11              # Orbital-projected DOS
NWRITE = 2               # Medium verbosity
```

### DFT+U Addition (for 3d/4f elements)

```bash
# ========== DFT+U Parameters ==========
LDAU = True              # Enable DFT+U
LDAUTYPE = 2             # Dudarev formulation
LDAUL = -1 2 -1          # l-quantum numbers (O=off, Cr=d, Pd=off)
LDAUU = 0 4 0            # U values (eV)
LDAUJ = 0 0 0            # J values (not used in Dudarev)
LMAXMIX = 4              # Mix d-orbitals in charge density
```

**Common U values:**
| Element | U (eV) | Reference |
|---------|--------|-----------|
| Cr      | 3-4    | Cr oxides |
| Mn      | 4-5    | Mn oxides |
| Fe      | 4-5    | Fe oxides |
| Co      | 3-4    | Co oxides |
| Ni      | 5-6    | Ni oxides |
| V       | 3-4    | V oxides  |

---

## MAGMOM Format

### Basic Rules

**1. Number of components:**
- Non-collinear: 3 components per atom
- Total values: 3 × N_atoms

**2. Order must match POSCAR:**
```
POSCAR line 6: O  Cr  Pd
POSCAR line 7: 36 18  18

MAGMOM order:
  36 O atoms  (mx, my, mz) × 36
  18 Cr atoms (mx, my, mz) × 18
  18 Pd atoms (mx, my, mz) × 18
```

**3. Shorthand notation:**
```bash
# Explicit
MAGMOM = 0 0 0  0 0 0  ...  # Very long

# Shorthand for non-magnetic atoms
MAGMOM = 36*0  ...          # 36 atoms with mx=my=mz=0
                            # Total: 108 values (36×3)

# Shorthand for repeated patterns
MAGMOM = 3*0  3*5  3*0      # 3 atoms: (0,0,0), (5,0,0), (0,0,0)
```

### Calculating Initial MAGMOM

**Method 1: Literature/Experimental Data**

Example: PdCrO₂ from Takatsu et al. (PRB 89, 104408, 2014)

Each Cr layer has 3 sublattices (A, B, C) with 120° rotation:
```python
import numpy as np

# Model 4 parameters (from paper)
alpha = [31, 44, 31, 44, 31, 44]  # degrees
phi = [17, 16, 17, 16, 17, 16]    # degrees

for n, (a, p) in enumerate(zip(alpha, phi)):
    a_rad = np.radians(a)
    p_rad = np.radians(p)
    gamma = np.radians(90)  # 120° rotation angle

    # Sublattice A
    S_A = [np.cos(gamma)*np.cos(p_rad),
           np.sin(p_rad),
           np.sin(gamma)*np.cos(p_rad)]

    # Sublattice B (120° rotation)
    S_B = [-np.cos(p_rad)/2 - np.sqrt(3)/2*np.sin(p_rad),
           np.sqrt(3)/2*np.cos(p_rad) - np.sin(p_rad)/2,
           np.sin(gamma)*np.cos(p_rad)]

    # Sublattice C (240° rotation)
    S_C = [-np.cos(p_rad)/2 + np.sqrt(3)/2*np.sin(p_rad),
           -np.sqrt(3)/2*np.cos(p_rad) - np.sin(p_rad)/2,
           np.sin(gamma)*np.cos(p_rad)]

    # Scale by moment magnitude and alpha
    xi = 3.0  # Approximate Cr moment (μB)
    for S in [S_A, S_B, S_C]:
        mx = xi * np.cos(a_rad) * S[0]
        my = xi * np.cos(a_rad) * S[1]
        mz = xi * np.sin(a_rad) * S[2]
        print(f"{mx:.6f}  {my:.6f}  {mz:.6f} \\")
```

**Method 2: Simple Ferromagnet**

```bash
# All moments along z-axis
MAGMOM = 10*0 \          # 10 O (non-magnetic)
         4*0 4*0 4*5 \    # 4 Fe (5 μB along z)
         6*0              # 6 more O
```

**Method 3: Simple Antiferromagnet**

```bash
# NiO: alternating Ni spins
MAGMOM = 0 0 2  0 0 -2  0 0 0  0 0 0  # Ni↑ Ni↓ O O
```

### Non-magnetic Atoms

For atoms with no magnetic moment (O, Pd, etc.):
```bash
# Method 1: Explicit zeros
MAGMOM = 0 0 0  0 0 0  ...

# Method 2: Shorthand
MAGMOM = 36*0            # 36 atoms, 108 values total

# Method 3: Single zero per atom (VASP expands automatically)
MAGMOM = 36*0            # VASP interprets as mx=my=mz=0
```

---

## Practical Workflow

### Step 1: Prepare Structure

```bash
# Ensure POSCAR has:
# - Correct atom order
# - Atom count matching MAGMOM
# - Proper symmetry breaking if needed
```

### Step 2: Calculate MAGMOM

```bash
# Use script (e.g., generate_magmom.py)
python generate_magmom.py > MAGMOM.txt

# Copy to INCAR
# MAGMOM = <paste from MAGMOM.txt>
```

### Step 3: Validate Setup

```bash
# Check file consistency
grep "ions per type" OUTCAR
# Should match POSCAR line 7

# Check MAGMOM count
# For 72 atoms: 72 × 3 = 216 values
```

### Step 4: Run SCF

```bash
# Submit job
sbatch submit.slurm

# Monitor convergence
tail -f OSZICAR
grep "mag=" OSZICAR
```

### Step 5: Check Results

```bash
# Verify convergence
grep "reached required accuracy" OUTCAR

# Check final magnetic moments
grep "magnetization (x,y,z)" OUTCAR | head -50
```

---

## Convergence Tips

### Issue 1: Magnetic Moments Collapse to Zero

**Symptoms:**
```
OSZICAR:
N  E_diel      E_ion       E_tot       mag
1 -123.456    -234.567    -358.023    54.321
...
50 -123.499   -234.599    -358.098     0.003  # Moments collapsed!
```

**Solutions:**

**A. Increase initial moments**
```bash
MAGMOM = 36*0 \
    0.5  0.3  1.9 \    # Scaled up by 2×
    ...
```

**B. Break symmetry**
```bash
ISYM = 0               # Turn off all symmetry
LCHARG = True
NUPDOWN = -1           # Let VASP determine total moment
```

**C. Use better initial guess**
```bash
# Read from previous calculation
ISTART = 1
ICHARG = 1
```

### Issue 2: Oscillating Convergence

**Symptoms:**
```
OSZICAR:
...
45  -358.098   1.234   # Energy oscillating
46  -358.102   0.891
47  -358.099   1.189
48  -358.103   0.923
```

**Solutions:**

**A. Adjust mixing parameters**
```bash
AMIX = 0.1             # Reduce from 0.4
BMIX = 0.5             # Reduce from 3.0
AMIN = 0.01
```

**B. Try different algorithm**
```bash
ALGO = All             # More stable but slower
# or
ALGO = VeryFast        # Might help difficult cases
```

**C. Increase mixing history**
```bash
MAXMIX = 80            # Increase from 40
```

### Issue 3: Slow Convergence

**Symptoms:**
- Takes > 100 electronic steps per ionic step
- Small energy changes

**Solutions:**

**A. Skip initial diagonalization**
```bash
NELMDL = -10           # Skip first 10 steps
```

**B. Optimize parallelization**
```bash
NPAR = 8               # Test different values
KPAR = 2               # Parallelize k-points
NCORE = 4              # Cores per band
```

**C. Adjust precision dynamically**
```bash
# Start with lower precision
PREC = Normal
ALGO = Fast

# Then tighten if needed
EDIFF = 1e-5           # Initial run
EDIFF = 1e-7           # Final run
```

---

## Analysis

### Extract Magnetic Moments

**From OUTCAR:**

```bash
# Total magnetization
grep "number of electron" OUTCAR | tail -1

# Per-atom moments (3D vectors)
grep "magnetization (x,y,z)" OUTCAR -A 80

# Example output:
#  ion      S_x      S_y      S_z
#    1    0.251    0.151    0.956
#    2   -0.480   -0.622    0.622
```

**Using Python (pymatgen):**

```python
from pymatgen.io.vasp import Outcar
import numpy as np

# Read OUTCAR
outcar = Outcar("OUTCAR")

# Get magnetization (list of [mx, my, mz] for each atom)
magmom = outcar.magnetization

# Calculate magnitudes
magnitudes = [np.linalg.norm(m) for m in magmom]

# Print magnetic atoms (|m| > 0.1 μB)
for i, (mag, m) in enumerate(zip(magnitudes, magmom)):
    if mag > 0.1:
        print(f"Atom {i+1:3d}: |m| = {mag:6.3f} μB, "
              f"({m[0]:7.3f}, {m[1]:7.3f}, {m[2]:7.3f})")
```

### Visualize with VESTA

**Generate MCIF file:**

```python
# Use generate_mcif.py script
python generate_mcif.py POSCAR

# Output: magnetic_structure.mcif
```

**Open in VESTA:**
1. File → Open → magnetic_structure.mcif
2. Edit → Atoms → Show magnetic moments
3. Adjust arrow size: Objects → Properties → Magnetic moments

### Check Convergence History

```python
import matplotlib.pyplot as plt
import re

# Parse OSZICAR
energies, mags = [], []
with open("OSZICAR", "r") as f:
    for line in f:
        if line.strip() and not line.startswith("N"):
            parts = line.split()
            energies.append(float(parts[2]))  # E0
            if "mag=" in line:
                mags.append(float(parts[-1]))

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

ax1.plot(energies)
ax1.set_ylabel("Energy (eV)")
ax1.set_xlabel("Electronic Step")
ax1.grid(True)

ax2.plot(mags)
ax2.set_ylabel("Total Magnetization (μB)")
ax2.set_xlabel("Electronic Step")
ax2.grid(True)

plt.tight_layout()
plt.savefig("convergence.png", dpi=300)
```

---

## Troubleshooting

### Problem: "VERY BAD NEWS! internal error in subroutine PRICEL"

**Cause:** Symmetry issues with magnetic structure

**Solution:**
```bash
ISYM = 0               # Disable symmetry completely
```

### Problem: "WARNING: random wavefunctions but no delay for mixing"

**Cause:** Starting from scratch with MAGMOM

**Solution:** (Usually harmless, but can add)
```bash
NELMDL = -5            # Delay mixing for 5 steps
```

### Problem: Magnetic moments differ significantly from initial MAGMOM

**Is this bad?** Not necessarily!
- Initial MAGMOM is just a guess
- VASP finds self-consistent solution
- Check if result makes physical sense

**When to worry:**
- All moments → 0 (unintended non-magnetic solution)
- Wrong magnetic order (FM instead of AFM)

**Solutions:**
- Try different initial configurations
- Increase MAGMOM magnitude
- Use NUPDOWN to constrain total moment

### Problem: Different runs converge to different magnetic states

**Cause:** Multiple local minima in energy landscape

**Solution:**
```bash
# Method 1: Fix total moment
NUPDOWN = 10           # Constrain N_up - N_down

# Method 2: Try multiple initial configurations
MAGMOM = ...           # Config A (run 1)
MAGMOM = ...           # Config B (run 2)
# Compare final energies
```

### Problem: SOC calculation much slower than collinear

**Expected behavior:**
- SOC calculations are 2-4× slower
- Wavefunction size doubles (4 components vs 2)

**Optimization:**
```bash
# Reduce k-points if possible
KPOINTS: 6 6 4  →  4 4 2

# Optimize parallelization
NPAR = 16              # Larger NPAR for SOC
NCORE = 2              # Smaller NCORE

# Use less strict convergence initially
EDIFF = 1e-5           # Quick test run
EDIFF = 1e-7           # Production run
```

---

## Advanced Topics

### Constrained Magnetic Calculations

Fix magnetic moment direction while optimizing magnitude:

```bash
I_CONSTRAINED_M = 2    # Constrain magnetic moments
LAMBDA = 10            # Penalty parameter
M_CONSTR = mx my mz    # Constrain direction per atom (in POSCAR order)
```

### Non-collinear + Hybrid Functionals

```bash
LHFCALC = True         # HSE06
HFSCREEN = 0.2
LSORBIT = True         # SOC
ALGO = All             # Required for hybrid+SOC
```

**Warning:** Extremely expensive! Consider PBE0 instead of HSE06.

### Magnetic Anisotropy Energy (MAE)

Calculate energy difference between different magnetization directions:

```bash
# Run 1: Moments along z
SAXIS = 0 0 1
MAGMOM = 0 0 5  ...    # All moments along z

# Run 2: Moments along x
SAXIS = 1 0 0
MAGMOM = 5 0 0  ...    # All moments along x

# MAE = E(x) - E(z)
```

---

## References

1. **VASP Manual:**
   - [Non-collinear calculations](https://www.vasp.at/wiki/index.php/Non-collinear_calculations)
   - [LSORBIT](https://www.vasp.at/wiki/index.php/LSORBIT)
   - [MAGMOM](https://www.vasp.at/wiki/index.php/MAGMOM)

2. **Literature:**
   - Hobbs et al., PRB 62, 11556 (2000): "Fully unconstrained noncollinear magnetism"
   - Takatsu et al., PRB 89, 104408 (2014): PdCrO₂ magnetic structure

3. **Tools:**
   - [pymatgen](https://pymatgen.org/): Analysis and structure manipulation
   - [VESTA](https://jp-minerals.org/vesta/): Visualization with magnetic moments

---

**Last updated:** 2025-01-30
**VASP version:** 6.x
