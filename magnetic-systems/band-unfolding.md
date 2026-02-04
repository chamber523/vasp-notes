# Band Structure Unfolding for Supercells

Complete guide for unfolding band structures from supercell calculations back to the primitive cell Brillouin zone.

## Table of Contents

1. [Why Unfold?](#why-unfold)
2. [Theory Overview](#theory-overview)
3. [VASP Setup](#vasp-setup)
4. [Unfolding with easyunfold](#unfolding-with-easyunfold)
5. [Alternative: BandUP](#alternative-bandup)
6. [Practical Examples](#practical-examples)
7. [Troubleshooting](#troubleshooting)

---

## Why Unfold?

### The Problem

When calculating band structures for **supercells** (magnetic ordering, defects, alloys), the resulting band structure is "folded":

**Primitive cell (1×1×1):**
```
Simple band structure in primitive Brillouin zone
Clear band dispersion, easy to interpret
```

**Supercell (2×2×2):**
```
Band structure is folded 8 times (2³)
Many overlapping bands (backfolded bands)
Difficult to identify primitive cell features
```

**Visual comparison:**
```
Primitive cell bands:     Supercell bands (folded):
      E                         E
      |  ╱─╲                    |  ╱─┬─╲
      | ╱   ╲                   | ╱ ┌┴┐ ╲
      |╱     ╲                  |╱  │││  ╲
      └────────> k              └───┴┴┴───> k
      Clear dispersion          Messy, hard to interpret
```

### The Solution: Band Unfolding

**Unfolding** maps supercell bands back to the primitive cell BZ:
- Identifies which supercell bands correspond to primitive cell states
- Calculates spectral weights (band intensity)
- Produces clean band structure comparable to primitive cell

**Use cases:**
- Magnetic supercells (like PdCrO₂ 2×2×2 supercell)
- Defect calculations
- Alloy/doping calculations (SQS structures)
- Comparing supercell results with pristine system

---

## Theory Overview

### Brillouin Zone Folding

**Primitive cell:**
- Lattice vectors: **a**, **b**, **c**
- Reciprocal lattice vectors: **b₁**, **b₂**, **b₃**
- BZ size: (2π/a) × (2π/b) × (2π/c)

**Supercell (2×2×2):**
- Lattice vectors: 2**a**, 2**b**, 2**c**
- Reciprocal lattice: **b₁**/2, **b₂**/2, **b₃**/2
- BZ size: (π/a) × (π/b) × (π/c)  ← Half in each direction!

**Supercell BZ is 8× smaller** → primitive cell BZ is folded 8 times into supercell BZ

### Transformation Matrix

**Direct space transformation:**
```
Supercell = M × Primitive cell

M = [[a₁₁, a₁₂, a₁₃],
     [a₂₁, a₂₂, a₂₃],
     [a₃₁, a₃₂, a₃₃]]
```

**Example: PdCrO₂ magnetic supercell**
```python
M = [[2, 1, 0],    # âₛᶜ = 2â + b̂
     [1, 2, 0],    # b̂ₛᶜ = â + 2b̂
     [0, 0, 2]]    # ĉₛᶜ = 2ĉ
```

**Reciprocal space transformation:**
```
k_supercell = M^T × k_primitive + G

where G is reciprocal lattice vector
```

### Spectral Weight

**Unfolding assigns spectral weight to each band:**
- Weight = 1: Band is pure primitive cell character
- Weight = 0: Band is artifact of supercell
- 0 < Weight < 1: Partial character (due to symmetry breaking)

**Calculation:**
```
P(k, n) = |⟨ψₙᵏ_SC | φₘᴷ_PC⟩|²

where:
  ψₙᵏ_SC = Supercell Bloch function
  φₘᴷ_PC = Primitive cell Bloch function
```

---

## VASP Setup

### Critical: LWAVE = True

**INCAR requirement:**
```bash
LWAVE = True           # REQUIRED: Write WAVECAR
LCHARG = True          # Recommended: Write CHGCAR
```

**Why LWAVE = True?**
- WAVECAR contains plane-wave coefficients
- Unfolding tools need these to calculate projections
- File is large (>1 GB for 72 atoms) but necessary

**DO NOT:**
```bash
LWAVE = False          # ✗ Unfolding will fail
LWAVE = .FALSE.        # ✗ Same
# (commented) LWAVE    # ✗ Defaults to False after SCF
```

### Band Structure Calculation

**Step 1: SCF Calculation (charge density)**

```bash
cd 01_scf/

INCAR:
  ISTART = 0
  ICHARG = 2
  LCHARG = True        # Write CHGCAR
  LWAVE = True         # Can be False for SCF (saves space)

KPOINTS:
  Automatic
  0
  Monkhorst-Pack
  8 8 1              # Dense mesh for charge density
```

**Step 2: Band Structure (for unfolding)**

```bash
cd 03_band_unfold/

# Copy converged charge density
cp ../01_scf/CHGCAR .

INCAR:
  ISTART = 1           # Read WAVECAR
  ICHARG = 11          # Read CHGCAR (non-self-consistent)
  LWAVE = True         # ✓ CRITICAL for unfolding
  LCHARG = False       # Save space
  NSW = 0              # Static
  IBRION = -1

KPOINTS:
  Line mode            # High-symmetry path
  40                   # Points per segment
  Reciprocal
  0.0000  0.0000  0.0   ! Γ
  0.5000  0.0000  0.0   ! M

  0.5000  0.0000  0.0   ! M
  0.3333  0.3333  0.0   ! K

  0.3333  0.3333  0.0   ! K
  0.0000  0.0000  0.0   ! Γ
```

**Important notes:**
- Use **supercell BZ** k-points for band calculation
- Will unfold to **primitive BZ** later
- Choose high-symmetry path appropriate for primitive cell structure

### High-Symmetry Paths

**For hexagonal/trigonal systems (like PdCrO₂):**

**Primitive cell path:** Γ-M-K-Γ-A-L-H-A

**Supercell path (folded BZ):**
- Same labels but different coordinates
- Use primitive cell symmetry, not supercell

**Example coordinates (fractional):**
```
Γ = (0.0000, 0.0000, 0.0000)
M = (0.5000, 0.0000, 0.0000)
K = (0.3333, 0.3333, 0.0000)
A = (0.0000, 0.0000, 0.5000)
L = (0.5000, 0.0000, 0.5000)
H = (0.3333, 0.3333, 0.5000)
```

**Tool to find high-symmetry points:**
- [SeeK-path](https://www.materialscloud.org/work/tools/seekpath)
- [pymatgen](https://pymatgen.org/): `HighSymmKpath`

---

## Unfolding with easyunfold

### Installation

```bash
# Using pip
pip install easyunfold

# Or from source
git clone https://github.com/SMTG-UCL/easyunfold.git
cd easyunfold
pip install .
```

### Complete Workflow

**Directory structure:**
```
calculation/
├── 01_scf/              # Converged SCF
│   └── CHGCAR
└── 03_band_unfold/      # Band + unfolding
    ├── POSCAR
    ├── POTCAR
    ├── KPOINTS
    ├── INCAR           (LWAVE = True)
    ├── CHGCAR          (from SCF)
    └── [After VASP runs...]
        ├── WAVECAR
        ├── EIGENVAL
        └── vasprun.xml
```

**Step 1: Generate k-points**

```bash
cd 03_band_unfold/

# Specify transformation matrix
easyunfold generate --matrix "[[2,1,0],[1,2,0],[0,0,2]]" \
                    --kpoint-path "GMKG"

# Output:
# - KPOINTS_band: High-symmetry path in supercell BZ
# - easyunfold.json: Metadata for unfolding
```

**Arguments:**
- `--matrix`: Transformation from primitive → supercell
- `--kpoint-path`: High-symmetry labels (Γ=G, M, K, etc.)
- `--nkpts`: Points per segment (default: 40)

**Step 2: Run VASP**

```bash
# Use generated KPOINTS_band or your own
cp KPOINTS_band KPOINTS

# Run VASP
sbatch ../submit_template.slurm

# Wait for convergence...
```

**Step 3: Unfold bands**

```bash
# Perform unfolding
easyunfold unfold calculate

# Output:
# - unfold_result.json: Unfolded eigenvalues and weights
# - unfold.csv: Human-readable format
```

**Step 4: Plot**

```bash
# Generate band structure plot
easyunfold unfold plot

# Output:
# - band_structure.png: Unfolded bands with weights
```

**Customization:**
```bash
# Plot with energy limits
easyunfold unfold plot --emin -5 --emax 5 --fermi 0

# Plot with custom colormap
easyunfold unfold plot --cmap viridis

# Export data for custom plotting
easyunfold unfold plot --no-plot --output bands.dat
```

### Advanced Options

**Spin-polarized calculations:**
```bash
# Unfold both spins
easyunfold unfold calculate --spin both

# Plot separately
easyunfold unfold plot --spin up
easyunfold unfold plot --spin down
```

**Non-collinear/SOC:**
```bash
# Automatically detected from OUTCAR
easyunfold unfold calculate

# If needed, specify explicitly
easyunfold unfold calculate --soc
```

**Different WAVECAR file:**
```bash
easyunfold unfold calculate --wavecar path/to/WAVECAR
```

---

## Alternative: BandUP

### Installation

```bash
# Download from GitHub
git clone https://github.com/band-unfolding/bandup.git
cd bandup

# Compile (requires Fortran compiler)
make

# Add to PATH
export PATH=$PATH:/path/to/bandup/bin
```

### Workflow

**Step 1: Prepare input**

Create `bandup.in`:
```
&input_from_user
  PW_code = 'VASP'
  input_file_prim_cell = 'POSCAR_primitive'
  input_file_super_cell = 'POSCAR'
  transformation_matrix(1,:) = 2, 1, 0
  transformation_matrix(2,:) = 1, 2, 0
  transformation_matrix(3,:) = 0, 0, 2
  kpoints_file = 'KPOINTS_prim_cell'
/
```

**Step 2: Run BandUP**

```bash
# Generate unfolding
bandup

# Output:
# - unfolded_EBS.dat: Energy bands with weights
```

**Step 3: Plot**

```python
# Use provided plotting script
python plot_unfolded_bands.py unfolded_EBS.dat
```

---

## Practical Examples

### Example 1: PdCrO₂ Magnetic Supercell

**System:**
- Primitive cell: 3 atoms (O, Cr, Pd)
- Supercell: 72 atoms (2×2×2 magnetic supercell)
- Transformation: [[2,1,0],[1,2,0],[0,0,2]]

**INCAR (03_band_unfold/INCAR):**
```bash
# Critical settings
LWAVE = True           # ← Must have!
ICHARG = 11            # Read CHGCAR
LSORBIT = True         # SOC enabled

# Standard band settings
ISTART = 1
IBRION = -1
NSW = 0
LORBIT = 11
NEDOS = 3001
```

**KPOINTS (Γ-M-K-Γ path):**
```
Line mode
40
Reciprocal
0.0000  0.0000  0.0   Γ
0.5000  0.0000  0.0   M

0.5000  0.0000  0.0   M
0.3333  0.3333  0.0   K

0.3333  0.3333  0.0   K
0.0000  0.0000  0.0   Γ
```

**Unfolding:**
```bash
# After VASP completes
cd 03_band_unfold/

# Method A: easyunfold
easyunfold unfold calculate --matrix "[[2,1,0],[1,2,0],[0,0,2]]"
easyunfold unfold plot --emin -3 --emax 3

# Method B: Manual using pymatgen
python unfold_bands.py
```

### Example 2: 2×2×1 Supercell (Surface Slab)

**System:**
- Surface slab with 2×2 lateral expansion
- No c-axis expansion
- Transformation: [[2,0,0],[0,2,0],[0,0,1]]

**Command:**
```bash
easyunfold generate --matrix "[[2,0,0],[0,2,0],[0,0,1]]" \
                    --kpoint-path "GMKG"
```

### Example 3: √3×√3 Supercell (60° rotation)

**System:**
- Rotated honeycomb supercell
- Common for graphene, h-BN
- Transformation: [[2,-1,0],[1,1,0],[0,0,1]]

**Command:**
```bash
easyunfold generate --matrix "[[2,-1,0],[1,1,0],[0,0,1]]" \
                    --kpoint-path "GMKG"
```

---

## Troubleshooting

### Problem: "WAVECAR not found" or "Cannot read WAVECAR"

**Solutions:**
```bash
# 1. Check LWAVE in INCAR
grep LWAVE INCAR
# Should be: LWAVE = True or LWAVE = .TRUE.

# 2. Check WAVECAR exists and is not empty
ls -lh WAVECAR
# Should be >100 MB for typical systems

# 3. Check VASP completed successfully
grep "reached required accuracy" OUTCAR
```

### Problem: Unfolded bands look wrong (no clear dispersion)

**Possible causes:**

**A. Wrong transformation matrix**
```bash
# Verify transformation
python -c "
from pymatgen.core import Structure
prim = Structure.from_file('POSCAR_primitive')
super = Structure.from_file('POSCAR_supercell')
print('Primitive lattice:', prim.lattice.matrix)
print('Supercell lattice:', super.lattice.matrix)
# Check if Supercell ≈ Matrix × Primitive
"
```

**B. Wrong k-point path**
- Using supercell symmetry instead of primitive cell
- Solution: Use primitive cell high-symmetry points

**C. Too few k-points**
```bash
# Increase points per segment
KPOINTS:
  40  →  80           # More points for smoother bands
```

### Problem: All spectral weights are very low (<0.5)

**Expected in some cases:**
- Strong symmetry breaking (magnetic ordering, defects)
- Means supercell bands deviate from primitive cell

**When to worry:**
- All weights uniformly low across all bands
- Suggests wrong transformation matrix or setup error

**Solutions:**
```bash
# 1. Verify supercell was built correctly
python make_supercell.py --verify

# 2. Check if magnetic order breaks translation symmetry severely
# (This can reduce spectral weights legitimately)

# 3. Try unfolding from collinear calculation first (test case)
```

### Problem: easyunfold fails with "Incompatible k-points"

**Cause:** K-point path doesn't match transformation

**Solution:**
```bash
# Let easyunfold generate k-points
easyunfold generate --matrix "..." --kpoint-path "GMKG"

# Use generated KPOINTS_band
cp KPOINTS_band KPOINTS
```

### Problem: Unfolding is very slow

**Expected:**
- Large WAVECAR files (>1 GB)
- Many k-points (>100)
- Can take 10-30 minutes

**Speed up:**
```bash
# 1. Reduce k-points (if resolution sufficient)
KPOINTS: 40/segment → 20/segment

# 2. Unfold only specific energy range
easyunfold unfold calculate --emin -5 --emax 5

# 3. Use fewer bands
NBANDS = 200  →  NBANDS = 100  (in INCAR)
```

---

## Custom Plotting

### Python Script (Matplotlib)

```python
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.electronic_structure.core import Spin
from easyunfold.unfold import UnfoldKSet

# Load unfolded data
unfold = UnfoldKSet.from_file("unfold_result.json")

# Get data
kpoints = unfold.kpoint_distances
energies = unfold.energies
weights = unfold.weights

# Plot
fig, ax = plt.subplots(figsize=(6, 8))

# Scatter plot with weights as color intensity
for i in range(energies.shape[1]):  # Loop over bands
    scatter = ax.scatter(kpoints, energies[:, i],
                        c=weights[:, i], s=20,
                        cmap='hot_r', vmin=0, vmax=1,
                        alpha=0.8)

# Add Fermi level
ax.axhline(0, color='gray', linestyle='--', linewidth=1)

# Add high-symmetry point labels
labels = ['Γ', 'M', 'K', 'Γ']
positions = [0.0, 1.0, 2.0, 3.0]  # Adjust based on path
for label, pos in zip(labels, positions):
    ax.axvline(pos, color='black', linewidth=0.5)
    ax.text(pos, ax.get_ylim()[1], label,
            ha='center', va='bottom', fontsize=12)

# Formatting
ax.set_ylabel("Energy (eV)", fontsize=12)
ax.set_ylim(-3, 3)
ax.set_xlim(kpoints.min(), kpoints.max())
ax.set_xticks([])

# Colorbar
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label("Spectral Weight", fontsize=12)

plt.tight_layout()
plt.savefig("unfolded_bands.png", dpi=300)
```

### Compare with Primitive Cell

```python
# Load both primitive and unfolded supercell bands
from pymatgen.io.vasp import BSVasprun

# Primitive cell band structure
bs_prim = BSVasprun("primitive/vasprun.xml").get_band_structure()

# Unfolded supercell
unfold = UnfoldKSet.from_file("supercell/unfold_result.json")

# Plot comparison
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Primitive cell (left)
bs_prim.plot(ax=ax1, show=False)
ax1.set_title("Primitive Cell")

# Unfolded supercell (right)
unfold.plot(ax=ax2)
ax2.set_title("Unfolded Supercell")

plt.tight_layout()
plt.savefig("comparison.png", dpi=300)
```

---

## Best Practices

1. **Always verify transformation matrix**
   - Check that supercell = M × primitive cell
   - Use `pymatgen` or `ASE` to verify

2. **Use sufficient k-points**
   - Minimum 20-40 points per segment
   - More for detailed features

3. **Check WAVECAR file**
   - Should be large (>100 MB typically)
   - Verify it's not corrupted: `strings WAVECAR | head`

4. **Compare with primitive cell**
   - If possible, run primitive cell calculation
   - Unfolded bands should match closely

5. **Document transformation**
   - Save transformation matrix in README
   - Note any special symmetry considerations

---

## References

1. **easyunfold:**
   - [GitHub](https://github.com/SMTG-UCL/easyunfold)
   - [Documentation](https://smtg-ucl.github.io/easyunfold/)
   - [Paper](https://doi.org/10.21105/joss.05974)

2. **BandUP:**
   - [GitHub](https://github.com/band-unfolding/bandup)
   - [Paper: Medeiros et al., PRB 89, 041407(R) (2014)](https://doi.org/10.1103/PhysRevB.89.041407)
   - [Paper: Medeiros et al., PRB 91, 041116(R) (2015)](https://doi.org/10.1103/PhysRevB.91.041116)

3. **Theory:**
   - Popescu & Zunger, PRB 85, 085201 (2012): "Extracting E versus k effective band structure"
   - Allen et al., PRB 87, 085322 (2013): "Band unfolding formalism"

---

**Last updated:** 2025-01-30
**Tested with:**
- easyunfold 0.3.3
- VASP 6.x
- pymatgen 2024.x
