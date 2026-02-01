# Calculation Directory

DFT calculations for PdCrO₂ non-collinear magnetic system using VASP.

## Directory Structure

```
calculation/
├── 01_scf/              # Self-consistent field calculation
│   ├── INCAR            # Input parameters (detailed comments)
│   ├── KPOINTS          # 8×8×1 k-point mesh
│   ├── POSCAR           # From ../model/03_POSCAR_supercell_6x
│   └── POTCAR           # Pseudopotentials (O, Cr, Pd)
├── 02_dos/              # Density of states calculation
│   ├── INCAR            # DOS-specific parameters
│   ├── KPOINTS          # 12×12×1 denser mesh
│   ├── POSCAR           # Same as SCF
│   ├── POTCAR           # Same as SCF
│   └── CHGCAR           # Copied from 01_scf/
├── 03_band_unfold/      # Band structure with unfolding
│   ├── INCAR            # Band-specific parameters
│   ├── KPOINTS          # High-symmetry path
│   ├── POSCAR           # Same as SCF
│   ├── POTCAR           # Same as SCF
│   └── CHGCAR           # Copied from 01_scf/
└── README.md            # This file
```

## Calculation Workflow

### Step 1: Prepare Files

```bash
# Copy POSCAR from model directory
cp ../model/03_POSCAR_supercell_6x 01_scf/POSCAR
cp ../model/03_POSCAR_supercell_6x 02_dos/POSCAR
cp ../model/03_POSCAR_supercell_6x 03_band_unfold/POSCAR

# Generate POTCAR (requires VASP pseudopotential library)
# Order must match POSCAR: O, Cr, Pd
cd 01_scf
cat $VASP_PP_PATH/potpaw_PBE/O/POTCAR \
    $VASP_PP_PATH/potpaw_PBE/Cr/POTCAR \
    $VASP_PP_PATH/potpaw_PBE/Pd/POTCAR > POTCAR
cd ..

# Copy POTCAR to other directories
cp 01_scf/POTCAR 02_dos/
cp 01_scf/POTCAR 03_band_unfold/
```

### Step 2: Run SCF Calculation

```bash
cd 01_scf
# Submit job (adjust for your cluster)
sbatch submit.slurm
# Or run directly
mpirun -np 64 vasp_std > vasp.out

# Check convergence
grep "E0=" OSZICAR | tail -1
grep "mag=" OSZICAR | tail -1
```

**Expected Runtime:** 2-6 hours (depends on hardware)

**Convergence Criteria:**
- Electronic: EDIFF = 1e-7 eV
- Energy change < 1 meV between steps
- Magnetic moments stabilized

**Key Output Files:**
- `OUTCAR` - Detailed output (magnetic moments, forces)
- `CHGCAR` - Charge density (for DOS/band)
- `WAVECAR` - Wave functions
- `vasprun.xml` - Complete run data

### Step 3: Run DOS Calculation

```bash
# Copy CHGCAR from SCF
cp 01_scf/CHGCAR 02_dos/

cd 02_dos
mpirun -np 64 vasp_std > vasp.out

# Quick check
grep "E-fermi" OUTCAR
```

**Expected Runtime:** 1-3 hours

**Key Output Files:**
- `DOSCAR` - Total and projected DOS
- `PROCAR` - Orbital-projected data

**Post-processing:**
```python
# Using vaspvis or pymatgen
from pymatgen.io.vasp import Vasprun
vr = Vasprun("vasprun.xml")
dos = vr.complete_dos
dos.get_gap()  # Get band gap
```

### Step 4: Run Band Structure

```bash
# Copy CHGCAR from SCF
cp 01_scf/CHGCAR 03_band_unfold/

cd 03_band_unfold
mpirun -np 64 vasp_std > vasp.out
```

**Expected Runtime:** 1-2 hours

**Key Output Files:**
- `WAVECAR` - For band unfolding (large file!)
- `EIGENVAL` - Eigenvalues along k-path
- `PROCAR` - Orbital projections

**Band Unfolding:**
```bash
# Using easyunfold (recommended)
pip install easyunfold

# Generate unfolding
easyunfold unfold calculate \
    --matrix "[[2,1,0],[1,2,0],[0,0,2]]" \
    --out-file unfold.json

# Plot
easyunfold unfold plot unfold.json --output band.png
```

## Detailed Parameter Explanations

### Critical Parameters

#### 1. LSORBIT = True (Spin-Orbit Coupling)
- **Purpose:** Enables non-collinear magnetism with SOC
- **Effect:** Spin is 3D vector (mx, my, mz)
- **Cost:** ~2× slower, 2× memory
- **When to use:** Essential for:
  - Non-collinear magnetic structures
  - Heavy elements (Pd has weak SOC)
  - Spin Hall effects, topological properties

#### 2. MAGMOM (Initial Magnetic Moments)
```
MAGMOM = 36*0 \         # O atoms: non-magnetic
    <18 Cr moments> \   # Cr: 3 components (mx, my, mz) × 18 atoms
    18*0                # Pd: non-magnetic
```
- **Format:** Each Cr has 3 values (mx, my, mz)
- **Source:** From model/MAGMOM.txt
- **Important:** Order must match POSCAR atom sequence
- **Tip:** Check final moments in OUTCAR: `grep "magnetization (x,y,z)" OUTCAR`

#### 3. DFT+U Parameters
```
LDAU = True            # Enable DFT+U
LDAUTYPE = 2           # Dudarev: U_eff = U - J
LDAUL = -1 2 -1        # Apply U to Cr 3d only
LDAUU = 0 4 0          # U = 4.0 eV for Cr
LDAUJ = 0 0.9 0        # J = 0 (not used in Dudarev)
```
- **Why U=4.0 eV:** Standard for Cr 3d states
- **Convergence:** Test U = 3.0, 4.0, 5.0 eV
- **Effect:** Opens band gap, localizes d-electrons
- **Check:** `grep "LDAU:" OUTCAR` for occupation matrix

#### 4. K-point Convergence
| Calculation | Mesh     | Reason                        |
|-------------|----------|-------------------------------|
| SCF         | 8×8×1    | Balance speed/accuracy        |
| DOS         | 12×12×1  | Tetrahedron method needs dense|
| Band        | Path     | High-symmetry lines           |

- **c-axis:** Only 1 k-point (supercell is 36 Å, negligible dispersion)
- **ab-plane:** Dense sampling for 2D electronic structure

#### 5. ICHARG Settings
| Value | Meaning                    | Use Case          |
|-------|----------------------------|-------------------|
| 2     | Superposition of atoms     | Fresh SCF         |
| 11    | Read CHGCAR, no update     | DOS, Band (faster)|
| 1     | Read CHGCAR, continue SCF  | Relaxation restart|

### Parallelization Guide

**NPAR, KPAR, NCORE relationship:**
```
Total cores ≈ NPAR × KPAR × (bands / NPAR)
```

**Recommended settings:**
| Cores | NPAR | KPAR | NCORE | Notes                  |
|-------|------|------|-------|------------------------|
| 32    | 16   | 2    | 8     | Small cluster          |
| 64    | 32   | 2    | 8     | Medium (recommended)   |
| 128   | 32   | 4    | 8     | Large cluster          |
| 256   | 64   | 4    | 8     | HPC                    |

**Tips:**
- KPAR = number of k-points that can be computed in parallel
- NCORE = cores per orbital (4-8 optimal for modern CPUs)
- NPAR divides bands among cores
- For SCF: Tune NPAR first, then KPAR
- For DOS/Band: Higher KPAR helps (more k-points)

## Common Issues and Solutions

### 1. SCF Not Converging
**Symptoms:** EDIFF not reached after NELM steps

**Solutions:**
```
# Increase mixing
BMIX = 5
AMIX = 0.2
AMIN = 0.01

# Or use better algorithm
ALGO = All

# Or start from smaller U
LDAUU = 0 3 0  # Instead of 4
```

### 2. Wrong Magnetic State
**Symptoms:** Magnetic moments flip or become zero

**Solutions:**
- Check MAGMOM initialization matches POSCAR order
- Increase NELM to 600
- Try different starting MAGMOM (add ±10% random noise)
- Check if structure is correct: `grep "magnetization (x,y,z)" OUTCAR | head -20`

### 3. Memory Issues
**Symptoms:** Killed by OOM (out of memory)

**Solutions:**
```
# Reduce memory
LREAL = Auto      # Use real-space projection
NPAR = <larger>   # More cores per band (less memory each)

# Or split calculation
# Run with smaller ENCUT first, then restart
```

### 4. Negative Frequency in OSZICAR
**Symptoms:** `WARNING: Sub-Space-Matrix is not hermitian`

**Solutions:**
- Usually harmless for non-collinear + SOC
- If persistent: Reduce BMIX
- Check POSCAR for symmetry issues

## Validation Checklist

After calculations, verify:

### SCF (01_scf/)
- [ ] `grep "reached required accuracy" OUTCAR` shows convergence
- [ ] Magnetic moments match expected values (±10%)
- [ ] Energy converged: `tail OSZICAR` shows E change < 1 meV
- [ ] CHGCAR file exists and > 100 MB

### DOS (02_dos/)
- [ ] Fermi level reasonable: `grep "E-fermi" OUTCAR`
- [ ] Band gap matches literature (if known)
- [ ] DOS smooth without spikes
- [ ] Cr 3d states prominent near Fermi level

### Band (03_band_unfold/)
- [ ] WAVECAR exists (large file, >1 GB)
- [ ] Band structure shows expected features
- [ ] Unfolding recovers primitive cell bands
- [ ] Gap/dispersion matches DOS

## Next Steps

1. **Analysis:**
   - Plot DOS with vaspvis/pymatgen
   - Unfold band structure
   - Extract magnetic anisotropy energy
   - Calculate Berry curvature (if needed)

2. **Doping Studies:**
   - Generate SQS structures
   - Modify POSCAR for dopants
   - Re-run calculations

3. **Advanced Calculations:**
   - Berry phase (LBERRY = True)
   - Optical properties (LOPTICS = True)
   - Phonons (IBRION = 5,6,7,8)

## References

- VASP Wiki: https://www.vasp.at/wiki/
- DFT+U: Dudarev et al., PRB 57, 1505 (1998)
- PdCrO2 magnetism: Takatsu et al., PRB 89, 104408 (2014)
- Band unfolding: easyunfold documentation
