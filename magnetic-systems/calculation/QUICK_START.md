# Quick Start Guide

Get your PdCrO₂ DFT calculations running in 5 minutes!

## Prerequisites

✓ VASP installed and accessible
✓ VASP pseudopotential library (potpaw_PBE)
✓ Access to computing cluster (or local machine with MPI)
✓ Completed model preparation (../model/)

## Step-by-Step Setup

### 1. Copy Structure Files (1 min)

```bash
cd calculation/01_scf

# Copy POSCAR from model directory
cp ../../model/03_POSCAR_supercell_6x POSCAR

# Verify
head -7 POSCAR
# Should show: O Cr Pd with 36 18 18 atoms
```

### 2. Generate POTCAR (1 min)

```bash
# Method A: Direct concatenation
cat $VASP_PP_PATH/potpaw_PBE/O/POTCAR \
    $VASP_PP_PATH/potpaw_PBE/Cr_pv/POTCAR \
    $VASP_PP_PATH/potpaw_PBE/Pd/POTCAR > POTCAR

# Method B: Using vaspkit
vaspkit
# → 103 (Generate POTCAR)
# → 1 (PBE)

# Verify POTCAR
grep TITEL POTCAR
# Should show: O, Cr_pv, Pd
```

**Important:** Element order must match POSCAR!

### 3. Verify Setup (30 sec)

```bash
# Run checker script
../check_setup.sh

# Should show all green checkmarks ✓
# Fix any warnings before proceeding
```

### 4. Adjust Submit Script (1 min)

```bash
# Copy template
cp ../submit_template.slurm submit.slurm

# Edit for your cluster
nano submit.slurm

# Adjust:
# - Partition name
# - Number of nodes/cores
# - Module loads
# - Email address
```

### 5. Submit Job (30 sec)

```bash
# Submit
sbatch submit.slurm

# Monitor
squeue -u $USER
tail -f vasp_*.out

# Check convergence
grep "E0=" OSZICAR | tail -5
```

## After SCF Converges

### Run DOS Calculation

```bash
cd ../02_dos

# Copy files
cp ../01_scf/POSCAR .
cp ../01_scf/POTCAR .
cp ../01_scf/CHGCAR .  # IMPORTANT!
cp ../01_scf/submit.slurm .

# Edit submit.slurm: change job-name to pdcro2_dos
nano submit.slurm

# Submit
sbatch submit.slurm
```

### Run Band Structure

```bash
cd ../03_band_unfold

# Copy files
cp ../01_scf/POSCAR .
cp ../01_scf/POTCAR .
cp ../01_scf/CHGCAR .  # IMPORTANT!
cp ../01_scf/submit.slurm .

# Edit submit.slurm: change job-name to pdcro2_band
nano submit.slurm

# Submit
sbatch submit.slurm
```

## Quick Checks

### Is SCF converged?

```bash
grep "reached required accuracy" OUTCAR
# Should return a line if converged

# Or check energy
grep "E0=" OSZICAR | tail -1
# Energy should stabilize
```

### Check magnetic moments

```bash
grep "magnetization (x,y,z)" OUTCAR | head -50
# Shows mx, my, mz for first atoms
# Cr atoms should have moments ~3 μB
```

### Estimate runtime

For PdCrO₂ (72 atoms, 8×8×1 k-points):
- **SCF:** 2-6 hours (64 cores)
- **DOS:** 1-3 hours (64 cores)
- **Band:** 1-2 hours (64 cores)

**Total:** ~1 day for complete workflow

## Troubleshooting

### Job dies immediately

```bash
# Check error file
cat vasp_*.err

# Common issues:
# 1. Wrong modules loaded
# 2. POTCAR missing/corrupted
# 3. Not enough memory
```

### SCF not converging

```bash
# Check OSZICAR
tail -20 OSZICAR

# If energy oscillating:
# → Increase BMIX in INCAR
# → Try ALGO = All

# If magnetic moments unstable:
# → Check MAGMOM order matches POSCAR
```

### Out of memory

```bash
# Reduce memory usage in INCAR:
# NPAR = <larger number>
# LREAL = Auto

# Or request more nodes
```

## Expected Results

### SCF (01_scf/)
- Energy: ~-500 to -600 eV (total)
- Cr magnetic moments: ~2-3 μB
- Runtime: 2-6 hours
- Key files: OUTCAR, CHGCAR (for DOS/band)

### DOS (02_dos/)
- Fermi level: ~5-7 eV (check OUTCAR)
- Band gap: ~0-2 eV (depends on U value)
- Strong Cr 3d features near E_F
- Key files: DOSCAR, PROCAR

### Band (03_band_unfold/)
- Typical metal/semiconductor bands
- Need unfolding for clean primitive cell bands
- Key files: WAVECAR (for unfolding), EIGENVAL

## Next Steps

1. **Analyze results:**
   ```bash
   # Extract data
   grep "E-fermi" */OUTCAR

   # Plot DOS
   python plot_dos.py  # (create this script)
   ```

2. **Unfold bands:**
   ```bash
   cd 03_band_unfold
   easyunfold unfold calculate --matrix "[[2,1,0],[1,2,0],[0,0,2]]"
   ```

3. **Validate:**
   - Compare with literature (if available)
   - Check magnetic moments make sense
   - Verify band gap/metallicity

4. **Document:**
   - Save key results (energies, gaps, moments)
   - Archive OUTCAR, DOSCAR, unfolded bands
   - Write notes on what worked/didn't work

## Get Help

- VASP Wiki: https://www.vasp.at/wiki/
- Email your cluster support for resource questions
- Check ../README.md for detailed parameter explanations
