# Magnetic Systems: Non-collinear DFT Workflow

Complete workflow for VASP non-collinear magnetic calculations, using PdCrO₂ delafossite as a reference system.

## Overview

This directory contains a production-ready workflow for calculating magnetic systems with:
- **Non-collinear magnetism** with spin-orbit coupling (SOC)
- **Complex magnetic structures** (18-sublattice ordering)
- **DFT+U corrections** for strongly correlated electrons
- **Band structure unfolding** for supercell calculations

**Reference system:** PdCrO₂ with 6-layer magnetic structure (Takatsu et al., PRB 89, 104408, 2014)

## Directory Structure

```
magnetic-systems/
├── README.md                      # This file - workflow overview
├── noncollinear-magnetism.md     # Theory and practical guide
├── band-unfolding.md             # Supercell band unfolding guide
├── model/                        # Structure preparation
│   ├── shift_z.py                 # Reorient cell (Cr to bottom)
│   ├── make_supercell.py          # Generate magnetic supercell
│   ├── generate_magmom.py         # Calculate initial magnetic moments
│   ├── generate_mcif.py           # VESTA visualization (mcif format)
│   ├── 01_POSCAR_unitcell         # Primitive cell
│   ├── 02_POSCAR_Cr-bottom        # Cr-terminated cell
│   ├── 03_POSCAR_supercell_6x     # 2×2×2 magnetic supercell
│   ├── MAGMOM.txt                 # Magnetic moment vectors
│   └── magnetic_structure.mcif    # VESTA visualization file
└── calculation/                  # DFT calculations
    ├── 01_scf/                    # Self-consistent field
    │   ├── INCAR                   # 103 lines with detailed comments
    │   └── KPOINTS                 # 8×8×1 mesh
    ├── 02_dos/                    # Density of states
    │   ├── INCAR                   # DOS-specific settings
    │   └── KPOINTS                 # 12×12×1 denser mesh
    ├── 03_band_unfold/            # Band structure + unfolding
    │   ├── INCAR                   # LWAVE=True for unfolding
    │   └── KPOINTS                 # Γ-M-K-Γ high-symmetry path
    ├── README.md                   # Detailed calculation guide
    ├── QUICK_START.md              # 5-minute setup guide
    ├── FILE_SUMMARY.md             # File descriptions
    ├── POTCAR_README.txt           # POTCAR generation
    ├── check_setup.sh              # Pre-flight validation script
    └── submit_template.slurm       # Job submission template
```

## Workflow Overview

### Phase 1: Structure Preparation (`model/`)

**Goal:** Transform primitive cell → magnetic supercell with proper termination

```bash
cd model/

# Step 1: Reorient cell (Cr layer to bottom)
python shift_z.py 01_POSCAR_unitcell 02_POSCAR_Cr-bottom

# Step 2: Generate supercell with transformation matrix
python make_supercell.py 02_POSCAR_Cr-bottom 03_POSCAR_supercell_6x

# Step 3: Generate magnetic moments from Takatsu model
python generate_magmom.py

# Step 4: Generate VESTA visualization
python generate_mcif.py 03_POSCAR_supercell_6x
```

**Key output:** `03_POSCAR_supercell_6x` (72 atoms: 36 O, 18 Cr, 18 Pd)

**Transformation matrix:**
```python
[[2, 1, 0],
 [1, 2, 0],
 [0, 0, 2]]
```

### Phase 2: DFT Calculations (`calculation/`)

**Goal:** SCF → DOS → Band structure with full documentation

#### 2.1 Self-Consistent Field (SCF)

```bash
cd calculation/01_scf/
cp ../../model/03_POSCAR_supercell_6x POSCAR

# Generate POTCAR (O, Cr_pv, Pd order)
cat $VASP_PP_PATH/potpaw_PBE/O/POTCAR \
    $VASP_PP_PATH/potpaw_PBE/Cr_pv/POTCAR \
    $VASP_PP_PATH/potpaw_PBE/Pd/POTCAR > POTCAR

# Validate setup
../check_setup.sh

# Submit job
sbatch ../submit_template.slurm
```

**Key INCAR parameters:**
```bash
LSORBIT = True        # Enable SOC (auto-enables non-collinear)
MAGMOM = 36*0 \       # O: non-magnetic
    0.250611  0.150583  0.956305 \  # Cr1: 3D moment vector
    [... 18 Cr magnetic moments ...]
    18*0              # Pd: non-magnetic

LDAU = True           # DFT+U correction
LDAUU = 0 4 0         # U = 4.0 eV for Cr only
```

#### 2.2 Density of States (DOS)

```bash
cd calculation/02_dos/
cp ../01_scf/{POSCAR,POTCAR,CHGCAR} .

# INCAR differences from SCF:
# ICHARG = 11    (read CHGCAR, non-self-consistent)
# ISMEAR = -5    (tetrahedron method)
# LORBIT = 11    (orbital-projected DOS)
# NEDOS = 3001   (high resolution)
```

#### 2.3 Band Structure + Unfolding

```bash
cd calculation/03_band_unfold/
cp ../01_scf/{POSCAR,POTCAR,CHGCAR} .

# CRITICAL: LWAVE = True (required for band unfolding)
sbatch ../submit_template.slurm

# After convergence, unfold bands
easyunfold unfold calculate --matrix "[[2,1,0],[1,2,0],[0,0,2]]"
```

## Key Features

### 1. Comprehensive Documentation

Every file is thoroughly documented:
- **INCAR files:** 100+ lines with line-by-line parameter explanations
- **KPOINTS files:** Rationale for mesh density, convergence notes
- **Scripts:** Detailed docstrings and usage examples
- **Guides:** README.md (detailed), QUICK_START.md (5 min), FILE_SUMMARY.md (reference)

### 2. Non-collinear Magnetism Setup

**18-sublattice magnetic structure:**
- 6 Cr layers (z-direction)
- 3 sublattices per layer (A, B, C in triangular lattice)
- Each Cr atom has unique 3D magnetic moment vector

**MAGMOM format:**
```bash
MAGMOM = mx1 my1 mz1  mx2 my2 mz2  ...  (3 components per atom)
```

**Computed from Takatsu Model 4:**
- Layer-dependent tilt angles α_n
- In-plane rotation angles φ_n
- Magnetic moment magnitudes ξ_n

### 3. DFT+U for Cr 3d Electrons

```bash
LDAU = True
LDAUTYPE = 2          # Dudarev formulation
LDAUL = -1 2 -1       # O=off, Cr=d, Pd=off
LDAUU = 0 4 0         # U = 4.0 eV for Cr
LMAXMIX = 4           # Mix d orbitals
```

**U value source:** Literature for Cr oxides (3-4 eV range)

### 4. Band Unfolding for Supercells

**Problem:** Supercell band structure is folded and messy

**Solution:** Unfold to primitive cell Brillouin zone

**Tools:** easyunfold or BandUP

**Requirements:**
- LWAVE = True in INCAR (writes WAVECAR)
- Transformation matrix from supercell construction
- High-symmetry k-path in primitive BZ (Γ-M-K-Γ)

### 5. Validation & Quality Control

**Pre-run checks (`check_setup.sh`):**
- File existence (INCAR, POSCAR, POTCAR, KPOINTS)
- POSCAR atom count consistency
- POTCAR element order matches POSCAR
- MAGMOM count validation
- INCAR parameter sanity checks

**Post-run analysis (`submit_template.slurm`):**
- Convergence detection
- Final energy extraction
- Magnetic moment summary
- Warning/error scanning

## Usage Scenarios

### For New Users

1. Read `calculation/QUICK_START.md` (5 minutes)
2. Run validation: `cd calculation/01_scf && ../check_setup.sh`
3. Submit first job: `sbatch ../submit_template.slurm`
4. Refer to `calculation/README.md` for troubleshooting

### For Experienced Users

1. Copy POSCAR, generate POTCAR
2. Quick validate with `check_setup.sh`
3. Customize parallelization (NPAR, KPAR, NCORE)
4. Submit with optimized resources

### For Similar Systems

**To adapt this workflow for other magnetic materials:**

1. **Structure preparation:**
   - Replace primitive POSCAR
   - Adjust supercell transformation matrix in `make_supercell.py`

2. **Magnetic structure:**
   - Modify `generate_magmom.py` for your magnetic model
   - Update layer count, sublattice count, spin directions

3. **DFT+U values:**
   - Update LDAUU in INCAR based on literature
   - Consider different U values for different magnetic elements

4. **K-point mesh:**
   - Adjust based on supercell size
   - Denser for smaller cells, sparser for larger

## Technical Details

### Computational Requirements

**For 72-atom PdCrO₂ supercell:**

| Calculation | Cores | Time    | Memory | Key Files         |
|-------------|-------|---------|--------|-------------------|
| SCF         | 64    | 2-6 hrs | ~20 GB | CHGCAR (for DOS)  |
| DOS         | 64    | 1-3 hrs | ~15 GB | DOSCAR, PROCAR    |
| Band        | 64    | 1-2 hrs | ~15 GB | WAVECAR (>1 GB)   |

**Parallelization (INCAR):**
```bash
NPAR = 8              # Parallel bands
KPAR = 2              # Parallel k-points
NCORE = 4             # Cores per band
```

### Convergence Parameters

```bash
ENCUT = 400           # Plane-wave cutoff (eV)
EDIFF = 1e-7          # Electronic convergence (eV)
ALGO = Normal         # Davidson + RMM-DIIS
NELM = 200            # Max electronic steps
```

### K-point Convergence

| Calculation | Mesh    | Total k-pts | Purpose              |
|-------------|---------|-------------|----------------------|
| SCF         | 8×8×1   | 64          | Total energy         |
| DOS         | 12×12×1 | 144         | Smooth DOS curves    |
| Band        | 40/seg  | ~120        | High-symmetry path   |

**Rationale:**
- Dense in ab-plane (5.07 Å lattice)
- Sparse along c-axis (36.24 Å supercell)

## Related Topics

- **[noncollinear-magnetism.md](noncollinear-magnetism.md)** - Theory, MAGMOM format, troubleshooting
- **[band-unfolding.md](band-unfolding.md)** - Detailed unfolding guide with easyunfold examples

## References

1. **Takatsu et al., Phys. Rev. B 89, 104408 (2014)**
   - Experimental magnetic structure determination
   - Model 4 parameters used in this workflow

2. **VASP Wiki:**
   - [Non-collinear calculations](https://www.vasp.at/wiki/index.php/Non-collinear_calculations)
   - [LSORBIT](https://www.vasp.at/wiki/index.php/LSORBIT)
   - [DFT+U](https://www.vasp.at/wiki/index.php/LDAU)

3. **easyunfold:**
   - [GitHub Repository](https://github.com/SMTG-UCL/easyunfold)
   - [Documentation](https://smtg-ucl.github.io/easyunfold/)

## Troubleshooting

### Common Issues

**1. Magnetic moments collapse to zero**
- Check MAGMOM order matches POSCAR atom sequence
- Increase initial magnetic moment magnitudes
- Try ISYM = 0 to break symmetry

**2. SCF not converging**
- Increase NELM (200 → 500)
- Try ALGO = All or ALGO = VeryFast
- Adjust mixing: BMIX, AMIX

**3. POTCAR errors**
- Verify element order: `grep TITEL POTCAR`
- Must match POSCAR line 6: O Cr Pd
- Use Cr_pv (not Cr) for better semicore description

**4. Band unfolding fails**
- Ensure LWAVE = True in INCAR
- Check transformation matrix matches supercell construction
- Verify WAVECAR file exists and is not corrupted

### Quick Diagnostics

```bash
# Check convergence
grep "reached required accuracy" OUTCAR

# View final energy
grep "E0=" OSZICAR | tail -1

# Check magnetic moments
grep "magnetization (x,y,z)" OUTCAR | head -50

# Verify file sizes
ls -lh CHGCAR WAVECAR
```

## Contributing

To improve this workflow:
1. Test on different magnetic systems
2. Document parameter optimization strategies
3. Add analysis scripts for magnetic structure
4. Expand troubleshooting guide with new issues

## License

MIT License - Free to use and modify

---

**Last updated:** 2025-01-30
**VASP version:** 6.x
**Tested on:** SLURM-based HPC clusters
