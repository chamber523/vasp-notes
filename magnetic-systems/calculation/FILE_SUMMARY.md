# Calculation Directory File Summary

Complete overview of all files in the calculation directory.

## Directory Structure

```
calculation/
├── 01_scf/                      # Self-consistent field calculation
│   ├── INCAR                    # 103 lines, detailed parameter explanations
│   └── KPOINTS                  # 8×8×1 mesh with convergence notes
│
├── 02_dos/                      # Density of states calculation
│   ├── INCAR                    # DOS-specific settings (ICHARG=11, ISMEAR=-5)
│   └── KPOINTS                  # 12×12×1 denser mesh
│
├── 03_band_unfold/              # Band structure calculation
│   ├── INCAR                    # Band settings (LWAVE=True for unfolding)
│   └── KPOINTS                  # Γ-M-K-Γ high-symmetry path
│
├── README.md                    # Main documentation (300+ lines)
├── QUICK_START.md               # 5-minute setup guide
├── POTCAR_README.txt            # POTCAR generation instructions
├── check_setup.sh               # Validation script (executable)
├── submit_template.slurm        # SLURM job template (executable)
└── FILE_SUMMARY.md              # This file
```

## File Details

### Main Documentation

#### README.md (9.8 KB)
**Purpose:** Comprehensive guide for all calculations

**Contents:**
- Directory structure overview
- Complete workflow (SCF → DOS → Band)
- Detailed parameter explanations
  - LSORBIT (SOC) rationale
  - MAGMOM format and source
  - DFT+U setup (U=4.0 eV for Cr)
  - K-point convergence guide
  - ICHARG modes explained
- Parallelization guide (NPAR, KPAR, NCORE)
- Common issues & solutions
- Validation checklist
- Expected results & runtime
- References

**Key Sections:**
- Critical Parameters (SOC, MAGMOM, DFT+U, k-points, ICHARG)
- Parallelization table (32-256 cores)
- Troubleshooting (4 common issues)
- Validation checklist (3 stages)

#### QUICK_START.md (4.2 KB)
**Purpose:** Get calculations running in 5 minutes

**Contents:**
- Prerequisites checklist
- 5-step setup (1 min each)
- Quick verification commands
- Expected results summary
- Troubleshooting shortcuts

**Use When:** First time setup or need quick reference

#### POTCAR_README.txt (1.8 KB)
**Purpose:** POTCAR generation guide

**Contents:**
- 3 methods (VASP PP library, vaspkit, pymatgen)
- Element order verification
- Pseudopotential recommendations (Cr_pv vs Cr vs Cr_sv)
- Validation commands
- Important warnings (licensing, git)

### Input Files

#### 01_scf/INCAR (4.1 KB, 103 lines)
**Calculation Type:** Self-consistent field

**Key Parameters:**
```
ALGO = Normal          # Davidson + RMM-DIIS
ENCUT = 400           # Plane-wave cutoff
EDIFF = 1e-7          # Tight convergence
LSORBIT = True        # Non-collinear + SOC
MAGMOM = 36*0 + 18 Cr moments + 18*0
LDAU = True           # DFT+U
LDAUU = 0 4 0         # U = 4.0 eV for Cr
NPAR/KPAR/NCORE       # Parallelization
```

**Special Features:**
- Line-by-line comments (every parameter explained)
- MAGMOM from ../model/MAGMOM.txt
- Optimized for 64-core jobs

#### 01_scf/KPOINTS (0.8 KB)
**Mesh:** 8×8×1 (Monkhorst-Pack)

**Rationale:**
- Dense in ab-plane (5.07 Å cell)
- Sparse along c (36.24 Å supercell)
- Converged for total energy

**Includes:**
- Convergence test template
- Recommendations for other calculations

#### 02_dos/INCAR (4.0 KB)
**Calculation Type:** Density of states

**Differences from SCF:**
```
ICHARG = 11           # Read CHGCAR (non-SCF)
ISMEAR = -5           # Tetrahedron method
LORBIT = 11           # Orbital projection
NEDOS = 3001          # High DOS resolution
LCHARG = False        # Don't write CHG*
LWAVE = False         # Don't write WAVECAR
```

**Prerequisites:**
- Converged CHGCAR from 01_scf/

**Output:**
- DOSCAR (total + projected DOS)
- PROCAR (orbital-decomposed)

#### 02_dos/KPOINTS (0.7 KB)
**Mesh:** 12×12×1 (denser than SCF)

**Why Denser:**
- Required for ISMEAR=-5
- Smoother DOS curves
- Better gap estimation

#### 03_band_unfold/INCAR (4.3 KB)
**Calculation Type:** Band structure

**Differences from DOS:**
```
ISMEAR = 0            # Gaussian (not tetrahedron)
LWAVE = True          # REQUIRED for unfolding
```

**Critical Note:**
- LWAVE=True creates large WAVECAR (>1 GB)
- Needed for easyunfold/BandUP

#### 03_band_unfold/KPOINTS (1.3 KB)
**Path:** Γ → M → K → Γ (hexagonal symmetry)

**Format:** Line mode (40 points/segment)

**Includes:**
- Unfolding instructions (easyunfold)
- Transformation matrix reminder
- K-point density recommendations

### Utility Scripts

#### check_setup.sh (4.5 KB, executable)
**Purpose:** Pre-flight validation

**Checks:**
1. Required files exist (INCAR, POSCAR, POTCAR, KPOINTS)
2. POSCAR atom count vs coordinates
3. POTCAR count vs elements
4. INCAR settings (SOC, DFT+U, MAGMOM)
5. KPOINTS mesh reasonableness
6. Resource estimation

**Usage:**
```bash
cd 01_scf
../check_setup.sh
# Shows ✓/✗/⚠ for each check
```

**Output:** Color-coded (green/yellow/red)

#### submit_template.slurm (4.8 KB, executable)
**Purpose:** Universal job submission template

**Features:**
- Pre-run file checks
- Key parameter display
- 3 MPI launch methods (srun/mpirun/mpiexec)
- Post-run analysis (convergence, energy, moments)
- Warning detection
- Optional cleanup

**Customization Points:**
- Partition name (#SBATCH --partition)
- Node/core counts (#SBATCH --nodes/tasks)
- Module loads (module load ...)
- Email (#SBATCH --mail-user)

**Includes:**
- Comments for Intel vs GNU environments
- Hybrid MPI+OpenMP support
- Runtime estimates

## Key Features Across All Files

### 1. Comprehensive Comments
Every INCAR parameter has inline explanation:
```
ENCUT = 400     # Plane-wave cutoff: 400 eV
                # Sufficient for O, Cr, Pd (check POTCAR ENMAX)
```

### 2. Cross-Referencing
Files reference each other:
- INCAR → mentions MAGMOM from ../model/
- KPOINTS → references convergence tests
- README → points to specific line numbers

### 3. Universal Applicability
Templates work for:
- Different cluster systems (SLURM examples)
- Various VASP versions (6.x compatible)
- Different parallelization schemes

### 4. Error Prevention
Multiple validation layers:
- check_setup.sh (before job)
- submit_template.slurm (during job)
- README troubleshooting (after job)

## Usage Workflow

### For First-Time Users:
1. Read QUICK_START.md (5 min)
2. Run check_setup.sh (30 sec)
3. Submit job (30 sec)
4. Refer to README for issues

### For Experienced Users:
1. Copy POSCAR, generate POTCAR
2. Quick verify with check_setup.sh
3. Submit with custom parameters
4. Use README as reference

### For Debugging:
1. Check README troubleshooting section
2. Run check_setup.sh for validation
3. Grep OUTCAR for specific warnings
4. Compare with expected results

## File Size Summary

| File | Size | Lines | Purpose |
|------|------|-------|---------|
| README.md | 9.8 KB | 300+ | Main documentation |
| QUICK_START.md | 4.2 KB | 200+ | Quick guide |
| 01_scf/INCAR | 4.1 KB | 103 | SCF parameters |
| 02_dos/INCAR | 4.0 KB | 100 | DOS parameters |
| 03_band_unfold/INCAR | 4.3 KB | 105 | Band parameters |
| submit_template.slurm | 4.8 KB | 150+ | Job script |
| check_setup.sh | 4.5 KB | 180+ | Validation |
| POTCAR_README.txt | 1.8 KB | 60+ | POTCAR guide |
| *KPOINTS files | ~1 KB each | 20-30 | k-point meshes |

**Total documentation:** ~35 KB, 1000+ lines

## Updates and Maintenance

### When to Update:

1. **VASP version change:**
   - Check new parameters in INCAR
   - Update module loads in submit script

2. **New literature:**
   - Update U values if better references
   - Adjust MAGMOM if new magnetic structure

3. **Cluster changes:**
   - Update partition names
   - Adjust parallelization recommendations

4. **Bug fixes:**
   - Document in README troubleshooting
   - Update check_setup.sh validation

### Version Control:
- All files are plain text (git-friendly)
- POTCAR excluded (licensing)
- Track major changes in README

## Additional Notes

### What's NOT Included:
- POTCAR files (license restrictions)
- Output files (OUTCAR, CHGCAR, etc.)
- Analysis scripts (user-specific)
- Visualization tools (separate package)

### What Users Need to Provide:
- VASP license & installation
- Pseudopotential library
- Computing resources
- Specific POSCAR (from ../model/)

### Future Extensions:
- Relaxation calculations (IBRION=2)
- Berry phase calculations (LBERRY=True)
- Phonon calculations (IBRION=5,6,7,8)
- Advanced analysis scripts
