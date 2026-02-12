# Choosing the Right Parallelization Script

Quick guide to select the appropriate SLURM script for your VASP calculation.

## Decision Tree

```
Start
  │
  ├─ Using Hybrid Functional (HSE/PBE0)?
  │   ├─ YES → GPU Available?
  │   │   ├─ YES → gpu_hybrid_functional.slurm
  │   │   └─ NO  → cpu_large_system.slurm (expect long runtime)
  │   │
  │   └─ NO → How many atoms?
  │       ├─ < 50 atoms      → cpu_small_system.slurm
  │       ├─ 50-200 atoms    → cpu_medium_system.slurm
  │       └─ > 200 atoms     → cpu_large_system.slurm
```

## Script Selection Guide

### cpu_small_system.slurm

**Use for:**
- Small unit cells (< 50 atoms)
- Quick test calculations
- Band structure calculations (few k-points)
- Relaxation of molecules

**Resources:**
- 1 node, 32 MPI ranks
- 30 minutes (debug queue)

**Example systems:**
- Diamond unit cell (2 atoms)
- GaAs supercell (32 atoms)
- Small molecules

---

### cpu_medium_system.slurm

**Use for:**
- Medium-sized systems (50-200 atoms)
- Non-collinear magnetic calculations
- Systems with moderate k-points
- Standard DFT+U calculations

**Resources:**
- 2 nodes, 32 MPI ranks
- 12 hours (regular queue)

**Example systems:**
- PdCrO₂ supercell (72 atoms)
- Interface calculations (~100 atoms)
- Doped semiconductors

---

### cpu_large_system.slurm

**Use for:**
- Large systems (> 200 atoms)
- Dense k-point meshes (DOS calculations)
- Hybrid functionals on CPU
- Long MD simulations

**Resources:**
- 8 nodes, 128 MPI ranks
- 2 days (regular queue)

**Example systems:**
- Large supercells (500+ atoms)
- Interface systems with vacuum
- Disordered/amorphous systems
- HSE calculations on CPU

---

### gpu_standard.slurm

**Use for:**
- Standard DFT on GPU
- Medium systems with GPU available
- Testing GPU performance
- GGA/PBE calculations

**Resources:**
- 1 GPU node, 4 A100 GPUs
- 12 hours (regular queue)

**Example systems:**
- Co₂XY Heusler compounds
- Standard band/DOS calculations
- Any GGA calculation

**Performance:**
- ~5-10× faster than CPU for GGA
- Best for 50-500 atoms

---

### gpu_hybrid_functional.slurm

**Use for:**
- HSE/PBE0 hybrid functionals
- GW calculations
- Accurate band gaps
- Spin-orbit coupling + HSE

**Resources:**
- 2 GPU nodes, 8 A100 GPUs
- 24 hours (regular queue)

**Example systems:**
- HSE band structure
- GW quasiparticle energies
- Accurate optical properties

**Performance:**
- ~20-50× faster than CPU for HSE
- Essential for HSE on large systems

## System Size vs. Resources

| Atoms | k-points | RAM needed | Recommended Script |
|-------|----------|------------|-------------------|
| < 50  | < 20     | < 50 GB    | cpu_small_system  |
| 50-100 | 10-50   | 50-100 GB  | cpu_medium_system |
| 100-200 | 10-50  | 100-200 GB | cpu_medium_system |
| 200-500 | < 20   | 200-500 GB | cpu_large_system  |
| > 500 | < 10    | > 500 GB   | cpu_large_system (multiple nodes) |

**Add GPU for:**
- HSE/PBE0 calculations (any size)
- Systems where time > 4 hours on CPU

## Calculation Type Recommendations

### SCF (Self-Consistent Field)

```bash
# Small system: 1 node
cpu_small_system.slurm

# Medium system: 2 nodes
cpu_medium_system.slurm

# Large system: 4-8 nodes
cpu_large_system.slurm
```

**Adjust:** `KPAR` based on number of k-points

---

### DOS (Density of States)

```bash
# Many k-points needed → use larger nodes
cpu_medium_system.slurm  # For quick DOS
cpu_large_system.slurm   # For high-quality DOS
```

**INCAR settings:**
```bash
ISMEAR = -5     # Tetrahedron method
NEDOS = 3001    # Dense DOS grid
KPAR = 8        # High for many k-points
```

---

### Band Structure

```bash
# Few k-points on path → use small
cpu_small_system.slurm
```

**INCAR settings:**
```bash
KPAR = 1        # Only 1 for line-mode KPOINTS
NCORE = 8       # Higher NCORE for fewer k-points
```

---

### HSE (Hybrid Functional)

```bash
# Always prefer GPU if available
gpu_hybrid_functional.slurm

# If no GPU, use maximum CPU resources
cpu_large_system.slurm (8+ nodes)
```

**INCAR settings:**
```bash
LHFCALC = True
PRECFOCK = Fast  # For GPU
KPAR = 4-8       # Balance k-point and band parallelization
```

---

### Magnetic Calculations (Non-collinear)

```bash
# Non-collinear requires vasp_ncl
cpu_medium_system.slurm   # For medium systems
cpu_large_system.slurm    # For large systems
```

**INCAR settings:**
```bash
LSORBIT = True   # SOC + non-collinear
NCORE = 8        # Good balance for SOC
KPAR = 4         # Moderate k-point parallelization
```

---

### MD (Molecular Dynamics)

```bash
# Long simulations → use large time allocation
cpu_large_system.slurm (modify time to 7-00:00:00)
```

**INCAR settings:**
```bash
NSW = 1000       # Many MD steps
NCORE = 16       # Higher NCORE for fewer k-points in MD
KPAR = 1-2       # Usually 1 k-point (Gamma) for MD
```

## Modifying Scripts

### Increase Time Limit

```bash
# Change this line in the script:
#SBATCH -t 2-00:00:00    # 2 days

# To (example):
#SBATCH -t 7-00:00:00    # 7 days
```

### Adjust MPI Ranks and Threads

**Formula:**
```
Total cores = Nodes × cores_per_node
MPI ranks × OMP_NUM_THREADS ≈ Total cores

For Perlmutter CPU:
cores_per_node = 128
```

**Example:** Change from 2 nodes to 4 nodes
```bash
# Original (2 nodes):
#SBATCH -N 2
srun -n 32 -c 16 --cpu_bind=cores vasp_ncl

# New (4 nodes):
#SBATCH -N 4
srun -n 64 -c 16 --cpu_bind=cores vasp_ncl

# Update INCAR:
KPAR = 8   # Increase from 4 to 8
```

### Switch Between Standard and Non-collinear

```bash
# Standard VASP (no SOC):
srun -n 64 -c 16 --cpu_bind=cores vasp_std

# Non-collinear (with SOC):
srun -n 64 -c 16 --cpu_bind=cores vasp_ncl

# Gamma-point only (fast for large systems):
srun -n 64 -c 16 --cpu_bind=cores vasp_gam
```

## Common Mistakes

### ❌ Wrong: Too many nodes for small system

```bash
# 8 nodes for 10-atom system → waste resources, poor scaling
#SBATCH -N 8
```

**✓ Fix:** Use 1-2 nodes for small systems

---

### ❌ Wrong: MPI ranks × threads ≠ available cores

```bash
# 2 nodes = 256 cores, but 128 ranks × 4 threads = 512 cores needed
#SBATCH -N 2
export OMP_NUM_THREADS=4
srun -n 128 -c 8 vasp_std
```

**✓ Fix:** Ensure `n × OMP_NUM_THREADS ≤ N × 128`

---

### ❌ Wrong: INCAR KPAR doesn't match resources

```bash
# 32 MPI ranks but KPAR=16 → only 2 ranks per k-point
KPAR = 16
```

**✓ Fix:** `KPAR ≤ number of k-points` and `KPAR ≤ MPI ranks / 4`

---

### ❌ Wrong: GPU script without GPU VASP

```bash
#SBATCH -C gpu
module load vasp/6.4.3-cpu   # Wrong module!
```

**✓ Fix:** Use `module load vasp/6.4.3-gpu` for GPU nodes

---

### ❌ Wrong: HSE on CPU with small resources

```bash
# HSE on 1 CPU node → will take weeks
#SBATCH -N 1
LHFCALC = True
```

**✓ Fix:** Use GPU or 8+ CPU nodes for HSE

## Performance Checklist

Before submitting, verify:

- [ ] Node count matches system size
- [ ] `MPI ranks × OMP_NUM_THREADS ≤ Total cores`
- [ ] INCAR `KPAR` is reasonable (≤ # k-points)
- [ ] INCAR `NCORE` matches `OMP_NUM_THREADS`
- [ ] Time limit is sufficient (check similar calculations)
- [ ] Correct VASP binary (vasp_std / vasp_ncl / vasp_gam)
- [ ] Correct module (cpu vs gpu)
- [ ] Memory optimization if needed (LREAL, LWAVE, LCHARG)

## Quick Examples

### Example 1: PdCrO₂ Magnetic Calculation (72 atoms, SOC)

```bash
# Use: cpu_medium_system.slurm
# 2 nodes, 32 MPI ranks, 8 threads

INCAR:
  LSORBIT = True
  NCORE = 8
  KPAR = 4
  LDAU = True
  LDAUU = 0 4 0
```

### Example 2: HSE Band Gap for GaAs (8 atoms)

```bash
# Use: gpu_hybrid_functional.slurm
# 2 GPU nodes, 8 GPUs

INCAR:
  LHFCALC = True
  HFSCREEN = 0.2
  AEXX = 0.25
  NCORE = 1
  KPAR = 8
  PRECFOCK = Fast
```

### Example 3: Interface System (400 atoms, PBE)

```bash
# Use: cpu_large_system.slurm
# 8 nodes, 128 MPI ranks, 8 threads

INCAR:
  LREAL = Auto    # Speed up for large system
  NCORE = 8
  KPAR = 8
  IDIPOL = 3
  LDIPOL = True
```

### Example 4: DOS Calculation (100 atoms, dense k-mesh)

```bash
# Use: cpu_medium_system.slurm
# 2 nodes, modify for dense k-points

INCAR:
  ISMEAR = -5
  NEDOS = 3001
  KPAR = 8       # Parallelize k-points
  NCORE = 4      # Lower NCORE for many k-points

KPOINTS:
  Gamma-centered
  0
  Monkhorst-Pack
  12 12 12       # Dense mesh
```
