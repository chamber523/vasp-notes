# VASP Parallelization Guide

Comprehensive guide for parallel VASP calculations on HPC clusters with SLURM.

## Overview

VASP supports multiple parallelization schemes:
1. **MPI parallelization** - Distribute over k-points, bands, or FFT grid
2. **OpenMP parallelization** - Shared memory threading within nodes
3. **Hybrid MPI+OpenMP** - Combine both for optimal performance

## SLURM Job Script Anatomy

### Basic Structure

```bash
#!/bin/bash
#SBATCH -J job_name          # Job name
#SBATCH -A account_id        # Account/allocation
#SBATCH -N 2                 # Number of nodes
#SBATCH -C cpu               # Constraint (cpu/gpu)
#SBATCH -q queue_name        # Queue (regular/debug/premium)
#SBATCH -t 12:00:00          # Time limit (HH:MM:SS or D-HH:MM:SS)

module load vasp/6.4.3-cpu

export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

srun -n 32 -c 16 --cpu_bind=cores vasp_ncl
```

## SLURM Parameters Explained

### Resource Allocation

| Parameter | Description | Example |
|-----------|-------------|---------|
| `-N` | Number of compute nodes | `-N 4` (4 nodes) |
| `-n` | Total MPI tasks | `-n 64` (64 MPI ranks) |
| `-c` | CPUs per task (for OpenMP) | `-c 16` (16 threads/rank) |
| `-t` | Wall time limit | `-t 2-00:00:00` (2 days) |
| `-C` | Node constraint | `-C cpu` or `-C gpu` |
| `-q` | Queue/QoS | `-q regular` |

### Node Architecture (Perlmutter Example)

**CPU nodes:**
- 2× AMD EPYC 7763 (64 cores each, 128 cores/node)
- 256 GB DDR4 memory
- Each core supports 2 hardware threads

**GPU nodes:**
- 1× AMD EPYC 7763 (64 cores)
- 4× NVIDIA A100 GPUs (40 GB each)
- 256 GB DDR4 memory + 160 GB HBM2e (GPU)

## Parallelization Strategy

### Formula: Total CPUs

```
Total CPUs = N (nodes) × cores_per_node
MPI tasks (n) × CPUs per task (c) ≤ Total CPUs
```

### Example 1: Pure MPI (No OpenMP)

```bash
#SBATCH -N 2                 # 2 nodes × 128 cores = 256 cores
srun -n 256 -c 1 vasp_std    # 256 MPI ranks, 1 thread each
```

**VASP INCAR:**
```bash
NCORE = 1                    # 1 core per orbital
KPAR = 8                     # Parallelize over 8 k-point groups
```

### Example 2: Hybrid MPI+OpenMP

```bash
#SBATCH -N 2                 # 2 nodes
export OMP_NUM_THREADS=8     # 8 OpenMP threads per MPI task
srun -n 32 -c 16 --cpu_bind=cores vasp_ncl
```
- 2 nodes × 16 MPI ranks/node = 32 MPI tasks
- Each MPI task uses 8 OpenMP threads
- `-c 16` allocates 16 CPUs (leaves room for hyperthreading)

**VASP INCAR:**
```bash
NCORE = 8                    # 8 cores per orbital (matches OMP_NUM_THREADS)
KPAR = 4                     # Parallelize over 4 k-point groups
NPAR = 8                     # NPAR = NCORE (deprecated, use NCORE)
```

### Example 3: GPU Acceleration

```bash
#SBATCH -N 1                 # 1 GPU node
#SBATCH -C gpu
#SBATCH --gpus=4             # Use 4 GPUs
#SBATCH --gpu-bind=none

module load vasp/6.4.3-gpu

export OMP_NUM_THREADS=16
srun -n 4 -c 32 --gpu-bind=single:1 vasp_std
```
- 1 GPU node with 4 A100 GPUs
- 4 MPI ranks (1 per GPU)
- 16 OpenMP threads per rank

**VASP INCAR:**
```bash
NCORE = 1                    # GPU offloads most work
KPAR = 4                     # 1 k-point group per GPU
NSIM = 4                     # GPU-specific: simultaneous bands
```

## OpenMP Environment Variables

```bash
export OMP_NUM_THREADS=8     # Number of threads per MPI rank
export OMP_PLACES=threads    # Map threads to CPU threads
export OMP_PROC_BIND=spread  # Distribute threads across cores
```

**Options for `OMP_PROC_BIND`:**
- `spread` - Distribute threads evenly (recommended)
- `close` - Pack threads tightly
- `master` - All threads on same core as master

## VASP INCAR Parallelization Parameters

### Key Parameters

| Parameter | Description | Recommendation |
|-----------|-------------|----------------|
| `KPAR` | k-point parallelization | # of k-points or # of nodes (if many k-points) |
| `NCORE` | Cores per orbital | Match `OMP_NUM_THREADS`, or cores_per_node/NPAR |
| `NPAR` | Band parallelization (deprecated) | Use `NCORE` instead |
| `NSIM` | Bands computed simultaneously (GPU) | 4-8 for GPUs |

### Relationship Between Parameters

```
Total MPI ranks = KPAR × NPAR (approximately)
NPAR ≈ Total cores / KPAR / NCORE
```

**Example:**
- 64 MPI ranks, 8 OpenMP threads each
- `KPAR = 4` → 4 k-point groups
- `NCORE = 8` → 8 cores per orbital
- `NPAR = 64 / 4 / 8 = 2` (implicit)

## Example Scripts

### CPU: Small System (< 100 atoms)

```bash
#!/bin/bash
#SBATCH -J small_job
#SBATCH -A m3578
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -t 00:30:00

module load vasp/6.4.3-cpu

export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# 1 node, 32 MPI ranks, 4 threads each
srun -n 32 -c 8 --cpu_bind=cores vasp_std
```

**INCAR:**
```bash
NCORE = 4
KPAR = 8
```

### CPU: Medium System (100-500 atoms)

```bash
#!/bin/bash
#SBATCH -J medium_job
#SBATCH -A m3578
#SBATCH -N 4
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 12:00:00

module load vasp/6.4.3-cpu

export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# 4 nodes, 64 MPI ranks, 8 threads each
srun -n 64 -c 16 --cpu_bind=cores vasp_ncl
```

**INCAR:**
```bash
NCORE = 8
KPAR = 4
```

### CPU: Large System (> 500 atoms) or Many k-points

```bash
#!/bin/bash
#SBATCH -J large_job
#SBATCH -A m3578
#SBATCH -N 8
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 2-00:00:00

module load vasp/6.4.3-cpu

export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# 8 nodes, 128 MPI ranks, 8 threads each
srun -n 128 -c 16 --cpu_bind=cores vasp_ncl
```

**INCAR:**
```bash
NCORE = 8
KPAR = 8
```

### GPU: Hybrid Functional (HSE)

```bash
#!/bin/bash
#SBATCH -J hse_gpu
#SBATCH -A m3578
#SBATCH -N 2
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 24:00:00
#SBATCH --gpus-per-node=4

module load vasp/6.4.3-gpu

export OMP_NUM_THREADS=16
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# 2 nodes × 4 GPUs = 8 MPI ranks (1 per GPU)
srun -n 8 -c 32 --gpu-bind=single:1 vasp_std
```

**INCAR:**
```bash
NCORE = 1
KPAR = 8
NSIM = 8

# HSE settings
LHFCALC = True
HFSCREEN = 0.2
AEXX = 0.25
PRECFOCK = Fast
```

## Performance Optimization

### 1. Test Scaling

Always test your parallelization on a small test case:

```bash
# Test different KPAR values
for kpar in 1 2 4 8; do
    sed -i "s/KPAR = .*/KPAR = $kpar/" INCAR
    sbatch submit.slurm
done
```

### 2. Monitor Performance

Add timing to your scripts:

```bash
echo "Start: $(date)" >> timing.log
DATE1=$(date +%s)

srun -n 64 -c 16 --cpu_bind=cores vasp_ncl

DATE2=$(date +%s)
diff=$((DATE2-DATE1))
printf "TIME: %d DAYS %02d:%02d:%02d\n" \
  $((diff/86400)) $(((diff/3600)%24)) $(((diff/60)%60)) $((diff%60)) >> timing.log
```

### 3. Check Memory Usage

```bash
free -mh                     # Available memory
df -Th | grep shm            # Shared memory
```

### 4. Read VASP Output

Check VASP's parallelization report in `OUTCAR`:

```bash
grep "running on" OUTCAR
grep "k-points in BZ" OUTCAR
grep "KPAR" OUTCAR
```

## Common Issues

### Out of Memory (OOM)

**Symptoms:** Job killed with "Out of Memory" error

**Solutions:**
- Reduce `NCORE` to use more MPI ranks (distribute memory)
- Increase number of nodes
- Set `LREAL = Auto` for large systems
- Disable unnecessary outputs: `LWAVE = False`, `LCHARG = False`

### Poor Scaling

**Symptoms:** Job doesn't speed up with more nodes

**Solutions:**
- Reduce number of nodes if system is too small
- Increase `KPAR` if you have many k-points
- Ensure `NCORE` matches `OMP_NUM_THREADS`
- For small systems: use fewer, faster nodes instead of many slow nodes

### Wrong Core Binding

**Symptoms:** Inconsistent timings, poor performance

**Solutions:**
- Always use `--cpu_bind=cores` with `srun`
- Set OpenMP environment variables correctly
- Ensure `-c` value gives enough CPUs for OpenMP threads

## Quick Reference

### Perlmutter CPU Node (128 cores)

| Configuration | `-N` | `-n` | `-c` | `OMP_NUM_THREADS` | `NCORE` | Total Cores |
|---------------|------|------|------|-------------------|---------|-------------|
| Pure MPI      | 1    | 128  | 1    | 1                 | 1       | 128         |
| Hybrid (low)  | 1    | 64   | 4    | 2                 | 2       | 128         |
| Hybrid (med)  | 1    | 32   | 8    | 4                 | 4       | 128         |
| Hybrid (high) | 1    | 16   | 16   | 8                 | 8       | 128         |

### Recommended Settings by System Size

| System | Nodes | MPI ranks | Threads | KPAR | NCORE |
|--------|-------|-----------|---------|------|-------|
| < 50 atoms | 1 | 32 | 4 | 4-8 | 4 |
| 50-200 atoms | 2-4 | 32-64 | 8 | 4-8 | 8 |
| 200-500 atoms | 4-8 | 64-128 | 8 | 4-8 | 8 |
| > 500 atoms | 8-16 | 128-256 | 8 | 8-16 | 8 |

## References

- VASP Wiki: https://www.vasp.at/wiki/index.php/Category:Parallelization
- NERSC Perlmutter Documentation: https://docs.nersc.gov/systems/perlmutter/
- OpenMP Best Practices: https://www.openmp.org/
