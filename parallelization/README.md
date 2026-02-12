# VASP Parallelization Guide

Guide for setting up SLURM job scripts and VASP parallelization parameters.

## SLURM Job Script Structure

```bash
#!/bin/bash
#SBATCH -J test              # Job name
#SBATCH -A m3578             # Account/allocation
#SBATCH -N 2                 # Number of nodes
#SBATCH -C cpu               # Constraint (cpu or gpu)
#SBATCH -q regular           # Queue (regular/debug/premium)
#SBATCH -t 12:00:00          # Time limit (HH:MM:SS or D-HH:MM:SS)

module load vasp/6.4.3-cpu

# OpenMP settings
export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Run VASP
srun -n 32 -c 16 --cpu_bind=cores vasp_ncl
```

## Key SLURM Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `-N` | Number of nodes | `-N 2` |
| `-n` | Total MPI tasks | `-n 32` |
| `-c` | CPUs per task | `-c 16` |
| `-t` | Time limit | `-t 2-00:00:00` (2 days) |
| `-C` | Node type | `-C cpu` or `-C gpu` |
| `-q` | Queue | `-q regular` |

## OpenMP Environment Variables

```bash
export OMP_NUM_THREADS=8     # Threads per MPI task
export OMP_PLACES=threads    # Thread placement
export OMP_PROC_BIND=spread  # Thread binding (spread/close/master)
```

## VASP INCAR Parallelization Parameters

```bash
KPAR = 4        # k-point parallelization
NCORE = 8       # Cores per orbital (should match OMP_NUM_THREADS)
NPAR = 4        # Band parallelization (deprecated, use NCORE)
```

**Relationship:**
```
Total MPI ranks ≈ KPAR × NPAR
NPAR = Total MPI ranks / KPAR / NCORE
```

## Example 1: Magnetic Calculation (CPU)

From `magnetic-calculation/example-PdCrO2/calculation/scf_noncollinear/submit_vasp6.4.3_cpu.slurm`:

```bash
#!/bin/bash
#SBATCH -J test
#SBATCH -A m3578
#SBATCH -N 4                 # 4 nodes
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 2-00:00:00        # 2 days

module load vasp/6.4.3-cpu

# OpenMP settings
export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

free -mh
df -Th | grep shm

echo ""
cat INCAR
echo ""

echo "the start time is:"   $(date)  >> timing.log
DATE1=$(date +%s)

# 4 nodes × 16 MPI ranks/node = 64 total ranks
# Each rank uses 8 OpenMP threads
srun -n 64 -c 16 --cpu_bind=cores vasp_ncl

DATE2=$(date +%s)
echo "the end time is:"   $(date)   >> timing.log

diff=$((DATE2-DATE1))
printf "TIME COST: %d DAYS %02d:%02d:%02d" \
$((diff/86400)) $(((diff/3600)%24)) $(((diff/60)%60)) $(($diff %60)) >> timing.log
echo -e "\n\n" >> timing.log
```

**Corresponding INCAR:**
```bash
NCORE = 8       # Matches OMP_NUM_THREADS
KPAR = 4        # 4 k-point groups
# Total: 64 MPI ranks / 4 KPAR / 8 NCORE = 2 NPAR (implicit)
```

## Example 2: Band Alignment Calculation (CPU)

From `band-alignment-calculation/example-PbSe_Al2Se3_Al/submit.slurm`:

```bash
#!/bin/bash
#SBATCH -J test
#SBATCH -A m3578
#SBATCH -N 8                 # 8 nodes for large interface
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 2-00:00:00

module load vasp/6.4.3-cpu

export OMP_NUM_THREADS=1     # Pure MPI (no OpenMP)
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

free -mh
df -Th | grep shm

echo ""
cat INCAR
echo ""

echo "the start time is:"   $(date)  >> timing.log
DATE1=$(date +%s)

# 8 nodes × 16 MPI ranks/node = 128 total ranks
# No OpenMP threading (OMP_NUM_THREADS=1)
srun -n 128 -c 16 --cpu_bind=cores vasp_ncl

DATE2=$(date +%s)
echo "the end time is:"   $(date)   >> timing.log

diff=$((DATE2-DATE1))
printf "TIME COST: %d DAYS %02d:%02d:%02d" \
$((diff/86400)) $(((diff/3600)%24)) $(((diff/60)%60)) $(($diff %60)) >> timing.log
echo -e "\n\n" >> timing.log
```

**Note:** Pure MPI parallelization (no OpenMP threading)

## Example 3: HSE Calculation (CPU)

From `gw-calculations/Co2NbSn/band_GXWKGL/submit.slurm`:

```bash
#!/bin/bash
#SBATCH -J test
#SBATCH -A m3578
#SBATCH -N 4
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 2-00:00:00

module load vasp/6.4.3-cpu

export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

free -mh
df -Th | grep shm

echo ""
cat INCAR
echo ""

echo "the start time is:"   $(date)  >> timing.log
DATE1=$(date +%s)

srun -n 64 -c 16 --cpu_bind=cores vasp_ncl

DATE2=$(date +%s)
echo "the end time is:"   $(date)   >> timing.log

diff=$((DATE2-DATE1))
printf "TIME COST: %d DAYS %02d:%02d:%02d" \
$((diff/86400)) $(((diff/3600)%24)) $(((diff/60)%60)) $(($diff %60)) >> timing.log
echo -e "\n\n" >> timing.log
```

**Corresponding INCAR:**
```bash
NCORE = 8
KPAR = 8        # More k-point parallelization for HSE

# HSE settings
LHFCALC = True
HFSCREEN = 0.2
AEXX = 0.25
PRECFOCK = Normal
```

## CPU vs GPU

### CPU Job

```bash
#SBATCH -C cpu
module load vasp/6.4.3-cpu

srun -n 64 -c 16 --cpu_bind=cores vasp_ncl
```

**INCAR:**
```bash
NCORE = 8
KPAR = 4
```

### GPU Job

```bash
#SBATCH -C gpu
#SBATCH --gpus-per-node=4    # 4 A100 GPUs per node

module load vasp/6.4.3-gpu

export OMP_NUM_THREADS=16

# 1 MPI rank per GPU
srun -n 4 -c 32 --gpu-bind=single:1 vasp_std
```

**INCAR:**
```bash
NCORE = 1       # GPU handles parallelization
KPAR = 4        # 1 k-point group per GPU
NSIM = 4        # GPU-specific: simultaneous bands

# For HSE on GPU
PRECFOCK = Fast # Critical for GPU performance
```

## Calculating Resources

### For Perlmutter CPU Nodes

Each CPU node has 128 cores (2 × AMD EPYC 7763).

**Formula:**
```
Total cores = N (nodes) × 128
MPI ranks (-n) × OMP_NUM_THREADS ≤ Total cores
```

**Examples:**

| Nodes | MPI ranks | Threads | Total cores used | Configuration |
|-------|-----------|---------|------------------|---------------|
| 2     | 32        | 8       | 256              | Hybrid MPI+OpenMP |
| 4     | 64        | 8       | 512              | Hybrid MPI+OpenMP |
| 8     | 128       | 1       | 128              | Pure MPI |

### For Perlmutter GPU Nodes

Each GPU node has 4 A100 GPUs and 64 CPU cores.

**Typical configuration:**
- 1 MPI rank per GPU
- 4 MPI ranks per node
- 16 OpenMP threads per rank

## Adding Timing to Scripts

```bash
echo "the start time is:"   $(date)  >> timing.log
DATE1=$(date +%s)

srun -n 64 -c 16 --cpu_bind=cores vasp_ncl

DATE2=$(date +%s)
echo "the end time is:"   $(date)   >> timing.log

diff=$((DATE2-DATE1))
printf "TIME COST: %d DAYS %02d:%02d:%02d" \
$((diff/86400)) $(((diff/3600)%24)) $(((diff/60)%60)) $(($diff %60)) >> timing.log
echo -e "\n\n" >> timing.log
```

This creates a `timing.log` file with start/end times and total duration.

## Common Configurations

### Standard DFT Calculation
```bash
#SBATCH -N 2-4
export OMP_NUM_THREADS=8
srun -n 32-64 -c 16 --cpu_bind=cores vasp_ncl

INCAR:
  NCORE = 8
  KPAR = 4-8
```

### Non-collinear Magnetic Calculation
```bash
#SBATCH -N 4
export OMP_NUM_THREADS=8
srun -n 64 -c 16 --cpu_bind=cores vasp_ncl

INCAR:
  LSORBIT = True
  NCORE = 8
  KPAR = 4
  MAGMOM = ...  # 3 components per atom
```

### Hybrid Functional (HSE)
```bash
# CPU (slower but accessible)
#SBATCH -N 4-8
export OMP_NUM_THREADS=8
srun -n 64-128 -c 16 --cpu_bind=cores vasp_ncl

# GPU (much faster)
#SBATCH -C gpu
#SBATCH --gpus-per-node=4
srun -n 4-8 -c 32 --gpu-bind=single:1 vasp_std

INCAR:
  LHFCALC = True
  HFSCREEN = 0.2
  AEXX = 0.25
  PRECFOCK = Normal (CPU) or Fast (GPU)
  NCORE = 8 (CPU) or 1 (GPU)
  KPAR = 4-8
```

## Troubleshooting

### Out of Memory
- Increase number of nodes to distribute memory
- Reduce `NCORE` to use more MPI ranks
- Set `LREAL = Auto` for large systems
- Disable unnecessary outputs: `LWAVE = False`, `LCHARG = False`

### Poor Performance
- Check VASP output: `grep "running on" OUTCAR`
- Ensure `NCORE` matches `OMP_NUM_THREADS`
- Verify `KPAR` ≤ number of k-points
- Use `--cpu_bind=cores` with `srun`

### Job Fails Immediately
- Check module is loaded: `module list`
- Verify VASP binary exists: `which vasp_ncl`
- Check file permissions: `ls -l INCAR POSCAR`

## References

- VASP Parallelization Wiki: https://www.vasp.at/wiki/index.php/Category:Parallelization
- NERSC Perlmutter Documentation: https://docs.nersc.gov/systems/perlmutter/
