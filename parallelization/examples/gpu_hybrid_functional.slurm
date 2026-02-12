#!/bin/bash
#SBATCH -J hse_gpu
#SBATCH -A m3578
#SBATCH -N 2                 # 2 GPU nodes for HSE
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 24:00:00
#SBATCH --gpus-per-node=4    # 4 GPUs per node = 8 total

module load vasp/6.4.3-gpu

# OpenMP settings (HSE benefits from more threads)
export OMP_NUM_THREADS=16
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# GPU-specific settings
export MPICH_GPU_SUPPORT_ENABLED=1

# Check available resources
free -mh
nvidia-smi

echo ""
echo "=== INCAR Settings ==="
cat INCAR
echo ""

# Timing
echo "Start time: $(date)" >> timing.log
DATE1=$(date +%s)

# Run VASP
# 2 nodes × 4 GPUs/node = 8 MPI ranks total
# Each rank uses 16 OpenMP threads and 1 GPU
srun -n 8 -c 32 --gpu-bind=single:1 vasp_std

# Calculate elapsed time
DATE2=$(date +%s)
echo "End time: $(date)" >> timing.log
diff=$((DATE2-DATE1))
printf "TIME COST: %d DAYS %02d:%02d:%02d\n" \
  $((diff/86400)) $(((diff/3600)%24)) $(((diff/60)%60)) $((diff%60)) >> timing.log
echo -e "\n" >> timing.log

# Recommended INCAR settings for HSE on GPU:
# NCORE = 1
# KPAR = 8        # 1 k-point group per GPU
# NSIM = 8        # More simultaneous bands for HSE
#
# HSE06 settings:
# LHFCALC = True
# HFSCREEN = 0.2
# AEXX = 0.25
# PRECFOCK = Fast # Critical for GPU performance
# TIME = 0.4
# NKRED = 2       # Memory optimization for GPU
#
# If GPU memory errors occur:
# NKREDX = 2      # Reduce k-points in X direction
# NKREDY = 2      # Reduce k-points in Y direction
# NKREDZ = 1      # Keep Z direction full
