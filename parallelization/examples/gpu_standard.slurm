#!/bin/bash
#SBATCH -J gpu_job
#SBATCH -A m3578
#SBATCH -N 1                 # 1 GPU node
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 12:00:00
#SBATCH --gpus-per-node=4    # Use all 4 A100 GPUs

module load vasp/6.4.3-gpu

# OpenMP settings (more threads for GPU offloading)
export OMP_NUM_THREADS=16
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# GPU-specific settings
export CUDA_VISIBLE_DEVICES=0,1,2,3
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
# 1 node × 4 GPUs = 4 MPI ranks (1 per GPU)
# Each rank uses 16 OpenMP threads
# --gpu-bind=single:1 binds 1 GPU per MPI rank
srun -n 4 -c 32 --gpu-bind=single:1 vasp_std

# Calculate elapsed time
DATE2=$(date +%s)
echo "End time: $(date)" >> timing.log
diff=$((DATE2-DATE1))
printf "TIME COST: %d DAYS %02d:%02d:%02d\n" \
  $((diff/86400)) $(((diff/3600)%24)) $(((diff/60)%60)) $((diff%60)) >> timing.log
echo -e "\n" >> timing.log

# Recommended INCAR settings for GPU:
# NCORE = 1       # GPU handles most parallelization
# KPAR = 4        # 1 k-point group per GPU
# NSIM = 4        # GPU-specific: bands computed simultaneously
#
# For hybrid functionals (HSE):
# LHFCALC = True
# PRECFOCK = Fast # Fast FFT for GPUs
# NKRED = 2       # Reduce k-points for HF part (GPU memory)
