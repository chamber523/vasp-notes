#!/bin/bash
#SBATCH -J large_cpu
#SBATCH -A m3578
#SBATCH -N 8                 # 8 nodes for large systems or many k-points
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 2-00:00:00        # 2 days

module load vasp/6.4.3-cpu

# OpenMP settings (optimized for large systems)
export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Check available resources
free -mh
df -Th | grep shm
#ulimit -s unlimited         # Uncomment if stack overflow errors occur

echo ""
echo "=== INCAR Settings ==="
cat INCAR
echo ""

# Timing
echo "Start time: $(date)" >> timing.log
DATE1=$(date +%s)

# Run VASP
# 8 nodes × 128 cores/node = 1024 total cores
# 128 MPI ranks × 8 threads = 1024 cores used
# -c 16: allocate 16 CPUs per task
srun -n 128 -c 16 --cpu_bind=cores vasp_ncl

# Calculate elapsed time
DATE2=$(date +%s)
echo "End time: $(date)" >> timing.log
diff=$((DATE2-DATE1))
printf "TIME COST: %d DAYS %02d:%02d:%02d\n" \
  $((diff/86400)) $(((diff/3600)%24)) $(((diff/60)%60)) $((diff%60)) >> timing.log
echo -e "\n" >> timing.log

# Recommended INCAR settings for this configuration:
# NCORE = 8       # Matches OMP_NUM_THREADS
# KPAR = 8        # Parallelize over 8 k-point groups (if available)
#                 # For systems with many k-points, can increase to 16
#
# Memory optimization (if OOM errors):
# LREAL = Auto    # Real-space projection for large systems
# LWAVE = False   # Don't write WAVECAR
# LCHARG = False  # Don't write CHGCAR
