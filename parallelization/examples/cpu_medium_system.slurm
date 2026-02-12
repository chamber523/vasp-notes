#!/bin/bash
#SBATCH -J medium_cpu
#SBATCH -A m3578
#SBATCH -N 2                 # 2 nodes for medium systems
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 12:00:00          # 12 hours

module load vasp/6.4.3-cpu

# OpenMP settings (balanced MPI+OpenMP)
export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Check available resources
free -mh
df -Th | grep shm

echo ""
echo "=== INCAR Settings ==="
cat INCAR
echo ""

# Timing
echo "Start time: $(date)" >> timing.log
DATE1=$(date +%s)

# Run VASP
# 2 nodes × 128 cores/node = 256 total cores
# 32 MPI ranks × 8 threads = 256 cores used
# -c 16: allocate 16 CPUs per task (2× threads for hyperthreading)
srun -n 32 -c 16 --cpu_bind=cores vasp_ncl

# Calculate elapsed time
DATE2=$(date +%s)
echo "End time: $(date)" >> timing.log
diff=$((DATE2-DATE1))
printf "TIME COST: %d DAYS %02d:%02d:%02d\n" \
  $((diff/86400)) $(((diff/3600)%24)) $(((diff/60)%60)) $((diff%60)) >> timing.log
echo -e "\n" >> timing.log

# Recommended INCAR settings for this configuration:
# NCORE = 8       # Matches OMP_NUM_THREADS
# KPAR = 4        # Good for systems with 4-16 k-points
# NPAR = 8        # Automatically calculated as total_ranks/KPAR/NCORE
