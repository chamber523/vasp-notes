#!/bin/bash
#SBATCH -J small_cpu
#SBATCH -A m3578
#SBATCH -N 1                 # 1 node for small systems
#SBATCH -C cpu
#SBATCH -q debug             # Use debug queue for testing
#SBATCH -t 00:30:00          # 30 minutes

module load vasp/6.4.3-cpu

# OpenMP settings (modest threading for small systems)
export OMP_NUM_THREADS=4
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
# 1 node × 128 cores = 128 total cores
# 32 MPI ranks × 4 threads = 128 cores used
# -c 8: allocate 8 CPUs per task (leaves headroom)
srun -n 32 -c 8 --cpu_bind=cores vasp_std

# Calculate elapsed time
DATE2=$(date +%s)
echo "End time: $(date)" >> timing.log
diff=$((DATE2-DATE1))
printf "TIME COST: %d DAYS %02d:%02d:%02d\n" \
  $((diff/86400)) $(((diff/3600)%24)) $(((diff/60)%60)) $((diff%60)) >> timing.log
echo -e "\n" >> timing.log

# Recommended INCAR settings for this configuration:
# NCORE = 4       # Matches OMP_NUM_THREADS
# KPAR = 8        # Good for systems with 8-32 k-points
