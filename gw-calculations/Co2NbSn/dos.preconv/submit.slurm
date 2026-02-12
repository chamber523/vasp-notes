#!/bin/bash
#SBATCH -J test
#SBATCH -A m3578
#SBATCH -N 4
#SBATCH -C gpu
#SBATCH -G 16
#SBATCH -q regular
#SBATCH --exclusive
#SBATCH -t 1-00:00:00
##SBATCH -t 00:30:00

module load vasp/6.4.3-gpu

export NCCL_IGNORE_CPU_AFFINITY=1

#OpenMP settings:
export OMP_NUM_THREADS=1 # this can be increased to 16, might improve things. Need to fo own testing
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

free -mh
df -Th | grep shm
#ulimit -s unlimited

echo ""
cat INCAR
echo ""

echo "the start time is:"   $(date)  >> timing.log
DATE1=$(date +%s)


srun -n 16 -c 32 --cpu_bind=cores --gpu-bind=none -G 16 vasp_ncl


DATE2=$(date +%s)
echo "the end time is:"   $(date)   >> timing.log

diff=$((DATE2-DATE1))
printf "TIME COST: %d DAYS %02d:%02d:%02d" \
$((diff/86400)) $(((diff/3600)%24)) $(((diff/60)%60)) $(($diff %60)) >> timing.log
echo -e "\n\n" >> timing.log

bash /global/homes/z/zefengc/soft/DeltaSpin/scripts/energy_force.sh
