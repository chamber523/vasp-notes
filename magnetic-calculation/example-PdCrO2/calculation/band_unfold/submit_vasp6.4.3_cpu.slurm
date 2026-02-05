#!/bin/bash
#SBATCH -J test
#SBATCH -A m3578
#SBATCH -N 4
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -t 30:00
##SBATCH -t 2-00:00:00

module load vasp/6.4.3-cpu

#OpenMP settings:
export OMP_NUM_THREADS=8 # this can be increased to 16, might improve things. Need to fo own testing
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


srun -n 64 -c 16 --cpu_bind=cores vasp_ncl


DATE2=$(date +%s)
echo "the end time is:"   $(date)   >> timing.log

diff=$((DATE2-DATE1))
printf "TIME COST: %d DAYS %02d:%02d:%02d" \
$((diff/86400)) $(((diff/3600)%24)) $(((diff/60)%60)) $(($diff %60)) >> timing.log
echo -e "\n\n" >> timing.log

