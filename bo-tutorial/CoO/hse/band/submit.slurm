#!/bin/bash
#SBATCH -J test
#SBATCH -A m3578
#SBATCH -N 4
#SBATCH -C cpu
##SBATCH -G 8
#SBATCH -q regular
#SBATCH -t 2-00:00:00
##SBATCH -t 00:30:00


module load vasp/6.4.3-cpu

#OpenMP settings:
export OMP_NUM_THREADS=8 # this can be increased to 16, might improve things. Need to fo own testing
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

srun -n 64 -c 16 --cpu-bind=cores vasp_ncl
