#!/bin/bash
#SBATCH -J test
#SBATCH -A m3578
#SBATCH -N 4
#SBATCH -C cpu
##SBATCH -G 8
#SBATCH -q regular
##SBATCH -t 2-00:00:00
#SBATCH -t 02:30:00
#SBATCH -o BO_%j.out


module load vasp/6.4.3-cpu

#OpenMP settings:
export OMP_NUM_THREADS=8 # this can be increased to 16, might improve things. Need to fo own testing
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Name of the conda environment to activate
CONDA_ENV_NAME=bo

# Load the Conda module
source /global/common/software/nersc/pe/conda/24.10.0/Miniforge3-24.7.1-0/etc/profile.d/conda.sh || { echo "Module 'conda' could not be loaded"; exit 1; }

# Activate the desired Conda environment
conda activate "$CONDA_ENV_NAME" || { echo "Conda environment '$CONDA_ENV_NAME' could not be activated"; exit 1; }

echo ""
cat input.json
echo ""

bo_dftu

