#!/bin/bash
#SBATCH --job-name=pdcro2_scf       # Job name
#SBATCH --partition=normal          # Partition/queue name
#SBATCH --nodes=2                   # Number of nodes
#SBATCH --ntasks-per-node=32        # Tasks per node
#SBATCH --cpus-per-task=1           # CPUs per task
#SBATCH --time=24:00:00             # Max runtime (HH:MM:SS)
#SBATCH --output=vasp_%j.out        # Standard output file
#SBATCH --error=vasp_%j.err         # Standard error file
#SBATCH --mail-type=END,FAIL        # Email notifications
#SBATCH --mail-user=your@email.com  # Your email

# ============================================================
# VASP Job Submission Script for PdCrO2 Calculations
# Adjust parameters based on your cluster configuration
# ============================================================

# ==================== Environment Setup ====================
# Load required modules (adjust for your cluster)
module purge
module load intel/2021.4              # Intel compiler
module load impi/2021.4               # Intel MPI
module load vasp/6.4.3                # VASP (adjust version)

# Or for GNU environment:
# module load gcc/11.2.0
# module load openmpi/4.1.1
# module load vasp/6.4.3-gnu

# Set number of threads (optional, for hybrid MPI+OpenMP)
export OMP_NUM_THREADS=1

# ==================== Job Information ====================
echo "======================================================"
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Total tasks: $SLURM_NTASKS"
echo "Working directory: $(pwd)"
echo "======================================================"
echo ""

# ==================== Pre-run Checks ====================
# Check required files exist
echo "Checking input files..."
for file in INCAR POSCAR POTCAR KPOINTS; do
    if [ ! -f "$file" ]; then
        echo "ERROR: $file not found!"
        exit 1
    fi
    echo "  ✓ $file found ($(wc -l < $file) lines)"
done
echo ""

# Display INCAR settings
echo "Key INCAR parameters:"
grep -E "^[[:space:]]*(ENCUT|EDIFF|LSORBIT|LDAU)" INCAR | sed 's/^/  /'
echo ""

# Display KPOINTS
echo "K-point mesh:"
head -4 KPOINTS | tail -1 | sed 's/^/  /'
echo ""

# ==================== Run VASP ====================
echo "Starting VASP calculation..."
echo "======================================================"

# Choose VASP executable:
# - vasp_std: Standard version (most common)
# - vasp_gam: Gamma-point only (faster for Gamma-only)
# - vasp_ncl: Non-collinear + SOC (use if std doesn't support SOC)

# Method 1: Using srun (SLURM native)
srun vasp_std

# Method 2: Using mpirun (alternative)
# mpirun -np $SLURM_NTASKS vasp_std

# Method 3: Specific for Intel MPI
# mpiexec.hydra -np $SLURM_NTASKS vasp_std

echo "======================================================"
echo "VASP finished at: $(date)"
echo ""

# ==================== Post-run Analysis ====================
echo "Post-run checks..."

# Check if calculation completed
if grep -q "reached required accuracy" OUTCAR 2>/dev/null; then
    echo "  ✓ Calculation CONVERGED"
else
    echo "  ✗ Calculation DID NOT CONVERGE or crashed"
fi

# Extract key results
if [ -f OUTCAR ]; then
    echo ""
    echo "Final energy:"
    grep "energy  without entropy" OUTCAR | tail -1 | sed 's/^/  /'

    echo ""
    echo "Final magnetic moments (first 3 Cr atoms):"
    grep "magnetization (x,y,z)" OUTCAR | head -40 | tail -3 | sed 's/^/  /'

    echo ""
    echo "Total runtime:"
    grep "LOOP+:" OUTCAR | tail -1 | sed 's/^/  /'
fi

# Check for warnings
if [ -f OUTCAR ]; then
    echo ""
    echo "Warnings (if any):"
    grep -i "WARNING\|ERROR" OUTCAR | head -5 | sed 's/^/  /'
fi

# File sizes
echo ""
echo "Output file sizes:"
ls -lh OUTCAR CHGCAR WAVECAR 2>/dev/null | awk '{print "  " $9 ": " $5}'

echo ""
echo "======================================================"
echo "Job completed at: $(date)"
echo "======================================================"

# ==================== Optional: Cleanup ====================
# Uncomment to save disk space after confirming results

# Remove large files if converged
# if grep -q "reached required accuracy" OUTCAR; then
#     echo "Removing large temporary files..."
#     rm -f WAVECAR CHG  # Keep CHGCAR for DOS/band
# fi

# Compress output
# gzip OUTCAR vasprun.xml

# ==================== Exit ====================
exit 0
