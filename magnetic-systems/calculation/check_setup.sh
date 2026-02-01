#!/bin/bash
# ============================================================
# Quick Setup Checker for VASP Calculations
# Validates INCAR, KPOINTS, POSCAR, POTCAR before submission
# ============================================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "======================================================"
echo "VASP Calculation Setup Checker"
echo "======================================================"
echo ""

# ==================== Check Files Exist ====================
echo "1. Checking required files..."
all_files_exist=true

for file in INCAR POSCAR POTCAR KPOINTS; do
    if [ -f "$file" ]; then
        echo -e "  ${GREEN}✓${NC} $file exists"
    else
        echo -e "  ${RED}✗${NC} $file MISSING!"
        all_files_exist=false
    fi
done

if [ "$all_files_exist" = false ]; then
    echo -e "\n${RED}ERROR: Missing required files!${NC}"
    exit 1
fi

echo ""

# ==================== Check POSCAR ====================
echo "2. Checking POSCAR..."

# Get atom counts
natoms_line=$(sed -n '7p' POSCAR)
natoms_total=0
for n in $natoms_line; do
    natoms_total=$((natoms_total + n))
done

echo "  Total atoms: $natoms_total"
echo "  Composition: $(sed -n '6p' POSCAR): $natoms_line"

# Check if coordinates match
coord_lines=$(tail -n +9 POSCAR | grep -c "^[[:space:]]*[0-9]")
if [ "$coord_lines" -eq "$natoms_total" ]; then
    echo -e "  ${GREEN}✓${NC} Coordinate count matches"
else
    echo -e "  ${YELLOW}⚠${NC}  WARNING: $coord_lines coordinates but $natoms_total atoms"
fi

echo ""

# ==================== Check POTCAR ====================
echo "3. Checking POTCAR..."

# Count POTCARs
npotcar=$(grep -c "TITEL" POTCAR)
nelements=$(sed -n '6p' POSCAR | wc -w)

echo "  Elements in POSCAR: $nelements"
echo "  POTCARs found: $npotcar"

if [ "$npotcar" -eq "$nelements" ]; then
    echo -e "  ${GREEN}✓${NC} POTCAR count matches"
    echo "  POTCAR elements:"
    grep "TITEL" POTCAR | sed 's/^/    /'
else
    echo -e "  ${RED}✗${NC} POTCAR count mismatch!"
fi

echo ""

# ==================== Check INCAR ====================
echo "4. Checking INCAR settings..."

# SOC check
if grep -q "LSORBIT.*=.*True" INCAR; then
    echo -e "  ${GREEN}✓${NC} SOC enabled (LSORBIT=True)"
else
    echo -e "  ${YELLOW}⚠${NC}  SOC not enabled"
fi

# DFT+U check
if grep -q "LDAU.*=.*True" INCAR; then
    echo -e "  ${GREEN}✓${NC} DFT+U enabled"
    u_values=$(grep "LDAUU" INCAR | sed 's/.*=//')
    echo "    U values: $u_values"
else
    echo -e "  ${YELLOW}⚠${NC}  DFT+U not enabled"
fi

# MAGMOM check
if grep -q "MAGMOM" INCAR; then
    magmom_line=$(grep -A 100 "MAGMOM" INCAR | grep -v "^#" | tr '\n' ' ' | tr '\' ' ')
    magmom_count=$(echo "$magmom_line" | grep -o "[0-9]\+\*[0-9]" | sed 's/\*/×/' | tr '\n' ' ')
    echo "  MAGMOM initialization: $magmom_count"

    # Count actual moment values
    moment_count=$(echo "$magmom_line" | tr ' ' '\n' | grep -c "^-\?[0-9]")
    expected_moments=$((natoms_total * 3))  # 3 components per atom for SOC

    if [ "$moment_count" -ge "$natoms_total" ]; then
        echo -e "  ${GREEN}✓${NC} MAGMOM count reasonable"
    else
        echo -e "  ${YELLOW}⚠${NC}  Check MAGMOM (found $moment_count values)"
    fi
else
    echo -e "  ${YELLOW}⚠${NC}  No MAGMOM specified"
fi

# Convergence
encut=$(grep "ENCUT" INCAR | sed 's/.*=//' | awk '{print $1}')
ediff=$(grep "EDIFF" INCAR | sed 's/.*=//' | awk '{print $1}')
echo "  ENCUT = $encut eV"
echo "  EDIFF = $ediff eV"

echo ""

# ==================== Check KPOINTS ====================
echo "5. Checking KPOINTS..."

kpoint_type=$(sed -n '3p' KPOINTS | tr '[:lower:]' '[:upper:]')
kpoint_mesh=$(sed -n '4p' KPOINTS)

echo "  Type: $kpoint_type"
echo "  Mesh: $kpoint_mesh"

# Estimate total k-points for Monkhorst-Pack
if [[ $kpoint_type == *"MONKHORST"* ]]; then
    read -r k1 k2 k3 <<< "$kpoint_mesh"
    total_kpts=$((k1 * k2 * k3))
    echo "  Total k-points: $total_kpts"

    if [ "$total_kpts" -lt 8 ]; then
        echo -e "  ${YELLOW}⚠${NC}  Very coarse mesh - check convergence"
    elif [ "$total_kpts" -gt 1000 ]; then
        echo -e "  ${YELLOW}⚠${NC}  Very dense mesh - will be slow"
    else
        echo -e "  ${GREEN}✓${NC} Reasonable mesh"
    fi
fi

echo ""

# ==================== Estimate Resources ====================
echo "6. Resource estimation..."

# Estimate memory (rough)
mem_per_atom_mb=5  # MB per atom per core (very rough)
estimated_mem_gb=$(echo "scale=1; $natoms_total * $mem_per_atom_mb * 64 / 1024" | bc)

echo "  Atoms: $natoms_total"
echo "  Estimated memory (64 cores): ~${estimated_mem_gb} GB"

if [ -f "OSZICAR" ]; then
    last_time=$(grep "LOOP+" OSZICAR | tail -1 | awk '{print $(NF-1)}')
    echo "  Previous run time: $last_time sec"
fi

echo ""

# ==================== Summary ====================
echo "======================================================"
echo "Setup Check Complete"
echo "======================================================"
echo ""
echo "Next steps:"
echo "  1. Review warnings above (if any)"
echo "  2. Copy POSCAR from ../model/ if needed"
echo "  3. Generate POTCAR (see POTCAR_README.txt)"
echo "  4. Submit job: sbatch submit.slurm"
echo ""
