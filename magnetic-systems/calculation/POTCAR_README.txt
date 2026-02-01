POTCAR Generation Instructions
==============================

POTCAR contains pseudopotentials for O, Cr, Pd (in this order)

Method 1: Using VASP Pseudopotential Library
---------------------------------------------
# Set VASP PP path (adjust to your installation)
export VASP_PP_PATH=/path/to/vasp/pseudopotentials/potpaw_PBE

# Generate POTCAR
cat $VASP_PP_PATH/O/POTCAR \
    $VASP_PP_PATH/Cr/POTCAR \
    $VASP_PP_PATH/Pd/POTCAR > POTCAR

Method 2: Using vaspkit
-----------------------
# In directory with POSCAR
vaspkit
# Select: 103 (Generate POTCAR)
# Select: 1 (PBE)

Method 3: Using pymatgen
------------------------
from pymatgen.io.vasp.inputs import Potcar
potcar = Potcar(['O', 'Cr', 'Pd'], functional='PBE')
potcar.write_file('POTCAR')

Important Notes:
----------------
1. Element order MUST match POSCAR:
   POSCAR line 6: O Cr Pd
   POSCAR line 7: 36 18 18
   → POTCAR order: O, Cr, Pd

2. Use PBE pseudopotentials (PAW)
   - O: O (standard)
   - Cr: Cr_pv (includes 3p as valence, recommended)
   - Pd: Pd (standard)

3. Alternative Cr pseudopotentials:
   - Cr: 6 valence electrons (3d5 4s1)
   - Cr_pv: 12 valence electrons (3p6 3d5 4s1) ← RECOMMENDED
   - Cr_sv: 24 valence electrons (3s2 3p6 3d5 4s1, overkill)

4. Check POTCAR after generation:
   grep TITEL POTCAR
   # Should show:
   #   TITEL  = PAW_PBE O 08Apr2002
   #   TITEL  = PAW_PBE Cr_pv 07Sep2000
   #   TITEL  = PAW_PBE Pd 04Jan2005

5. POTCAR cannot be distributed due to VASP license
   - Each user must generate from their own VASP installation
   - Never commit POTCAR to git repositories

6. Verification:
   grep ENMAX POTCAR
   # Should show cutoff energies for each element
   # O: ~400 eV
   # Cr_pv: ~270 eV
   # Pd: ~250 eV
   # Use ENCUT = 1.3 × max(ENMAX) ≈ 520 eV (we use 400 eV for speed)
