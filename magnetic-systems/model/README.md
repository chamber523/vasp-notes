# Model Directory

Structure preparation and magnetic configuration for PdCrO₂ calculations.

## Structure Files

- `01_POSCAR_unitcell` - Primitive cell from Materials Project (Pd-bottom, 12 atoms)
- `02_POSCAR_Cr-bottom` - Unit cell with Cr at bottom layer (z-shifted by 0.5c)
- `03_POSCAR_supercell_6x` - 6× supercell (72 atoms, 3× in ab-plane, 2× along c)

## Generated Files

- `MAGMOM.txt` - Magnetic moment vectors for 18 Cr atoms (3 components per line)
- `magnetic_structure.mcif` - Magnetic CIF for VESTA visualization

## Scripts

### 1. shift_z.py
Shift atomic positions along z-axis by fractional amount.

```bash
python shift_z.py <POSCAR_input> <shift_fraction> [output_name]

# Example: shift by half cell (Pd-bottom → Cr-bottom)
python shift_z.py 01_POSCAR_unitcell 0.5 02_POSCAR_Cr-bottom
```

### 2. make_supercell.py
Generate supercell from transformation matrix.

```bash
# Interactive mode (enter matrix row by row)
python make_supercell.py <POSCAR_input> [output_name]

# Preset mode for PdCrO₂ (2 1 0 / 1 2 0 / 0 0 2)
python make_supercell.py 02_POSCAR_Cr-bottom --pdcro2
```

**Transformation matrix format:**
- Row 1: `a_M = [integer coefficients of a, b, c]`
- Row 2: `b_M = [integer coefficients of a, b, c]`
- Row 3: `c_M = [integer coefficients of a, b, c]`

**Example for 6× expansion:**
```
Row 1: 2 1 0    # a_M = 2a + b
Row 2: 1 2 0    # b_M = a + 2b
Row 3: 0 0 2    # c_M = 2c
```

### 3. generate_magmom.py
Generate magnetic moments for PdCrO₂ 6-layer structure.

```bash
python generate_magmom.py [output_name]

# Default output: MAGMOM.txt
python generate_magmom.py
```

**Theory:** Based on model 4 from H. Takatsu et al., Phys. Rev. B 89, 104408 (2014)
- 6 layers with alternating magnetic structures
- 3 sublattices (A, B, C) per layer
- Non-collinear 120° structure with canting

**Output format:** 18 lines, each with Mx My Mz components

### 4. generate_mcif.py
**Universal** magnetic CIF generator for VESTA visualization.

Works for **any magnetic system** (Cr, Fe, Mn, Co, Ni, rare earths, etc.)

```bash
python generate_mcif.py [options]

# For PdCrO2 (specify Cr as magnetic element)
python generate_mcif.py -p 03_POSCAR_supercell_6x -m MAGMOM.txt -e Cr -o vesta.mcif

# Auto-detect magnetic element (works when only one element matches MAGMOM count)
python generate_mcif.py -p POSCAR -m MAGMOM.txt

# For other magnetic systems (e.g., Fe-based compounds)
python generate_mcif.py -p POSCAR_Fe -m MAGMOM_Fe.txt -e Fe
```

**Options:**
- `-p, --poscar <file>` : POSCAR file (default: POSCAR)
- `-m, --magmom <file>` : MAGMOM file (default: MAGMOM.txt)
- `-o, --output <file>` : Output mCIF file (default: magnetic_structure.mcif)
- `-e, --element <elem>` : Magnetic element symbol (e.g., Cr, Fe, Mn)
  - Auto-detects if only one element count matches MAGMOM count
  - **Required** when multiple elements have same atom count

**Key Features:**
- **Universal:** Works for any magnetic material
- **Auto-detection:** Automatically identifies magnetic element when unambiguous
- **Flexible:** Handles collinear and non-collinear magnetic structures
- **VESTA-compatible:** Generates standard mCIF format

**Requirements:** `pymatgen` (auto-installed)

## Complete Workflow

```bash
# Step 1: Obtain and prepare primitive cell
# (Download from Materials Project → 01_POSCAR_unitcell)

# Step 2: Adjust to Cr-bottom configuration
python shift_z.py 01_POSCAR_unitcell 0.5 02_POSCAR_Cr-bottom

# Step 3: Create 6× supercell
python make_supercell.py 02_POSCAR_Cr-bottom --pdcro2

# Step 4: Generate magnetic moments
python generate_magmom.py

# Step 5: Create magnetic structure for visualization
python generate_mcif.py -p 03_POSCAR_supercell_6x -m MAGMOM.txt -e Cr

# Step 6: Visualize in VESTA
vesta magnetic_structure.mcif
```

## Output Summary

After completing the workflow:
- `03_POSCAR_supercell_6x` → Use for VASP calculations
- `MAGMOM.txt` → Use to set MAGMOM in INCAR
- `magnetic_structure.mcif` → Use for visualization and validation

## Next Steps

Move to `../calculation/` directory to:
1. Set up INCAR with magnetic parameters
2. Prepare KPOINTS for DFT calculations
3. Configure POTCAR for O, Cr, Pd
4. Run VASP calculations
