#!/usr/bin/env python3
"""
Generate magnetic CIF (mCIF) file for VESTA visualization

Universal script for combining POSCAR structure with magnetic moments.
Works for any magnetic system (transition metals, rare earths, etc.).

Usage:
    python generate_mcif.py [options]

Options:
    -p, --poscar <file>   : POSCAR file (default: POSCAR)
    -m, --magmom <file>   : MAGMOM file (default: MAGMOM.txt)
    -o, --output <file>   : Output mCIF file (default: magnetic_structure.mcif)
    -e, --element <elem>  : Magnetic element (default: auto-detect from MAGMOM count)

Examples:
    # Auto-detect magnetic element
    python generate_mcif.py -p POSCAR -m MAGMOM.txt

    # Specify magnetic element explicitly
    python generate_mcif.py -p POSCAR -m MAGMOM.txt -e Cr

    # For PdCrO2
    python generate_mcif.py -p 03_POSCAR_supercell_6x -m MAGMOM.txt -o vesta.mcif

Requirements:
    pymatgen (auto-installed if needed)
"""

import sys
import numpy as np

try:
    from pymatgen.core import Structure
except ImportError:
    print("Error: pymatgen is required")
    print("Install with: pip install pymatgen")
    sys.exit(1)

def parse_args():
    """Parse command line arguments"""
    args = {
        'poscar': 'POSCAR',
        'magmom': 'MAGMOM.txt',
        'output': 'magnetic_structure.mcif',
        'element': None  # Auto-detect if None
    }

    i = 1
    while i < len(sys.argv):
        if sys.argv[i] in ['-p', '--poscar']:
            args['poscar'] = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] in ['-m', '--magmom']:
            args['magmom'] = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] in ['-o', '--output']:
            args['output'] = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] in ['-e', '--element']:
            args['element'] = sys.argv[i + 1]
            i += 2
        else:
            i += 1

    return args

def read_magmom(filename):
    """Read magnetic moments from file"""
    mag_flat = np.fromstring(open(filename).read(), sep=" ")
    magmoms = mag_flat.reshape(-1, 3)
    return magmoms

def assign_magmoms_to_structure(structure, magmoms, magnetic_element=None):
    """
    Assign magnetic moments to specified element atoms

    Parameters:
    -----------
    structure : pymatgen Structure
        Crystal structure
    magmoms : np.ndarray
        Magnetic moment vectors (N x 3)
    magnetic_element : str or None
        Element symbol for magnetic atoms (e.g., 'Cr', 'Fe', 'Mn')
        If None, auto-detect based on MAGMOM count

    Returns:
    --------
    structure : pymatgen Structure
        Structure with magnetic moments assigned
    """

    # Count atoms by element
    element_counts = {}
    for site in structure:
        elem = site.species_string
        element_counts[elem] = element_counts.get(elem, 0) + 1

    print(f"Structure composition: {element_counts}")

    # Auto-detect magnetic element if not specified
    if magnetic_element is None:
        # Find element whose count matches magmom count
        n_magmoms = len(magmoms)
        candidates = [elem for elem, count in element_counts.items()
                      if count == n_magmoms]

        if len(candidates) == 0:
            raise ValueError(
                f"No element found with {n_magmoms} atoms matching MAGMOM count. "
                f"Element counts: {element_counts}. "
                f"Use -e option to specify magnetic element."
            )
        elif len(candidates) > 1:
            raise ValueError(
                f"Multiple elements with {n_magmoms} atoms: {candidates}. "
                f"Use -e option to specify which element is magnetic."
            )

        magnetic_element = candidates[0]
        print(f"Auto-detected magnetic element: {magnetic_element}")
    else:
        print(f"Using specified magnetic element: {magnetic_element}")

    # Find magnetic element atom indices
    mag_indices = []
    for i, site in enumerate(structure):
        if site.species_string == magnetic_element:
            mag_indices.append(i)

    if len(mag_indices) != len(magmoms):
        raise ValueError(
            f"Mismatch: {len(mag_indices)} {magnetic_element} atoms "
            f"but {len(magmoms)} magnetic moments"
        )

    print(f"Assigning {len(magmoms)} magnetic moments to {magnetic_element} atoms")

    # Assign magmoms to all sites (zero for non-magnetic)
    mag_idx = 0
    for i, site in enumerate(structure):
        if i in mag_indices:
            site.properties["magmom"] = magmoms[mag_idx]
            mag_idx += 1
        else:
            site.properties["magmom"] = np.array([0.0, 0.0, 0.0])

    return structure

def main():
    args = parse_args()

    print("="*60)
    print("Generating magnetic CIF file")
    print("="*60)

    # Read POSCAR
    print(f"\nReading structure: {args['poscar']}")
    structure = Structure.from_file(args['poscar'])
    print(f"  Total atoms: {len(structure)}")

    # Read MAGMOM
    print(f"\nReading magnetic moments: {args['magmom']}")
    magmoms = read_magmom(args['magmom'])
    print(f"  Magnetic moment vectors: {len(magmoms)}")

    # Assign magnetic moments to structure
    print("\nAssigning magnetic moments to atoms...")
    structure = assign_magmoms_to_structure(structure, magmoms, args['element'])

    # Write mCIF
    print(f"\nWriting mCIF: {args['output']}")
    structure.to(fmt="mcif", filename=args['output'])

    print("\n" + "="*60)
    print(f"SUCCESS: mCIF file created: {args['output']}")
    print("="*60)
    print("\nVisualization in VESTA:")
    print(f"  vesta {args['output']}")
    print("\nIn VESTA:")
    print("  1. File → Open → Select the mCIF file")
    print("  2. Properties → Magnetic moments")
    print("  3. Enable magnetic moment display")
    print("  4. Adjust arrow size and style")

if __name__ == "__main__":
    main()
