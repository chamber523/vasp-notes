#!/usr/bin/env python3
"""
Universal POSCAR supercell generation script (pure Python)

Interactive transformation matrix input

Usage:
    python make_supercell.py POSCAR_input [output_name]

After running, you will be prompted to enter the transformation matrix
row by row (3 integers per row).

Example:
    python make_supercell.py POSCAR_Cr_bottom

    Then input:
    Enter transformation matrix (3 lines, 3 integers each):
    Row 1: 2 1 0
    Row 2: 1 2 0
    Row 3: 0 0 2

Quick mode (preset transformation matrix):
    python make_supercell.py POSCAR_input --pdcro2
    (Automatically uses PdCrO2 standard supercell: 2 1 0 / 1 2 0 / 0 0 2)
"""

import sys

def mat_mul(A, B):
    """Matrix multiplication A @ B"""
    if isinstance(B[0], list):
        result = [[sum(A[i][k] * B[k][j] for k in range(len(B)))
                   for j in range(len(B[0]))] for i in range(len(A))]
    else:
        result = [sum(A[i][k] * B[k] for k in range(len(B))) for i in range(len(A))]
    return result

def mat_det(M):
    """Calculate 3x3 matrix determinant"""
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
            M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
            M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))

def mat_inv(M):
    """Calculate 3x3 matrix inverse"""
    det = mat_det(M)
    if abs(det) < 1e-10:
        raise ValueError("Matrix is singular, cannot invert")

    adj = [
        [M[1][1]*M[2][2] - M[1][2]*M[2][1], M[0][2]*M[2][1] - M[0][1]*M[2][2], M[0][1]*M[1][2] - M[0][2]*M[1][1]],
        [M[1][2]*M[2][0] - M[1][0]*M[2][2], M[0][0]*M[2][2] - M[0][2]*M[2][0], M[0][2]*M[1][0] - M[0][0]*M[1][2]],
        [M[1][0]*M[2][1] - M[1][1]*M[2][0], M[0][1]*M[2][0] - M[0][0]*M[2][1], M[0][0]*M[1][1] - M[0][1]*M[1][0]]
    ]
    return [[adj[i][j] / det for j in range(3)] for i in range(3)]

def mat_transpose(M):
    """Matrix transpose"""
    return [[M[j][i] for j in range(len(M))] for i in range(len(M[0]))]

def vec_add(a, b):
    """Vector addition"""
    return [a[i] + b[i] for i in range(len(a))]

def vec_norm(v):
    """Vector norm"""
    return sum(x**2 for x in v) ** 0.5

def read_poscar(filename):
    """Read POSCAR file"""
    with open(filename, 'r') as f:
        lines = f.readlines()

    comment = lines[0].strip()
    scale = float(lines[1].strip())

    lattice = []
    for i in range(2, 5):
        vec = [float(x) * scale for x in lines[i].split()]
        lattice.append(vec)

    elements = lines[5].split()
    counts = [int(x) for x in lines[6].split()]

    coord_line = 7
    if lines[coord_line].strip()[0].upper() in ['S']:
        coord_line = 8
    coord_type = lines[coord_line].strip()[0].upper()

    total_atoms = sum(counts)
    coords = []
    for i in range(coord_line + 1, coord_line + 1 + total_atoms):
        coords.append([float(x) for x in lines[i].split()[:3]])

    return {
        'comment': comment,
        'lattice': lattice,
        'elements': elements,
        'counts': counts,
        'coord_type': coord_type,
        'coords': coords
    }

def write_poscar(filename, data):
    """Write POSCAR file"""
    with open(filename, 'w') as f:
        f.write(f"{data['comment']}\n")
        f.write("   1.0\n")

        for vec in data['lattice']:
            f.write(f"    {vec[0]:19.16f} {vec[1]:19.16f} {vec[2]:19.16f}\n")

        f.write("   " + "   ".join(data['elements']) + "\n")
        f.write("   " + "   ".join(str(c) for c in data['counts']) + "\n")

        f.write(f"{data['coord_type']}artesian\n")
        for coord in data['coords']:
            f.write(f"  {coord[0]:19.16f} {coord[1]:19.16f} {coord[2]:19.16f}\n")

def make_supercell(poscar_data, P):
    """Generate supercell"""
    P_T = mat_transpose(P)
    new_lattice = mat_mul(P_T, poscar_data['lattice'])

    n_cells = int(round(abs(mat_det(P))))

    print(f"\nSupercell size: {n_cells}x unit cell")
    print(f"Atoms: {sum(poscar_data['counts'])} -> {sum(poscar_data['counts']) * n_cells}")

    P_inv = mat_inv(P)
    translations = []

    max_val = max(max(abs(x) for x in row) for row in P)
    max_range = int(max_val) + 2

    for i in range(-max_range, max_range + 1):
        for j in range(-max_range, max_range + 1):
            for k in range(-max_range, max_range + 1):
                t = [float(i), float(j), float(k)]
                frac = mat_mul(P_inv, t)
                if all(-1e-10 <= f < 1.0 - 1e-10 for f in frac):
                    translations.append(t)

    print(f"Generated {len(translations)} translation vectors")

    lattice_T = mat_transpose(poscar_data['lattice'])
    new_coords = []

    for coord in poscar_data['coords']:
        for t in translations:
            shift = mat_mul(lattice_T, t)
            new_coord = vec_add(coord, shift)
            new_coords.append(new_coord)

    new_counts = [c * n_cells for c in poscar_data['counts']]

    return {
        'comment': poscar_data['comment'] + f" (supercell {n_cells}x)",
        'lattice': new_lattice,
        'elements': poscar_data['elements'],
        'counts': new_counts,
        'coord_type': poscar_data['coord_type'],
        'coords': new_coords
    }

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    input_file = sys.argv[1]
    poscar_data = read_poscar(input_file)

    print(f"Read structure: {input_file}")
    print(f"Atoms: {sum(poscar_data['counts'])}")
    print(f"Elements: {' '.join(f'{e}({n})' for e, n in zip(poscar_data['elements'], poscar_data['counts']))}")

    # Display original lattice
    print(f"\nOriginal lattice vectors:")
    for i, vec in enumerate(poscar_data['lattice']):
        label = ['a', 'b', 'c'][i]
        length = vec_norm(vec)
        print(f"  {label} = [{vec[0]:8.4f}, {vec[1]:8.4f}, {vec[2]:8.4f}] Å  (|{label}| = {length:.4f} Å)")

    # Check for preset mode
    if len(sys.argv) > 2 and sys.argv[2] == '--pdcro2':
        print("\nUsing PdCrO2 preset mode")
        P = [[2, 1, 0], [1, 2, 0], [0, 0, 2]]
        print("Transformation matrix:")
        for i, row in enumerate(P):
            print(f"  Row {i+1}: {row[0]:2d} {row[1]:2d} {row[2]:2d}")
        output = input_file + '_supercell_6x'
    else:
        # Interactive transformation matrix input
        print("\nEnter transformation matrix (Transformation Matrix):")
        print("Enter 3 integers per row, separated by spaces")
        print("Example: 2 1 0")
        print()

        P = []
        for i in range(3):
            while True:
                try:
                    line = input(f"Row {i+1}: ").strip()
                    row = [int(x) for x in line.split()]
                    if len(row) != 3:
                        print("  Error: Need exactly 3 integers, please re-enter")
                        continue
                    P.append(row)
                    break
                except ValueError:
                    print("  Error: Please enter integers, please re-enter")
                except EOFError:
                    print("\nInput interrupted")
                    sys.exit(1)

        # Generate output filename
        det = abs(mat_det(P))
        output = input_file + f'_supercell_{int(det)}x'

    # Custom output name
    if len(sys.argv) > 2 and sys.argv[-1] not in ['--pdcro2']:
        output = sys.argv[-1]

    # Display transformation matrix info
    print(f"\nTransformation matrix P:")
    for row in P:
        print(f"  [{row[0]:2d}, {row[1]:2d}, {row[2]:2d}]")

    # Generate supercell
    supercell = make_supercell(poscar_data, P)
    write_poscar(output, supercell)

    print(f"\nOutput file: {output}")
    print(f"\nNew lattice vectors:")
    for i, vec in enumerate(supercell['lattice']):
        label = ['a', 'b', 'c'][i]
        length = vec_norm(vec)
        print(f"  {label}_M = [{vec[0]:8.4f}, {vec[1]:8.4f}, {vec[2]:8.4f}] Å  (|{label}_M| = {length:.4f} Å)")

    print("\n✓ Supercell generation complete!")

if __name__ == "__main__":
    main()
