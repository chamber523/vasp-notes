#!/usr/bin/env python3
"""
POSCAR z-axis shift tool

Shifts atomic positions along the z-axis by a fractional amount of the c lattice parameter.
Useful for adjusting layered materials (e.g., changing from Pd-bottom to Cr-bottom).

Usage:
    python shift_z.py POSCAR_input shift_fraction [output_name]

Parameters:
    POSCAR_input    : Input POSCAR file
    shift_fraction  : Fractional shift along z-axis (0-1)
                      e.g., 0.5 = shift up by half cell
    output_name     : Optional output filename (default: input_shifted)

Examples:
    python shift_z.py POSCAR_primitive 0.5
    python shift_z.py POSCAR_primitive 0.5 POSCAR_Cr_bottom
"""

import sys

# Read POSCAR
with open(sys.argv[1], 'r') as f:
    lines = f.readlines()

shift_frac = float(sys.argv[2])
output = sys.argv[3] if len(sys.argv) > 3 else f"{sys.argv[1]}_shifted"

# Parse lattice parameters
scale = float(lines[1].strip())
lattice = []
for i in range(2, 5):
    lattice.append([float(x) * scale for x in lines[i].split()])

c_length = lattice[2][2]  # z-direction length
shift_angstrom = shift_frac * c_length

print(f"Cell c-axis length: {c_length:.4f} Å")
print(f"Shift fraction: {shift_frac:.4f}")
print(f"Shift distance: {shift_angstrom:.4f} Å")

# Find coordinate starting line
coord_start = 8
if lines[7].strip()[0].upper() in ['S', 'D']:  # Selective dynamics
    coord_start = 9

# Read atom counts
counts = [int(x) for x in lines[6].split()]
total_atoms = sum(counts)

# Shift atomic coordinates
with open(output, 'w') as f:
    # Write header lines
    for i in range(coord_start):
        f.write(lines[i])

    # Shift and write coordinates
    for i in range(coord_start, coord_start + total_atoms):
        coords = [float(x) for x in lines[i].split()[:3]]
        coords[2] = (coords[2] + shift_angstrom) % c_length  # Periodic boundary
        # Handle boundary values (e.g., 18.122 -> 0.000)
        if abs(coords[2] - c_length) < 1e-6:
            coords[2] = 0.0
        f.write(f"  {coords[0]:19.16f} {coords[1]:19.16f} {coords[2]:19.16f}\n")

print(f"\nOutput file: {output}")
