#!/usr/bin/env python3
"""
Generate magnetic moments for PdCrO2 supercell

Based on model 4 from H. Takatsu et al., Phys. Rev. B 89, 104408 (2014)

Usage:
    python generate_magmom.py [output_name]

Parameters:
    output_name : Optional output filename (default: MAGMOM.txt)

Output:
    MAGMOM.txt - Magnetic moment vectors for 18 Cr atoms (6 layers x 3 sublattices)
    Format: Mx My Mz (one line per Cr atom)
"""

import sys
import numpy as np

def generate_pdcro2_magmoms():
    """
    Generate magnetic moments for PdCrO2 6-layer structure

    Returns 18 magnetic moment vectors (6 layers x 3 sublattices A, B, C)
    Each vector is normalized to unit length
    """

    deg = np.pi / 180  # degree to radian

    # Parameters from Figure 5(d) of the reference paper
    alpha_deg = [31, 44, 31, 44, 31, 44]  # α_n for each layer
    phi_deg = [17, 16, 17, 16, 17, 16]     # φ_n for each layer
    gamma_deg = [0, 0, 0, 0, 0, 0]         # γ_n = 0 (spin plane ⊥ z-axis)

    # ξ_n = +1 (even layers), -1 (odd layers)
    xi_n = [1 if n % 2 == 0 else -1 for n in range(6)]

    # Unit vectors
    x_hat = np.array([1, 0, 0])
    y_hat = np.array([0, 1, 0])
    z_hat = np.array([0, 0, 1])

    def e_alpha(alpha_rad):
        """Direction vector in α_n plane"""
        return np.cos(alpha_rad) * x_hat + np.sin(alpha_rad) * y_hat

    def e_alpha_prime(alpha_rad):
        """Perpendicular direction in α_n plane"""
        return np.cos(alpha_rad + np.pi/2) * x_hat + np.sin(alpha_rad + np.pi/2) * y_hat

    magnetic_moments = []

    for n in range(6):
        alpha_n = alpha_deg[n] * deg
        phi_n = phi_deg[n] * deg
        gamma_n = gamma_deg[n] * deg
        xi = xi_n[n]

        # Direction vectors
        e_a = e_alpha(alpha_n)
        e_a_p = e_alpha_prime(alpha_n)

        # A_n sublattice spin direction
        S_A = (np.cos(gamma_n) * np.cos(phi_n) * z_hat +
               np.sin(phi_n) * e_a +
               np.sin(gamma_n) * np.cos(phi_n) * e_a_p)

        # B_n sublattice: φ + ξ * 120°
        phi_B = phi_n + xi * 2 * np.pi / 3
        S_B = (np.cos(gamma_n) * np.cos(phi_B) * z_hat +
               np.sin(phi_B) * e_a +
               np.sin(gamma_n) * np.cos(phi_B) * e_a_p)

        # C_n sublattice: φ - ξ * 120°
        phi_C = phi_n - xi * 2 * np.pi / 3
        S_C = (np.cos(gamma_n) * np.cos(phi_C) * z_hat +
               np.sin(phi_C) * e_a +
               np.sin(gamma_n) * np.cos(phi_C) * e_a_p)

        magnetic_moments.append(S_A)
        magnetic_moments.append(S_B)
        magnetic_moments.append(S_C)

    return np.array(magnetic_moments)

def main():
    output_file = sys.argv[1] if len(sys.argv) > 1 else "MAGMOM.txt"

    print("Generating PdCrO2 magnetic moments...")
    print("Reference: H. Takatsu et al., Phys. Rev. B 89, 104408 (2014)")

    # Generate magnetic moments
    magmoms = generate_pdcro2_magmoms()

    print(f"\nGenerated {len(magmoms)} magnetic moment vectors")
    print("Structure: 6 layers × 3 sublattices (A, B, C)")

    # Write to file
    with open(output_file, 'w') as f:
        for i, m in enumerate(magmoms):
            f.write(f"{m[0]:12.6f} {m[1]:12.6f} {m[2]:12.6f}\n")

    print(f"\nOutput: {output_file}")
    print("\nFirst 3 moments (Layer 0: A, B, C):")
    for i in range(3):
        print(f"  {i}: [{magmoms[i][0]:7.4f}, {magmoms[i][1]:7.4f}, {magmoms[i][2]:7.4f}]")

    print("\nLast 3 moments (Layer 5: A, B, C):")
    for i in range(15, 18):
        print(f"  {i}: [{magmoms[i][0]:7.4f}, {magmoms[i][1]:7.4f}, {magmoms[i][2]:7.4f}]")

    print("\nUse with VASP:")
    print("  Copy these vectors to INCAR as:")
    print("  MAGMOM = 36*0 18*<insert_values> 18*0")
    print("  (36 O atoms, 18 Cr atoms, 18 Pd atoms)")

if __name__ == "__main__":
    main()
