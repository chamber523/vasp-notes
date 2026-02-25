# VASP INCAR Tags Reference

## What is INCAR?

INCAR is the central input file in VASP that controls all calculation parameters. It determines:

- **Calculation type**: structural relaxation, self-consistent field (SCF), density of states (DOS), band structure, hybrid functional (HSE), etc.
- **Algorithm selection**: electronic step iteration algorithm, charge density mixing scheme, etc.
- **Precision control**: energy cutoff, convergence criteria, k-point smearing, etc.
- **Output files**: whether to write CHGCAR, WAVECAR, PROCAR, etc.

Each line in INCAR follows the format `TAG = value  # comment`. Tags not specified will use VASP default values.

---

## Sections Overview

INCAR tags are organized by calculation type into the following sections:

| Section | Purpose |
|---------|---------|
| [General](#general) | Fundamental parameters shared across all calculation types |
| [Relaxation](#relaxation) | Structural relaxation of atomic positions and/or cell shape |
| [SCF](#scf) | Self-consistent field calculation to converge the charge density |
| [DOS](#dos) | Density of states calculation |
| [Band](#band) | Band structure calculation |
| [SOC](#soc) | Spin-orbit coupling calculation |
| [HSE](#hse) | HSE hybrid functional calculation |

---

## General

Universal parameters applicable to all calculation types.

```
ALGO = Fast     # Mixture of Davidson and RMM-DIIS algos
PREC = N        # Normal precision
EDIFF = 1e-5    # Convergence criteria for electronic convergence
NELM = 500      # Max number of electronic steps
ENCUT = 400     # Plane-wave cutoff energy (eV)
LASPH = True    # Include non-spherical contributions from gradient corrections
NBANDS = 12     # Number of bands to include in the calculation
BMIX = 3        # Mixing parameter for convergence
AMIN = 0.01     # Minimum mixing parameter for convergence
SIGMA = 0.05    # Width of smearing (eV)
```

### Tag Reference

| Tag | Common Values | Description |
|-----|---------------|-------------|
| `ALGO` | `Fast` / `Normal` / `All` | `Fast` is recommended for most systems; magnetic or SOC calculations may require `All` |
| `PREC` | `N` / `A` / `H` | Normal / Accurate / High — higher precision is slower |
| `EDIFF` | `1e-5` ~ `1e-6` | Use `1e-5` for relaxations, `1e-6` for high-accuracy calculations |
| `ENCUT` | 400–600 | Must be kept consistent across all calculations for valid energy comparisons |
| `LASPH` | `True` | Important for systems with d or f electrons |
| `SIGMA` | 0.05–0.2 | Use larger values for metals; smaller values for insulators/semiconductors |

---

## Relaxation

Structural relaxation calculation. Optimizes atomic positions and/or the unit cell shape and volume until forces and stresses fall below the specified convergence criteria.

```
ICHARG = 2      # Generate CHG* from a superposition of atomic charge densities
ISMEAR = 0      # Fermi smearing
LCHARG = False  # Does not write the CHG* files
LWAVE = False   # Does not write the WAVECAR
IBRION = 2      # Ionic relaxation algorithm (conjugate gradient)
ISIF = 3        # Relax ions, cell shape, and cell volume
NSW = 50        # Maximum number of ionic steps
EDIFFG = -0.01  # Convergence criteria for ionic relaxation (eV/Å, negative = force-based)
POTIM = 0.5     # Step size for ionic motion (scaling factor for CG/RMM, time step for MD)
```

### Tag Reference

| Tag | Common Values | Description |
|-----|---------------|-------------|
| `IBRION` | `2` / `1` / `0` | `2` = conjugate gradient (recommended for most relaxations); `1` = quasi-Newton (RMM-DIIS, good for small displacements); `0` = molecular dynamics |
| `ISIF` | `2` / `3` / `4` / `7` | Controls which degrees of freedom are relaxed (see table below) |
| `NSW` | 50–500 | Maximum ionic steps; set to `0` for a static single-point calculation |
| `EDIFFG` | `-0.01` ~ `-0.05` | Negative value = force convergence (eV/Å); positive value = energy convergence (eV) |
| `POTIM` | `0.3`–`0.5` | Step size scaling for CG/RMM; reduce if relaxation diverges |

### ISIF Values

| ISIF | Forces | Stress tensor | Ions | Cell shape | Cell volume |
|------|--------|---------------|------|------------|-------------|
| `2` | Yes | Calculated | Yes | No | No |
| `3` | Yes | Yes | Yes | Yes | Yes |
| `4` | Yes | Yes | Yes | Yes | No |
| `7` | Yes | Yes | No | No | Yes |

> **Note**: Always check the final forces and stress in OUTCAR after relaxation. If NSW is reached before convergence, restart from the last CONTCAR.

---

## SCF

Self-consistent field calculation. Generates a converged charge density (CHGCAR) used as input for DOS and band structure calculations.

```
ICHARG = 2      # Generate CHG* from a superposition of atomic charge densities
ISMEAR = 0      # Fermi smearing
LCHARG = True   # Write the CHG* files
LWAVE = False   # Does not write the WAVECAR
LREAL = Auto    # Automatically chooses real/reciprocal space for projections
```

### Tag Reference

| Tag | Common Values | Description |
|-----|---------------|-------------|
| `ICHARG` | `2` / `1` / `11` | `2` = start from atomic superposition; `1` = read CHGCAR and continue; `11` = fix charge density |
| `ISMEAR` | `0` / `1` / `-5` | `0` = Fermi (semiconductors); `1` = Methfessel-Paxton (metals); `-5` = tetrahedron (DOS) |
| `LCHARG` | `True` / `False` | Must be `True` in SCF to output CHGCAR for subsequent calculations |
| `LWAVE` | `True` / `False` | Only write if the wavefunction is needed (e.g., band unfolding) |
| `LREAL` | `Auto` / `False` | `False` for small cells (accurate); `Auto` for large cells (faster) |

---

## DOS

Density of states calculation. Reads the pre-converged CHGCAR from SCF and fixes the charge density while computing eigenvalues on a dense k-mesh.

```
ICHARG = 11     # Calculate eigenvalues from preconverged CHGCAR
ISMEAR = -5     # Tetrahedron method with Blochl corrections
LCHARG = False  # Does not write the CHG* files
LWAVE = False   # Does not write the WAVECAR files
LORBIT = 11     # Projected data (lm-decomposed PROCAR)
NEDOS = 3001    # Number of grid points sampled for the DOS
EMIN = -3.7     # Minimum energy for the DOS plot (eV)
EMAX = 10.3     # Maximum energy for the DOS plot (eV)
```

### Tag Reference

| Tag | Common Values | Description |
|-----|---------------|-------------|
| `ICHARG` | `11` | Must be `11` for DOS/band calculations to keep charge density fixed |
| `ISMEAR` | `-5` | Tetrahedron method gives the most accurate DOS lineshape; requires sufficient k-points |
| `LORBIT` | `10` / `11` | `10` = lm-decomposed without full projection; `11` = full lm-decomposed PROCAR for PDOS |
| `NEDOS` | 2000–5000 | Default is 301; increase for smoother DOS curves |
| `EMIN` / `EMAX` | system-dependent | Setting a reasonable range improves sampling density |

---

## Band

Band structure calculation. Computes eigenvalues along a high-symmetry k-path using the fixed charge density from SCF.

```
ICHARG = 11     # Calculate eigenvalues from preconverged CHGCAR
ISMEAR = 0      # Fermi smearing
LCHARG = False  # Does not write the CHG* files
LWAVE = False   # Does not write the WAVECAR files (set True for band unfolding)
LORBIT = 11     # Projected data (lm-decomposed PROCAR)
```

### Tag Reference

| Tag | Common Values | Description |
|-----|---------------|-------------|
| `ISMEAR` | `0` | Band calculations use a k-path, so the tetrahedron method (`-5`) cannot be used |
| `LWAVE` | `False` / `True` | Set to `True` when performing supercell band unfolding |
| `LORBIT` | `11` | Outputs PROCAR for fat band analysis with tools like pyprocar or vaspkit |

> **Note**: The KPOINTS file for band calculations must use line mode (high-symmetry path), which is different from the k-mesh used in SCF.

---

## SOC

Spin-orbit coupling (SOC) calculation. Required for heavy-element systems and topological materials.

```
LSORBIT = True    # Turn on spin-orbit coupling
MAGMOM = 6*0      # Magnetic moment for each atom (3 components per atom: mx my mz)
```

### Tag Reference

| Tag | Common Values | Description |
|-----|---------------|-------------|
| `LSORBIT` | `True` | Enables SOC; automatically activates non-collinear magnetism (`LNONCOLLINEAR = True`) |
| `MAGMOM` | `N*0` or `0 0 m` | Under SOC, each atom requires 3 components (x, y, z magnetic moment) |
| `SAXIS` | `0 0 1` | Spin quantization axis direction; default is the z-axis |

> **Note**: SOC calculations typically require a pre-converged non-SOC SCF calculation. Use its CHGCAR as the starting charge density (`ICHARG = 1`).

---

## HSE

HSE hybrid functional calculation. Provides more accurate band gaps than GGA by including a fraction of exact Hartree-Fock exchange.

```
LHFCALC = True    # Determines if a hybrid functional is used
HFSCREEN = 0.2    # Range-separation parameter (Å⁻¹)
AEXX = 0.25       # Fraction of exact exchange to be used
PRECFOCK = Fast   # Increases the speed of HSE calculations
```

### Tag Reference

| Tag | Common Values | Description |
|-----|---------------|-------------|
| `LHFCALC` | `True` | Must be set to enable Hartree-Fock exchange |
| `HFSCREEN` | `0.2` | Standard HSE06 value; `0` gives PBE0; larger values converge to PBE |
| `AEXX` | `0.25` | Standard HSE06 value (25%); can be tuned to match experimental band gaps |
| `PRECFOCK` | `Fast` / `Accurate` | `Fast` reduces computational cost; use `Accurate` for high-precision work |
| `ALGO` | `All` / `Damped` | `All` is recommended for HSE to avoid convergence issues |

> **Note**: HSE is far more expensive than GGA. It is best practice to first complete a GGA structural relaxation, then run HSE as a single-point calculation on the relaxed structure.

---

## Typical Calculation Workflow

```
Relaxation  (IBRION=2, ISIF=3, NSW=50)
      |
      v
SCF  (ICHARG=2, LCHARG=True)
      |
      |---> DOS  (ICHARG=11, dense k-mesh)
      |
      |---> Band (ICHARG=11, k-path line mode)
      |
      |---> SOC  (LSORBIT=True, ICHARG=1 from SCF CHGCAR)
      |
      |---> HSE  (LHFCALC=True, single-point on relaxed structure)
```
