# GW Calculations (HSE Hybrid Functional)

Electronic structure calculations using HSE06 hybrid functional for accurate band gaps and DOS.

## Example: Co₂XY Heusler Compounds

```
gw-calculations/
├── Co2NbSn/                        # Example: Co₂NbSn
│   ├── POSCAR                       # Unit cell structure
│   ├── dos.preconv/                 # Pre-convergence with HSE
│   │   ├── INCAR                    # HSE SCF settings
│   │   ├── KPOINTS                  # k-mesh for SCF
│   │   ├── POSCAR
│   │   └── submit.slurm             # Job script
│   │
│   ├── dos/                         # Density of states
│   │   ├── INCAR                    # DOS calculation (ALGO=NONE)
│   │   ├── KPOINTS                  # Dense k-mesh
│   │   ├── POSCAR
│   │   └── submit.slurm
│   │
│   └── band_GXWKGL/                 # Band structure
│       ├── INCAR                    # Band calculation
│       ├── KPOINTS                  # High-symmetry path
│       ├── POSCAR
│       └── submit.slurm
│
├── Co2TiSn/                        # Additional examples
├── Co2ZrAl/
└── dos_band_plot.ipynb             # Visualization notebook
```

## Workflow

### 1. Pre-convergence Calculation

Generate charge density with HSE06:

```bash
cd Co2NbSn/dos.preconv/
sbatch submit.slurm
```

**Key INCAR settings:**
```bash
# Electronic convergence
ALGO = All
EDIFF = 1e-8
NELM = 500

# HSE06 hybrid functional
LHFCALC = True         # Enable hybrid functional
HFSCREEN = 0.2         # HSE screening parameter
AEXX = 0.25            # 25% exact exchange
PRECFOCK = Normal      # Speed up HSE

# SCF settings
ISTART = 0
ICHARG = 2
ISMEAR = 0             # Gaussian smearing
LCHARG = True          # Write CHGCAR
LWAVE = True           # Write WAVECAR
```

### 2. DOS Calculation

Calculate density of states using pre-converged charge density:

```bash
cd ../dos/
cp ../dos.preconv/CHGCAR .
sbatch submit.slurm
```

**Key INCAR settings:**
```bash
ALGO = NONE            # No SCF, read from CHGCAR
NELM = 1               # Single step
ISTART = 1             # Read WAVECAR
ICHARG = 0             # Read CHGCAR

ISMEAR = -5            # Tetrahedron method
NEDOS = 3001           # Dense DOS grid
LORBIT = 11            # Projected DOS (orbital-resolved)

LCHARG = False         # Don't write CHGCAR
LWAVE = False          # Don't write WAVECAR
```

### 3. Band Structure Calculation

Calculate band structure along high-symmetry path:

```bash
cd ../band_GXWKGL/
cp ../dos.preconv/CHGCAR .
sbatch submit.slurm
```

**KPOINTS format:**
```
High-symmetry path for FCC: Γ-X-W-K-Γ-L
40                         # Points per segment
Line-mode
reciprocal
0.0 0.0 0.0  ! Γ
0.5 0.0 0.5  ! X

0.5 0.0 0.5  ! X
0.5 0.25 0.75  ! W

0.5 0.25 0.75  ! W
0.375 0.375 0.75  ! K

0.375 0.375 0.75  ! K
0.0 0.0 0.0  ! Γ

0.0 0.0 0.0  ! Γ
0.5 0.5 0.5  ! L
```

### 4. Visualization

Use `dos_band_plot.ipynb` to:
1. Load DOSCAR and EIGENVAL files
2. Plot total and projected DOS (orbital-decomposed)
3. Plot band structure with spin-orbit coupling
4. Generate combined band+DOS figures

## Key Parameters

### HSE06 Hybrid Functional

HSE (Heyd-Scuseria-Ernzerhof) is a screened hybrid functional that provides more accurate band gaps than standard DFT:

```bash
LHFCALC = True         # Enable hybrid functional
HFSCREEN = 0.2         # Screening parameter (HSE06)
AEXX = 0.25            # 25% exact exchange (HSE06)
```

**Note:** HSE calculations are computationally expensive (~10-50× slower than PBE).

### Spin-Orbit Coupling (SOC)

```bash
LSORBIT = True         # Enable SOC
GGA_COMPAT = False     # Use correct kinetic energy
MAGMOM = 0 0 3 ...     # Initial magnetic moments (3 components per atom)
```

### Parallelization

HSE calculations benefit from k-point parallelization:

```bash
KPAR = 4-8             # Parallelize over k-points
NPAR = 4-8             # Parallelize over bands
NCORE = 1              # Cores per orbital
```

## Workflow Summary

```
1. dos.preconv/   →  Generate CHGCAR with HSE06 SCF
                     (ALGO=All, ICHARG=2, LCHARG=True)

2. dos/          →  Calculate DOS from CHGCAR
                     (ALGO=NONE, ICHARG=0, ISMEAR=-5, NEDOS=3001)

3. band_GXWKGL/  →  Calculate bands from CHGCAR
                     (ALGO=All, ICHARG=2, Line-mode KPOINTS)

4. Analysis      →  Plot DOS + bands (dos_band_plot.ipynb)
```

## Common Issues

### HSE Convergence

If HSE fails to converge:
- Increase `TIME` parameter (default 0.4)
- Adjust mixing parameters: `BMIX = 3`, `AMIN = 0.01`
- Try `ALGO = Damped` instead of `ALGO = All`

### Memory Usage

HSE calculations require significant memory:
- Use `PRECFOCK = Normal` (not `Fast`)
- Reduce `KPAR` if out-of-memory errors occur
- Consider using fewer k-points initially

## References

- HSE06: J. Chem. Phys. 118, 8207 (2003)
- VASP HSE manual: https://www.vasp.at/wiki/index.php/Hybrid_functionals
