# VASP Functional & vdW Correction INCAR Tags

VASP supports a range of van der Waals (vdW) corrections and meta-GGA functionals. This document covers the most common ones and their exact INCAR settings.

---

## Overview

| Method | Type | INCAR key | Typical use case |
|--------|------|-----------|-----------------|
| [TS](#ts-tkatchenko-scheffler) | Pairwise | `IVDW = 20` | Molecules, surfaces, general vdW |
| [DFT-D3(zero)](#dft-d3zero) | Pairwise | `IVDW = 11` | General purpose |
| [DFT-D3(BJ)](#dft-d3bj) | Pairwise | `IVDW = 12` | Preferred D3 variant |
| [MBD@rsSCS](#mbdrsscs-many-body-dispersion) | Many-body | `IVDW = 202` | Layered/molecular crystals |
| [rVV10](#rvv10) | Nonlocal | `LUSE_VDW + IVDW_NL = 2` | Nonlocal correlation |
| [r²SCAN+rVV10](#rscanrvv10) | Meta-GGA + Nonlocal | `METAGGA = R2SCAN + LUSE_VDW` | Accurate general-purpose |
| [optB86b-vdW](#optb86b-vdw) | Nonlocal | `GGA = MK + LUSE_VDW` | Adsorption, interfaces |
| [optB88-vdW](#optb88-vdw) | Nonlocal | `GGA = BO + LUSE_VDW` | Adsorption, interfaces |
| [optPBE-vdW](#optpbe-vdw) | Nonlocal | `GGA = OR + LUSE_VDW` | General |
| [vdW-DF](#vdw-df) | Nonlocal | `GGA = RE + LUSE_VDW` | Original Dion functional |
| [vdW-DF2](#vdw-df2) | Nonlocal | `GGA = ML + LUSE_VDW` | Improved Dion functional |

---

## TS (Tkatchenko-Scheffler)

Pairwise dispersion correction with charge-density-derived C₆ coefficients. Requires PAW datasets v52 or later.

```
IVDW        = 20       # Tkatchenko-Scheffler method
LVDW_EWALD  = .TRUE.  # Ewald summation for periodic systems (recommended)
```

**Optional parameters (defaults are usually sufficient):**

| Tag | Default | Description |
|-----|---------|-------------|
| `VDW_SR` | `0.94` | Scaling factor s_R for damping function |
| `VDW_S6` | `1.00` | Global scaling factor s₆ |
| `VDW_D` | `20.0` | Damping steepness |
| `VDW_RADIUS` | `50.0` Å | Cutoff radius for pairwise interactions |
| `LVDWSCS` | `.FALSE.` | Enable self-consistent screening (TS/SCS variant) |

> **Note**: `PREC = Accurate` is strongly recommended for TS calculations.

**Example INCAR (from TS → r²SCAN+rVV10 migration):**
```
# was:
IVDW       = 20
LVDW_EWALD = .TRUE.
```

---

## DFT-D3(zero)

Grimme's D3 correction with zero-damping. Simple and widely used.

```
IVDW = 11
```

**Optional parameters:**

| Tag | Default | Description |
|-----|---------|-------------|
| `VDW_S6` | functional-dependent | Dipole–dipole scaling |
| `VDW_S8` | functional-dependent | Dipole–quadrupole scaling |
| `VDW_SR` | functional-dependent | Radii scaling for dipole–dipole damping |
| `VDW_RADIUS` | `50.0` Å | Two-body interaction cutoff |
| `VDW_CNRADIUS` | `20.0` Å | Coordination number cutoff |

---

## DFT-D3(BJ)

D3 with Becke-Johnson damping. Generally preferred over zero-damping — avoids repulsive contributions at short range.

```
IVDW = 12
```

**Optional parameters:**

| Tag | Default | Description |
|-----|---------|-------------|
| `VDW_S6` | functional-dependent | Dipole–dipole scaling |
| `VDW_S8` | functional-dependent | Dipole–quadrupole scaling |
| `VDW_A1` | functional-dependent | Scaling of critical radii |
| `VDW_A2` | functional-dependent | Offset of critical radii |
| `VDW_RADIUS` | `50.0` Å | Two-body interaction cutoff |
| `VDW_CNRADIUS` | `20.0` Å | Coordination number cutoff |

> **Note**: Default parameters are set automatically for PBE, PBE0, RPBE, and revPBE. For other functionals, set `VDW_S6`, `VDW_S8`, `VDW_A1`, `VDW_A2` manually.

---

## MBD@rsSCS (Many-Body Dispersion)

Goes beyond pairwise additivity by including many-body screening and long-range correlation. Recommended for molecular crystals and layered materials.

```
IVDW = 202
```

**Optional parameters:**

| Tag | Default | Description |
|-----|---------|-------------|
| `VDW_SR` | `0.83` | Range-separation parameter β (functional-dependent) |
| `LSCSGRAD` | `.TRUE.` | Gradient computation for relaxations |
| `LVDWEXPANSION` | `.FALSE.` | Write 2–6 body contributions |

**Recommended `VDW_SR` per functional:**

| Functional | `VDW_SR` |
|-----------|---------|
| PBE | 0.83 |
| PBE0 | 0.85 |
| HSE06 | 0.85 |
| SCAN | 1.12 |

> **Note**: Use `PREC = Accurate`.

---

## rVV10

Nonlocal correlation functional (revised VV10). More accurate than pairwise methods for systems with nonlocal dispersion interactions.

```
GGA        = ML
LUSE_VDW   = .TRUE.
IVDW_NL    = 2
BPARAM     = 6.3
CPARAM     = 0.0093
LASPH      = .TRUE.
```

---

## r²SCAN+rVV10

Meta-GGA r²SCAN combined with rVV10 nonlocal correlation. Currently one of the most accurate general-purpose functional + vdW combinations in VASP. **Replaces TS** in the example below.

> **POTCAR requirement**: Meta-GGA functionals require POTCAR files that contain kinetic energy-density information. Verify before running:
> ```bash
> grep kinetic POTCAR
> ```
> The output must include `kinetic energy-density` and `mkinetic energy-density pseudized`. Almost all recent PAW potentials satisfy this, but older ones (e.g., `O_GW`) may not.

> **Expected VASP warning**: Running r²SCAN with standard PBE POTCARs will trigger:
> ```
> You enforced a specific xc type in the INCAR file but a different type was found in the POTCAR file.
> ```
> This is expected and harmless. There are no dedicated r²SCAN POTCARs — PBE POTCARs are the correct choice for Meta-GGA calculations. The xc type embedded in the POTCAR only affects the reference atomic calculation used to generate the PAW data, not the functional used in your actual run.

```
METAGGA    = R2SCAN    # Use r²SCAN meta-GGA functional
LUSE_VDW   = .TRUE.   # Enable nonlocal vdW correlation (rVV10)
BPARAM     = 11.95     # rVV10 b parameter (optimized for r²SCAN)
CPARAM     = 0.0093    # rVV10 c parameter (universal constant)
LASPH      = .TRUE.   # Required for meta-GGA
```

**Migrated INCAR (TS → r²SCAN+rVV10):**

```
# general
ALGO   = Fast
PREC   = Normal
EDIFF  = 1e-5
NELM   = 500
ENCUT  = 400
LASPH  = .TRUE.       # required for meta-GGA
BMIX   = 3
AMIN   = 0.01
SIGMA  = 0.05
ISMEAR = 0

# scf
ISTART = 0
ICHARG = 2
LCHARG = .FALSE.
LWAVE  = .FALSE.
LREAL  = Auto

# soc
LSORBIT     = .TRUE.
MAGMOM      = 561*0.0
GGA_COMPAT  = .FALSE.

# vdw — changed from TS (IVDW=20) to r²SCAN+rVV10
METAGGA  = R2SCAN
LUSE_VDW = .TRUE.
BPARAM   = 11.95
CPARAM   = 0.0093

# slab
IDIPOL = 3
DIPOL  = 0.50837 0.47477 0.52017
LDIPOL = .TRUE.
LREAL  = Auto
```

> **Note**: When switching from TS to r²SCAN+rVV10, remove `IVDW` and `LVDW_EWALD` entirely. The nonlocal correlation in r²SCAN+rVV10 is handled through `LUSE_VDW`, not `IVDW`.

---

## optB86b-vdW

Nonlocal functional optimized for adsorption energies. Generally gives slightly overbinding compared to optB88.

```
GGA      = MK
PARAM1   = 0.1234
PARAM2   = 1.0
AGGAC    = 0.0
LUSE_VDW = .TRUE.
LASPH    = .TRUE.
```

---

## optB88-vdW

```
GGA      = BO
PARAM1   = 0.1833333333
PARAM2   = 0.22
AGGAC    = 0.0
LUSE_VDW = .TRUE.
LASPH    = .TRUE.
```

---

## optPBE-vdW

```
GGA      = OR
AGGAC    = 0.0
LUSE_VDW = .TRUE.
LASPH    = .TRUE.
```

---

## vdW-DF

Original nonlocal functional by Dion et al. (2004).

```
GGA      = RE
AGGAC    = 0.0
LUSE_VDW = .TRUE.
LASPH    = .TRUE.
```

---

## vdW-DF2

Revised nonlocal functional by Lee et al. (2010).

```
GGA      = ML
ZAB_VDW  = -1.8867
AGGAC    = 0.0
LUSE_VDW = .TRUE.
LASPH    = .TRUE.
```

---

## References

- Tkatchenko & Scheffler, PRL 102, 073005 (2009)
- Grimme et al., JCP 132, 154104 (2010) — D3
- Grimme et al., JCC 32, 1456 (2011) — D3(BJ)
- Tkatchenko et al., PRL 108, 236402 (2012) — MBD
- Sabatini et al., PRB 87, 041108 (2013) — rVV10
- Ning et al., JPCL 13, 4029 (2022) — r²SCAN+rVV10
- Dion et al., PRL 92, 246401 (2004) — vdW-DF
- Klimeš et al., JPCM 22, 022201 (2010) — optB86b/optB88/optPBE
