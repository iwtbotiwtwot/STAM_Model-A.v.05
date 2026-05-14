# G120a -- GR Kerr QNM calibration (no STAM)

## Purpose
Establish a verified GR Kerr QNM reference table for (l=2, m=2, n=0) at the five target spins (a = 0, 0.3, 0.5, 0.7, 0.9).  No STAM substitution.  This is the calibration baseline that G120b's standalone solver must match to < 0.5% before STAM extension.

## Reference values
Source: qnm package (version 0.4.4, based on Leaver continued fraction + Berti's verified tabulation).

| a | Re(omega) | Im(omega) | |omega| | A_lm |
|---:|---:|---:|---:|---|
| 0.00 | 0.3736716844 | -0.0889623157 | 0.38411564 | 4.000000 |
| 0.30 | 0.4195266818 | -0.0877292719 | 0.42860129 | 3.653137 + 0.075127i |
| 0.50 | 0.4641230260 | -0.0856388350 | 0.47195783 | 3.342260 + 0.129196i |
| 0.70 | 0.5326002436 | -0.0807928732 | 0.53869333 | 2.903167 + 0.183158i |
| 0.90 | 0.6716142721 | -0.0648692359 | 0.67473976 | 2.109820 + 0.211124i |

**Passed: 5/5 spins.**

## Spin trends
- Re(omega) monotonically increases with spin (0.3737 at a=0.00 → 0.6716 at a=0.90).
- |Im(omega)| monotonically decreases with spin (0.0890 at a=0.00 → 0.0649 at a=0.90).

These trends match the expected Kerr GR behavior (faster prograde orbit at higher spin; less efficient ringdown decay due to longer-lived photon orbits near extremal).

## Verdict
PASS  --  qnm successfully provides verified Kerr (l=m=2, n=0) QNMs
  at all 5 target spins.  This is the GR calibration baseline.
  
  G120b's standalone Kerr QNM solver must reproduce these values to
  < 0.5% before any STAM substitution can be trusted.

## Files
- [scripts/G120a_gr_kerr_qnm_calibration.py](../scripts/G120a_gr_kerr_qnm_calibration.py)
- Reference JSON: `results/G120a_gr_kerr_qnm_reference.json`
- Reference CSV:  `results/G120a_gr_kerr_qnm_reference.csv`
- Plot: `plots/G120a_gr_kerr_qnm_calibration.png`

## Note on the < 0.5% gate
Because qnm IS the published reference (it uses Berti's verified tabulation), the 'error' against published values is intrinsically at the floating-point level.  The < 0.5% gate is what G120b's standalone Leaver / Sasaki-Nakamura / Detweiler solver must hit when calibrated against the G120a table.  This script establishes the table; G120b builds and tests the standalone solver.
