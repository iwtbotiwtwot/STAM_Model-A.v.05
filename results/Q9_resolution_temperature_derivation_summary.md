# STAM Model-A Q9: Path 1 Derivation of Boundary Temperature

## Purpose

The Q8 script postulated `k_B T = hbar c |grad A| / (4 pi)` and showed it reproduces Schwarzschild, Unruh, and de Sitter temperatures. Q9 asks how much of that follows from STAM-native premises *alone*, without importing the standard surface-gravity-to-temperature relation.

## Premises

- **P1 (basic-quantum):** Vacuum holds unresolved A-fluctuations of energy E with lifetime tau = hbar / E.
- **P2 (basic-quantum):** Bulk density of fluctuation states is rho(E) = E^2 / (2 pi^2 hbar^3 c^3) per unit volume per unit energy.
- **P3 (STAM-native):** Resolution rate per fluctuation = c |grad A|. This is the *only* rate scale that can be built from c, |grad A|, and no extra dimensional inputs.
- **P4 (STAM-native):** At an A=1 boundary, only outward-oriented resolutions complete (no committed substrate inward to anchor to).

## Predicted outgoing-flux spectrum (Path 1)

Combining P1-P4: rate of resolved fluctuations per unit volume per unit time per unit energy is
```text
S(E) = rho(E) * min(1, hbar c |grad A| / E) * E / hbar
     ~ E^3       for E < hbar c |grad A|   (P_res saturated)
     ~ E^2 |grad A|   for E > hbar c |grad A|   (P_res falling as 1/E)
```
Both branches are POWER LAWS that grow with E. There is no exponential cutoff. The Planckian distribution at T_Q8 = hbar c |grad A| / (4 pi k_B) has a peak near k_B T_Q8 and falls exponentially above; the Path 1 spectrum has neither feature.

## Case-by-case results

| case                   |   grad_A_per_m |    knee_E_J |   knee / (hbar c |grad A|) |   knee / (k_B T_Q8) |   integral up to 1e3 E_cross |   integral up to 1e6 E_cross |   ratio (1e6 / 1e3) |
|:-----------------------|---------------:|------------:|---------------------------:|--------------------:|-----------------------------:|-----------------------------:|--------------------:|
| Solar-mass black hole  |    0.0003386   | 1.07234e-29 |                    1.00173 |             12.5881 |                 66.5453      |                  6.65454e+10 |               1e+09 |
| 10^15 kg primordial BH |    6.73295e+11 | 2.13232e-14 |                    1.00173 |             12.5881 |                  1.04038e+63 |                  1.04039e+72 |               1e+09 |


The `ratio (1e6 / 1e3)` column shows that extending the integration upper limit by a factor of 1000 changes the integrated emission by orders of magnitude. This is the diagnostic that the spectrum has no UV cutoff in Path 1's premises — the total emission rate is not finite without additional input.

## What fell out, what didn't — honest report

- **Characteristic energy scale:** YES. The Path 1 spectrum has a structural feature (the knee where probability saturates) at E = hbar c |grad A|. This matches the Q8 thermal scale up to a factor of 4 pi; both quantities track |grad A| with the same dependence.
- **Thermal Planckian shape:** NO. The Path 1 spectrum is two power laws (E^3 below the knee, E^2 above), with no exponential suppression. It is not a thermal distribution.
- **A meaningful temperature:** NOT FROM PATH 1 ALONE. Without an exponential UV cutoff, no finite mean energy and no well-defined temperature can be extracted. The first-moment estimator depends entirely on where the integration grid is truncated, which means it has no physical meaning.
- **Prefactor 1/(4 pi):** NO. Path 1 cannot fix the prefactor because it cannot fix the temperature.

## What Path 1 actually accomplished

Path 1 successfully established two pieces of structure from STAM-native premises:
1. A characteristic energy scale `hbar c |grad A|` exists at any A=1 boundary, with no help from QFT-on-curved-spacetime.
2. The spectrum has a knee at that scale, where the resolution probability saturates.

What Path 1 *fails* to establish is the exponential UV cutoff that would make the spectrum thermal. This is exactly the structural feature that Path 2 would need to provide. The natural candidate is the strong-field redshift / clock-rate relation: fluctuations with very high local energy have very small asymptotic energy, and a proper accounting of the redshift between local and asymptotic frequencies should produce the missing exponential cutoff via the SU logarithmic divergence near the boundary.

## Plain-language summary

STAM-native ingredients (uncertainty + resolution rate + asymmetric boundary) are enough to identify the right *energy scale* for thermal emission at an A=1 boundary, but they are NOT enough to produce a thermal spectrum or to fix the 1/(4 pi) prefactor. The missing ingredient is whatever produces the exponential cutoff at high energies — most likely the strong-field redshift relation, which would naturally connect to Path 2 via the SU integral.

## Generated plots

- `plots/Q9_spectrum_Solar-mass_black_hole.png`
- `plots/Q9_spectrum_1015_kg_primordial_BH.png`
- `plots/Q9_knee_vs_kT_scaling.png`
