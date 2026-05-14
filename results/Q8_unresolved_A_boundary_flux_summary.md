# STAM Model-A Q8: Unresolved-A Boundary Flux Test

## Purpose

Evaluate a single STAM-native rule for thermal emission at A=1 resolution boundaries:

```text
k_B T = (1 / (4 pi)) * hbar * c * |grad A|_boundary
```

The asymmetry producing a net outward flux follows from the resolved/unresolved framing. A > 1 is unresolved by definition: resolutions oriented inward have no committed substrate to anchor to and cannot complete. Outward resolutions can. The boundary therefore emits a net outward flux without pair production or negative-energy partners.

## Mechanism summary

1. Vacuum fluctuations everywhere are unresolved A excursions.
2. The resolution rate scale at any boundary is set by the gradient |grad A|.
3. At a resolution boundary (A=1) the inward direction lacks committed substrate, so only outward resolutions complete.
4. The characteristic temperature is k_B T = (hbar c / (4 pi)) |grad A|.
5. The 1/(4 pi) prefactor is forced by the STAM gravity bridge g = (c^2/2) grad A together with the standard surface-gravity-to-temperature relation T = hbar kappa / (2 pi k_B c), since kappa = (c^2/2) |grad A| at the boundary. It is not a free fit.

## Schwarzschild (Hawking) case

|   M_over_Msun |           Rs_m |   grad_A_per_m |    T_STAM_K |   T_Hawking_K |   T_ratio |    L_STAM_W |   L_standard_W |   L_ratio |
|--------------:|---------------:|---------------:|------------:|--------------:|----------:|------------:|---------------:|----------:|
|         1e-06 |    0.00295334  |    338.6       | 0.0617007   |   0.0617007   |         1 | 9.00761e-17 |    9.00761e-17 |         1 |
|         0.001 |    2.95334     |      0.3386    | 6.17007e-05 |   6.17007e-05 |         1 | 9.00761e-23 |    9.00761e-23 |         1 |
|         1     | 2953.34        |      0.0003386 | 6.17007e-08 |   6.17007e-08 |         1 | 9.00761e-29 |    9.00761e-29 |         1 |
|         1e+06 |    2.95334e+09 |      3.386e-10 | 6.17007e-14 |   6.17007e-14 |         1 | 9.00761e-41 |    9.00761e-41 |         1 |
|         1e+10 |    2.95334e+13 |      3.386e-14 | 6.17007e-18 |   6.17007e-18 |         1 | 9.00761e-49 |    9.00761e-49 |         1 |


## Unruh case

|   a_m_per_s2 |   grad_A_per_m |         T_STAM_K |        T_Unruh_K |   T_ratio |
|-------------:|---------------:|-----------------:|-----------------:|----------:|
|        1     |     2.2253e-17 |      4.05501e-21 |      4.05501e-21 |         1 |
|        1e+10 |     2.2253e-07 |      4.05501e-11 |      4.05501e-11 |         1 |
|        1e+20 |  2225.3        |      0.405501    |      0.405501    |         1 |
|        1e+26 |     2.2253e+09 | 405501           | 405501           |         1 |


## de Sitter case

|   H_per_s |      R_dS_m |   grad_A_per_m |    T_STAM_K |   T_dS_standard_K |   T_ratio |
|----------:|------------:|---------------:|------------:|------------------:|----------:|
|  1e-20    | 2.99792e+28 |    6.67128e-29 | 1.21566e-32 |       1.21566e-32 |         1 |
|  2.27e-18 | 1.32067e+26 |    1.51438e-26 | 2.75955e-30 |       2.75955e-30 |         1 |
|  1e-15    | 2.99792e+23 |    6.67128e-24 | 1.21566e-27 |       1.21566e-27 |         1 |
|  1e-10    | 2.99792e+18 |    6.67128e-19 | 1.21566e-22 |       1.21566e-22 |         1 |


## Interpretation

All three temperature ratios T_STAM / T_standard are 1.0 to floating-point precision. This is one rule applied to three different boundary configurations, with the prefactor fixed by the STAM gravity bridge rather than by data fitting. The Hawking-style luminosity scaling L proportional to 1/M^2 follows from combining the temperature rule with Stefan-Boltzmann thermodynamics for a thermal source at that T.

## What is and is not derived

- Derived from the rule: temperature scale at any A=1 boundary; recovery of Hawking, Unruh, and de Sitter temperatures with the standard prefactor.
- Not derived here: the Planckian spectrum itself. The rule fixes the temperature scale; treating the emission as thermal at that T is an additional step that future STAM work would have to ground in resolution-event statistics.
- Stefan-Boltzmann luminosity uses photon-only (two-helicity) sigma. Full Hawking emission summed over species is higher by a degrees-of-freedom factor. The 1/M^2 scaling is the result being checked, not the absolute multiplicative prefactor.

## Caveats and open questions

- The de Sitter check assumes A_dS(r) = (r/R_dS)^2 with R_dS = c/H. STAM has not yet committed to a cosmological form for A; the Schwarzschild and Unruh applications do not depend on this choice.
- The mechanism predicts that any boundary in A produces an asymmetric outward flux at T = hbar c |grad A| / (4 pi k_B), regardless of whether the boundary is gravitational, kinematic, or cosmological. This unification is the strongest content of the rule.
- Where STAM could differ from standard QFT-on-curved-spacetime predictions: details of the spectrum (deviations from exact Planckian), behavior near very small Rs, and the scaling for non-spherical / multi-source A configurations.

## Generated plots

- `plots/Q8_T_vs_grad_A_unified.png`
- `plots/Q8_hawking_T_comparison.png`
- `plots/Q8_hawking_L_comparison.png`
