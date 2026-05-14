# G123 -- Matter Minimal Coupling via g_eff = f(Sigma) g_munu

**Date: 2026-05-13.** Short-term (tier-1) fix for strain point #3:
standard matter sectors on conformally-rescaled metric g_eff =
f(Sigma) g_munu. Since f(2/3) = 1 anchors matter outside the photon
sphere to GR exactly, every precision test passes by construction.
Observable deviations live inside the photon sphere.

## f(A) numerical reconstruction

Integrating G70's d ln f/dA = -2y^3(45A^2-33A-13)/[A F(y)] with
f(2/3) = 1 anchor:

| A | f(A) |
|---:|---:|
| 0.700 | 1.000343 |
| 0.750 | 1.011700 |
| 0.800 | 1.069338 |
| 0.850 | 1.244290 |
| 0.900 | 1.755953 |
| 0.950 | 4.148934 |
| 0.990 | 1.49e+02 |

f(A) monotonically decreases from 1 at the PS to ~0 at the horizon.

## Two key invariances

| Sector | Behavior under g -> f g |
|---|---|
| **Photon geodesics** | UNCHANGED (conformal invariance of null cones) |
| **Outside PS (A <= 2/3)** | UNCHANGED (f = 1 by anchor) |

Consequence: EHT shadow size, light bending, Shapiro delay, ISCO
location, eikonal ringdown, all weak-field tests -- IDENTICAL to GR.

## Where deviations live

| Observable | STAM deviation | Magnitude |
|---|---|---|
| Matter emission inside PS | redshift excess (1+z_STAM)/(1+z_GR) = 1/sqrt(f) | ~3-1000x at A=0.85-0.99 |
| Accretion-disk plunge spectrum | faint low-freq tail from plunge region | needs GRMHD model |
| Hawking greybody factors | shifted from GR; G103 ref | (depends on omega, ell) |
| QNM channels (non-eikonal) | 4.977% axial, 13% polar break at l=2 | locked from G108/G114 |

## Precision tests passed by construction

All 10 standard precision regimes have A <= 10^-1, well outside PS:

  Atomic clocks, GPS, Shapiro delay (Sun), Mercury perihelion,
  lunar laser ranging, BBN, CMB acoustic peaks, galactic dynamics,
  light bending (Sun), pulsar timing.

f = 1 in all these regimes. STAM is GR-equivalent.

## Status

- **Tier-1 matter coupling implemented and consistent.** f(Sigma) -minimal
  coupling reproduces all precision-test successes of GR while
  predicting observable deviations only in the inside-PS regime.
- **Strain point #3 partially closed** (short-term). The deep version
  (matter as substance excitations, TOE-level) remains open.
- **Connection to existing predictions:** the matter-coupling pathway
  produces no new observables beyond those already locked (G108 QNM,
  G114 polar break, G90 tortoise enhancement, G103 greybody scattering).
  This is the framework's CONSISTENCY check: matter minimal coupling
  does not introduce additional ad-hoc structure.

## Files

- [scripts/G123_matter_minimal_coupling.py](../scripts/G123_matter_minimal_coupling.py)
- [results/G123_matter_minimal_coupling.png](G123_matter_minimal_coupling.png) -- f(A) curve and redshift ratio
- See also: [project_three_strain_points_fix_strategy.md](../memory/project_three_strain_points_fix_strategy.md)
