# G17: PPN Parameter Calculation for Model-A
## Setup
Compute the parameterized post-Newtonian (PPN) parameters γ and β for Model-A's strong-field metric and compare to GR Schwarzschild and to current solar-system observational bounds. Non-distance test — probes weak-field deviations of Model-A from GR using local spacetime structure.

Model-A metric:
```text
g_tt = -(1 − A) c²
g_rr = 1 / [(1 − A)(1 − A²)²]
g_θθ = r²,  g_φφ = r² sin²θ
A = 2GM/(c²r)
```

PPN formalism requires isotropic coordinates. The coordinate transformation r → r̃ is solved symbolically (sympy) to high enough order to extract γ and β.
## Result
| Theory | γ (PPN) | β (PPN) | α₂ (g_ij U² coeff) |
|---|:---:|:---:|:---:|
| GR | 1 | 1 | 3/2 |
| Model-A | 1 | 1 | 11/2 |

## Reading
**γ (Cassini-bounded at 10⁻⁵):**
- GR: γ − 1 = 0
- Model-A: γ − 1 = 0

Model-A predicts γ identical to GR. Cassini bound is automatically satisfied. Light deflection and Shapiro delay match GR exactly at the level current measurements can distinguish.

**β (Mercury-bounded at 10⁻⁴):**
- GR: β − 1 = 0
- Model-A: β − 1 = 0

Model-A predicts β identical to GR. Mercury perihelion shift, lunar laser ranging, and other β-dependent solar-system tests match GR exactly.

**Conclusion on first/second-order PPN:** Model-A passes every current solar-system PPN constraint automatically. The framework's modified k(A) does not show up in any PPN-bounded observable.
## Where Model-A first differs from GR
Model-A's deviation from GR appears in the third-order coefficient of the spatial isotropic metric (α₂ in `g_ij = (1 + 2γu + α₂ u² + α₃ u³ + ...) δ_ij`):

- GR α₂ = 3/2
- Model-A α₂ = 11/2
- Difference: 4

This is **third-order PPN territory**, where solar-system constraints are essentially absent (no current measurement reaches this precision). Model-A is observationally indistinguishable from GR at all current solar-system tests, but the deviation is real and would show up in any future measurement that reaches third-order PPN sensitivity.
## Coordinate-dependent structure (worth noting)
In **Schwarzschild coordinates**, the second-order coefficient of `g_rr` (in U_Sch = GM/(c²r)) is:
- GR: 4
- Model-A: 12
- Ratio: 3

**Model-A's Schwarzschild-coordinate g_rr has exactly 3× the second-order coefficient of GR.** This is coordinate-dependent — in isotropic coordinates the ratio shifts away from exactly 3 — but the appearance of an integer factor of 3 in the natural Schwarzschild-coordinate expansion is suggestive in light of the framework's other 'thirds' structure (ISCO/photon-sphere/horizon at A = 1/3, 2/3, 1).

Whether this is the same '3' appearing structurally in different places is the open question we've been circling. This script doesn't answer it — it surfaces another instance and notes it.
## Implications
- **Model-A is solar-system-safe.** No current PPN test constrains the framework. Every test that has bounded γ or β to 10⁻⁵–10⁻⁴ precision is automatically passed.
- **Model-A is solar-system-indistinguishable from GR.** The framework's modified k(A) does not produce any second-order PPN deviation. The first signature appears at third-order, currently unconstrained.
- **The falsifiable corners are not in the solar system.** As the strong-field commitment memory already states, the genuine tests are LIGO QNM ringdown (G1's τ_MA/τ_GR = 1.80), EHT shadow structure beyond size, and any near-horizon Shapiro-class measurement at A → 1.
- **The 3× ratio in Schwarzschild g_rr** is coordinate-dependent but worth tracking. If a coordinate-invariant statement of this structure exists, it would be another '3' to add to the framework's list — alongside the strong-field thirds and the (4π × 3) decomposition of A_0. Whether they are the same 3 is open.
