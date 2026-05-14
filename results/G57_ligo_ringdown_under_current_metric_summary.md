# G57 — LIGO Ringdown Under Current k(A) Commitment

**Date: 2026-05-12 (continuation).** Test re-run under the strong-field metric committed in [strong_field_metric_commitment_2026_05_12.md](strong_field_metric_commitment_2026_05_12.md) and described in the current [README](../README.md). Replaces [G1](../scripts/G1_qnm_ringdown.py) for the spinless QNM ringdown prediction.

## What changed

G1 was written when the framework's strong-field metric was `k(A) = (1-A)(1-A^2)^2`. Under the current commitment:

```text
k(A) = (1 - A)                          for A <= 2/3       (outside PS, exact GR)
k(A) = 27(1-A)^3 (9A^2 - 10A + 3)       for 2/3 < A < 1    (final-shell quartic Hermite)
```

**At the photon sphere (A = 2/3):**

| Metric                | k(A_PS)              |
|---|---:|
| Schwarzschild         | 1/3 = 0.333333 |
| STAM obsolete (G1)    | 25/243 = 0.102881 |
| **STAM current**      | **1/3 = 0.333333**  *(= GR exactly)* |

Because both `h` and `k` are now identical to Schwarzschild at the photon sphere, the eikonal Lyapunov exponent `lambda^2 = h*k / (3 M^2)` is also identical to Schwarzschild's:

```text
lambda_Schwarzschild  = 0.192450 / M
lambda_STAM_obsolete  = 0.106917 / M   (was 5/9 smaller -> hence the historic 1.80)
lambda_STAM_current   = 0.192450 / M   (= Schwarzschild)
```

## Spinless eikonal prediction

| Metric             | tau_STAM / tau_GR (spinless, eikonal) |
|---|---:|
| G1 (obsolete)       | **1.8000**   (the historic 9/5 = 1.80) |
| Current             | **1.0000**   (exact GR) |

The framework's earlier signature prediction (1.80x longer ringdown) is **retired** under the current k(A). Spinless STAM produces the Schwarzschild eikonal QNM exactly because the metric is exactly GR at and outside the photon sphere.

More strongly: the quartic Hermite F(Sigma) was constructed with `F(2) = 1, F'(2) = 0, F''(2) = 0` (C^2 at PS), so `k`, `dk/dA`, and `d^2k/dA^2` all agree with the outside-PS GR form at A = 2/3. First-order WKB beyond eikonal is also GR-identical at the photon sphere. Sub-leading departures require higher-order WKB corrections or the inside-PS region.

## Predictions for representative LIGO events (spinless eikonal)

```text
Event                                 M (M_sun)      f (Hz)   tau_GR (ms)   tau_old (ms)   tau_new (ms)
----------------------------------------------------------------------------------------------------
GW150914 remnant                             62       200.5         3.174          5.714          3.174
GW170729 remnant                             80       155.4         4.096          7.373          4.096
GW190521 remnant                        1.4e+02        87.6          7.27          13.09           7.27
Stellar-mass BBH (typical)                   30       414.5         1.536          2.765          1.536
Intermediate-mass (~1000 M_sun)           1e+03        12.4          51.2          92.16           51.2
Sgr A*                                  4.3e+06         0.0     2.202e+05      3.963e+05      2.202e+05
M87*                                    6.5e+09         0.0     3.328e+08       5.99e+08      3.328e+08
```
Frequencies are identical across all three metrics because they depend only on `h`, which is unchanged. Damping time tau is what the metric modification of g_rr would have shifted; under the current k, it doesn't.

## Where the current k still differs from GR (inside the photon sphere)

Inside the photon sphere (A > 2/3), the current k drops sharply and reaches zero with order D = 3 at the horizon. Sampling:


| A | k_GR | k_STAM(current) | ratio |
|---:|---:|---:|---:|
| 0.70 | 0.3000 | 2.9889e-01 | 0.9963 |
| 0.75 | 0.2500 | 2.3730e-01 | 0.9492 |
| 0.80 | 0.2000 | 1.6416e-01 | 0.8208 |
| 0.85 | 0.1500 | 9.1353e-02 | 0.6090 |
| 0.90 | 0.1000 | 3.4830e-02 | 0.3483 |
| 0.95 | 0.0500 | 5.4759e-03 | 0.1095 |
| 0.99 | 0.0100 | 5.1864e-05 | 0.0052 |

Eikonal QNM is photon-sphere-localized and doesn't see this. What does see it:

- **Higher overtones (n >= 1)**: increasingly sensitive to the geometry just inside r_c. High-SNR overtone extraction in LIGO O5+ could probe this directly.

- **Higher-order WKB / full Regge-Wheeler**: integrates over the radial profile, picking up the inside-PS region. Needed for sharp sub-percent predictions on n=0 itself.

- **Late-inspiral chirp shape**: probes A close to 2/3 from the merger side.

- **LISA EMRI ringdowns**: probe intermediate A with high precision.

## Spinning case (LIGO observed remnants)

LIGO BBH remnants are spinning Kerr-like, typically a ~ 0.5-0.7. Per [HANDOFF_2026_05_12_strong_field.md](../memory/HANDOFF_2026_05_12_strong_field.md), the static-bubble + rotating-matter ontology gives:

```text
tau_STAM / tau_GR_Kerr  ~  0.87   at a ~ 0.67   (13% deficit)
```

This is within current LIGO ringdown precision (+/- 20-30%). It is a **directional prediction**, not a sharp number -- the explicit STAM-Kerr metric has not been derived. Two readings of the spinning extension (naive A_Kerr_BL vs static-bubble) currently give different predictions. Closing that gap is the next strong-field deliverable.

## Status summary

- **Spinless eikonal**: tau_STAM / tau_GR = 1.0000 (was 1.8000 under G1)

- **First-order WKB beyond eikonal**: also exact GR (C^2 smoothness at PS)

- **Spinning Kerr ringdown**: ~ 13% deficit estimate, contingent on STAM-Kerr extension

- **G1 result (1.80) is preserved** in [G1_qnm_ringdown.py](../scripts/G1_qnm_ringdown.py) as the historical record of the obsolete k(A) form

- **G53 candidate-n exploration** is also superseded; current k is not in the (1-A)(1-A^2)^n family but a piecewise form with explicit outside-PS = GR

## Implications

1. **Spinless LIGO consistency is now automatic.** The framework no longer has a sharp tau-ratio tension to defend against current LIGO ringdown bounds.

2. **The strong-field STAM-vs-GR wedge is pushed inside the photon sphere.** All photon-sphere-localized observables (n=0 ringdown frequency and damping in eikonal, EHT shadow, light bending up to the PS) are now GR-exact in the spinless case. STAM departures live inside r_c.

3. **The framework's main strong-field empirical handle is now Kerr.** Spinless ringdown can no longer distinguish STAM from GR; spinning ringdown (via the static-bubble + rotating-matter reading) is the testable channel. STAM-Kerr extension is gating the framework's near-term LIGO comparison.

## Files

- [scripts/G57_ligo_ringdown_under_current_metric.py](../scripts/G57_ligo_ringdown_under_current_metric.py)

- [plots/G57_k_comparison.png](../plots/G57_k_comparison.png)

- [plots/G57_ringdown_waveforms.png](../plots/G57_ringdown_waveforms.png)

- [plots/G57_inside_ps_departure.png](../plots/G57_inside_ps_departure.png)
