# G58 - Ledger configuration volume -> strong-field metric

**Date: 2026-05-12 (continuation).** Tests the ledger-configuration-volume derivation of the final-shell closure profile F(Sigma) and explores symphony-refinement (non-uniform alpha_i) variants.

## Setup

Under substance ontology + presentism, the ledger at final-shell depth y = Sigma - 2 is a continuous configuration of A-resolution. With D + 2 channels (D spatial + 2 horizon-pair) and uniform per-entry measure, the configuration-volume measure on the spatial-vs-horizon-pair split is:

```
g(y)  ~  y^(alpha_S - 1) * (1 - y)^(alpha_H - 1)
    alpha_S = sum of alpha_i over spatial channels   = D    (uniform)
    alpha_H = sum over horizon-pair channels         = 2    (uniform, two-face)
```

For D = 3 with uniform measure: Beta(3, 2) = 12 y^2 (1-y).

Survival: F(y) = 1 - I_y(3, 2) = **1 - 4 y^3 + 3 y^4**  (the quartic Hermite).

## Verified equivalence (C1)

Beta(3, 2) Dirichlet density and the quartic Hermite F match to numerical precision across the final shell. Ledger-configuration-volume gives the same p(y) and F(y) as the simplex argument.

## Candidates tested

| Candidate | a_S | a_H | PS smoothness | ord_{A=1} k | Verdict |
|---|---:|---:|---|---:|---|
| C1: uniform D=3, two-face | 3 | 2 | C^2 | 3 | **matches all constraints** |
| C2: pair-weighted (each x2) | 6 | 4 | C^5 | 5 | ord != D (gives 5) |
| C3: higher-alpha spatial | 4 | 2 | C^3 | 3 | **matches all constraints** |
| C4: higher-alpha horizon | 3 | 3 | C^2 | 4 | ord != D (gives 4) |
| C5: softer (cubic smoothstep) | 2 | 2 | C^1 | 3 | fails C^2 at PS |
| C6: thicker horizon | 3 | 4 | C^2 | 5 | ord != D (gives 5) |

## Joint constraint analysis

Three independent structural constraints:

- **SU shell-count**: ord_{A=1} k(A) = D

- **C^2 at PS**: alpha_S > 2 (so F''(0) = 0)

- **Two-face refinement**: alpha_H = 2 (geometrically, two horizon faces; independent of D)



From the third: ord_{A=1} k = alpha_H + 1 = 3 **automatically**, in any spatial dimension D.



From the first: ord = D required.



**Joint compatibility forces D = 3.** This is the structural derivation: D = 3 is the unique spatial dimension where SU shell-count and the two-face horizon refinement co-determine the same horizon closure order.



Combined with C^2 at PS: among uniform-measure candidates with integer alpha, only **Beta(3, 2)** satisfies all three constraints. The strong-field metric is structurally unique under these commitments.

## QNM consequence (spinless eikonal): robust at the photon sphere

For **all candidates** with normalized F (so F(0) = 1):

```
k(A = 2/3) = (1 - 2/3) * F(0) = (1/3) * 1 = 1/3   (exact Schwarzschild)
```
So **tau_STAM / tau_GR = 1.000 (spinless eikonal) for every ledger-measure candidate.** The choice of alpha doesn't shift the spinless ringdown prediction.



Confirms G57: spinless STAM ringdown = exact GR is robust to any sensible ledger-measure choice, not specific to the quartic Hermite.

## Where candidates differ: inside-PS metric shape

All differences live in the survival ratio F(y) = k_STAM / k_GR for A > 2/3. Sample at representative depths:


| A | C1 Beta(3,2) | C2 Beta(6,4) | C5 Beta(2,2) |
|---:|---:|---:|---:|
| 0.70 | 0.9963 | 0.9999 | 0.9720 |
| 0.80 | 0.8208 | 0.9006 | 0.6480 |
| 0.90 | 0.3483 | 0.2703 | 0.2160 |
| 0.99 | 0.0052 | 0.0001 | 0.0026 |

Distinguishing observables (none observationally settled yet):

- Higher overtones n >= 1: probe geometry just inside r_c

- Late-inspiral chirp: A approaching 2/3 from outside

- LISA EMRI ringdowns: sample intermediate A precisely

- Sub-leading WKB or full Regge-Wheeler on n=0

## What this does for the framework

1. **Ledger-config-volume = simplex measure under existing commitments.** Substance ontology + presentism naturally produce continuous Dirichlet measure; uniform per-entry weighting is the natural expression of presentism (no internal hierarchy of entry-types).

2. **Substance + presentism are now load-bearing for the strong-field metric**, not just for the quantum interpretation. Tighter framework integration.

3. **D = 3 is forced by SU shell-count + two-face refinement + Beta horizon exponent.** 'Why three spatial dimensions' shifts from cosmological/anthropic to framework-internal.

4. **Spinless ringdown = exact GR is a robust prediction**, not specific to the quartic. Any normalized ledger-measure candidate gives F(0) = 1 -> k(2/3) = 1/3 = Schwarzschild.

## Open follow-ups

- Whether **uniform per-entry measure (alpha_i = 1)** is itself derivable from presentism, or is a fourth independent commitment. Best reading: under presentism, the ledger has no internal hierarchy of entry-types, so uniform is natural -- but it could be stated explicitly.

- A **symphony-refinement** observation prediction for late-inspiral chirp or EMRI ringdown shape, distinguishing C1 from non-uniform-alpha candidates.

- Confirm the **SU shell-count argument is independent** of the horizon-pair argument (otherwise the D = 3 forcing dissolves into a tautology).

## Files

- [scripts/G58_ledger_configuration_volume.py](../scripts/G58_ledger_configuration_volume.py)

- [plots/G58_p_and_F_comparison.png](../plots/G58_p_and_F_comparison.png)

- [plots/G58_k_inside_PS.png](../plots/G58_k_inside_PS.png)

- [plots/G58_survival_ratio.png](../plots/G58_survival_ratio.png)
