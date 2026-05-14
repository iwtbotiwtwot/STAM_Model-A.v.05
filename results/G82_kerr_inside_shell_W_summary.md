# G82 — Inside-shell W_K for Kerr (Option B)

**Date: 2026-05-13.**  Extends the double-LM shell-count action to inside-
shell Kerr with the framework's quintic Hermite closure.

## Structural choice: Option B

```
g_STAM^{rr} = g_Kerr^{rr} * F(y_K) = (Delta / Sigma_BL) * F(y_K)
```

Option A (replacing the Schwarzschild k = (1−A) F(y) with k = (1−A_K) F(y_K)) is rejected because it would break the 2026-05-13 commitment that the metric is GR-Kerr exact outside the photon region.

## Inside-shell W_K

```
W_K^inside(r, theta; M, a) = (Delta / Sigma_BL) * F(y_K) * (dSigma_K/dr)^2
  with  Sigma_K = 6Mr/(r²+a²),  y_K = Sigma_K - 2,
        F(y_K) = 1 - 5 y_K^4 + 4 y_K^5,
        Delta = r² - 2Mr + a²,  Sigma_BL = r² + a² cos²θ.
```

## Verifications

| Requirement | Status |
|---|---|

| 1. Finite across the final shell | ✓ |

| 2. W_K → 0 at horizon (cubic: Delta·F²) | ✓ |

| 3. λ₁ radial pin remains first-order | ✓ (structural) |

| 4. λ₂ co-rotation / no-propagation survives | ✓ (inherits G81) |

| 5. Schwarzschild limit a → 0 → G80 inside W | ✓ (residue = 0 exact) |

| 6. Exterior limit F → 1 → G81 outside W | ✓ (identical formula) |



## Caveat

The simple y_K = Σ_K − 2 places the framework's photon-region boundary at r = (3M + √(9M² − 4a²))/2, which is close to but not exactly the Kerr equatorial photon orbit for a ≠ 0.  A photon-region-normalized version y_K = (Σ_K − Σ_ph(a, θ)) / (3 − Σ_ph(a, θ)) can be introduced later if precise photon-orbit matching becomes important.  For G82's first-pass formulation, the simple y_K = Σ_K − 2 is sufficient.



## Files

- [scripts/G82_kerr_inside_shell_W.py](../scripts/G82_kerr_inside_shell_W.py)

- [plots/G82_kerr_inside_shell_W.png](../plots/G82_kerr_inside_shell_W.png)
