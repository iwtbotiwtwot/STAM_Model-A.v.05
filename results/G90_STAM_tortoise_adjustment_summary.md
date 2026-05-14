# G90 — STAM tortoise-coordinate adjustment near horizon

**Date: 2026-05-13.**  Numerical computation of r*_GR(ε), r*_STAM(ε), Δr*(ε), and ratio r*_STAM/r*_GR at cutoffs ε = (r-2M)/(2M) = 10⁻¹, ..., 10⁻⁶, and verification of the local STAM traversal-adjustment formula 1/√F(y).

## Local adjustment formula (confirmed)

```
(dr*/dr)_STAM / (dr*/dr)_GR  =  1 / sqrt(F(y))
```

with F(y) = 1 − 5y⁴ + 4y⁵, y = 3A − 2.  Inside the photon sphere (y > 0), this adjustment factor is greater than 1; outside (y ≤ 0), F = 1 and the adjustment is unity (STAM ≡ GR).



## Cutoff table  (M = 1, r*(3M) = 0)

| ε | r/M | r*_GR | r*_STAM | Δr* | STAM/GR |
|---:|---:|---:|---:|---:|---:|

| 1e-01 | 2.200000 | -4.018876e+00 | -4.542922e+00 | -5.240466e-01 | 1.130396e+00 |

| 5e-02 | 2.100000 | -5.505170e+00 | -7.504801e+00 | -1.999631e+00 | 1.363228e+00 |

| 1e-02 | 2.020000 | -8.804046e+00 | -2.616053e+01 | -1.735648e+01 | 2.971421e+00 |

| 5e-03 | 2.010000 | -1.020034e+01 | -4.798445e+01 | -3.778411e+01 | 4.704201e+00 |

| 1e-03 | 2.002000 | -1.342722e+01 | -2.183449e+02 | -2.049176e+02 | 1.626136e+01 |

| 1e-04 | 2.000200 | -1.803419e+01 | -2.118141e+03 | -2.100106e+03 | 1.174514e+02 |

| 1e-05 | 2.000020 | -2.263954e+01 | -2.109423e+04 | -2.107159e+04 | 9.317432e+02 |

| 1e-06 | 2.000002 | -2.724472e+01 | -2.108333e+05 | -2.108061e+05 | 7.738501e+03 |



## Asymptotic scaling near horizon

- **GR**: r*_GR(ε) ~ 2M ln(ε) — logarithmic divergence.

- **STAM**: r*_STAM(ε) ~ -C/ε — power-law divergence (much faster).

- **Ratio**: r*_STAM/r*_GR → ∞ as ε → 0.

- **Difference**: |Δr*| ~ C/ε with C ≈ 0.105 (leading order from F ~ 10(1-y)²).



## Physical reading

Near the horizon, the STAM substance density approaches saturation and the quintic Hermite closure F(y) drives k(r) to zero **cubically** (matching the order-D zero of the SU shell-count commitment). This means signals traversing the final shell accumulate **exponentially more tortoise-coordinate path** than they would in GR — a structural near-horizon distance prediction.



## Observable consequences

- Hawking emission spectrum reads off the modified near-horizon r* structure.

- Ringdown late-time power-law tails sample the framework's deeper near-horizon region.

- Black-hole imaging / shadow may sense the modified radial gradient near r = 2M.



## Files

- [scripts/G90_STAM_tortoise_adjustment.py](../scripts/G90_STAM_tortoise_adjustment.py)

- [plots/G90_STAM_tortoise_adjustment.png](../plots/G90_STAM_tortoise_adjustment.png)
