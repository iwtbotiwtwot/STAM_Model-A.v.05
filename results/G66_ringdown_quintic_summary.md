# G66 - Ringdown Re-run with Quintic Hermite F(y)

**Date: 2026-05-13.** Updates [G57](G57_ligo_ringdown_under_current_metric_summary.md) which used the quartic Hermite, to use the **quintic Hermite** from the G58 ledger-channel commitment.

## What changed since G57

- **G57 used quartic Hermite**: `F = 1 - 4 y^3 + 3 y^4` (C^2 at PS, C^1 at horizon).

- **G58 ledger-channel commitment** gave quintic Hermite: `F = 1 - 5 y^4 + 4 y^5` (C^3 at PS, C^1 at horizon).

- Both have F(0) = 1, so at the photon sphere k = (1-2/3)*1 = 1/3 = exact Schwarzschild value. Eikonal Lyapunov identical.

## Result: ringdown prediction unchanged

```
tau_STAM / tau_GR = 1.000000 exactly (eikonal, spinless)
```

Verified for both quartic and quintic forms. The G58 structural commitment shift did not affect the LIGO observable.

## F derivatives at PS (y=0)

| Derivative | Quartic (G57) | Quintic (G58) | GR (k=1-A) |

|---|---:|---:|---:|

| F(0) | 1 | 1 | 1 |

| F'(0) | 0 | 0 | 0 |

| F''(0) | 0 | 0 | 0 |

| F'''(0) | -24 | 0 | 0 |

| F''''(0) | 72 | -120 | 0 |



Quintic matches GR through F'''(0); quartic only through F''(0). Both differ from GR at F''''(0).

## Practical sensitivity (LIGO)

- **Eikonal (leading)**: both forms give exact GR (no change from G57).

- **First-order WKB beyond eikonal**: both match GR through F''(0); no change.

- **Second-order WKB**: tiny difference (sub-percent) between quartic and quintic, both still close to GR. Below LIGO precision for n=0 fundamental.

- **Higher overtones (n >= 1)**: sample geometry deeper; can in principle distinguish quartic from quintic at very high SNR.

## Bottom line

The structural commitment shifted (quartic -> quintic) between G57 and G58, but **the spinless ringdown prediction is unchanged**: tau_STAM/tau_GR = 1 exactly in eikonal, both forms. G57's main result holds under the current commitment. Updated for the record.

## Files

- [scripts/G66_ringdown_quintic.py](../scripts/G66_ringdown_quintic.py)

- [plots/G66_F_comparison.png](../plots/G66_F_comparison.png)

- [plots/G66_k_quartic_vs_quintic.png](../plots/G66_k_quartic_vs_quintic.png)

- [plots/G66_ringdown_quintic.png](../plots/G66_ringdown_quintic.png)
