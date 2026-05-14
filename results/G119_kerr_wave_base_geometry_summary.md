# G119 -- Kerr STAM wave base (geometry only)

Replaces the broken Schwarzschild base of G118 with the proper Kerr
radial structure.  Geometry-only -- no Teukolsky physics, no QNM.

## Anchor table
| a | r_+ | r_pr_eq | r_polar | Sigma_ph_eq |
|---:|---:|---:|---:|---:|
| 0.00 | 2.000000 | 3.000000 | 3.000000 | 2.000000 |
| 0.10 | 1.994987 | 2.882194 | 2.995547 | 2.079245 |
| 0.20 | 1.979796 | 2.759193 | 2.982089 | 2.163183 |
| 0.30 | 1.953939 | 2.630026 | 2.959311 | 2.252044 |
| 0.50 | 1.866025 | 2.347296 | 2.883218 | 2.445185 |
| 0.70 | 1.714143 | 2.013334 | 2.757907 | 2.658736 |
| 0.90 | 1.435890 | 1.557855 | 2.559997 | 2.887669 |
| 0.95 | 1.312250 | 1.386281 | 2.492694 | 2.945070 |
| 0.99 | 1.141067 | 1.167642 | 2.430983 | 2.989498 |

## Imaginary / NaN check
- Passed: 9/9 sampled spins.
- Delta, F(y_K), Delta_STAM, and dr*/dr_STAM are real and positive
  throughout the inside-shell at all spins.

## Throat widths
| a | r_+ | r_pr_eq | Δr | r* GR width | r* STAM width | stretch |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 2.0000 | 3.0000 | 1.0000 | 22.6395 | 111131121.3037 | 4908718.9108 |
| 0.10 | 1.9950 | 2.8822 | 0.8872 | 22.3445 | 95408173.6766 | 4269879.2035 |
| 0.20 | 1.9798 | 2.7592 | 0.7794 | 22.1520 | 81900789.4828 | 3697210.5927 |
| 0.30 | 1.9539 | 2.6300 | 0.6761 | 22.0694 | 69971167.5355 | 3170511.0464 |
| 0.50 | 1.8660 | 2.3473 | 0.4813 | 22.3303 | 49141814.6028 | 2200676.8290 |
| 0.70 | 1.7141 | 2.0133 | 0.2992 | 23.6674 | 30459508.1734 | 1286981.4383 |
| 0.90 | 1.4359 | 1.5579 | 0.1220 | 29.7551 | 12159457.6130 | 408650.6826 |
| 0.95 | 1.3122 | 1.3863 | 0.0740 | 36.1287 | 7229957.2690 | 200116.9493 |
| 0.99 | 1.1411 | 1.1676 | 0.0266 | 62.1926 | 2494204.2003 | 40104.5422 |

## Verdict
Wave base is numerically sound at all sampled spins (a in [0, 0.99]).
Delta, F(y_K), Delta_STAM, and dr*/dr_STAM are real, finite, and
appropriately signed throughout the inside-shell.  No Schwarzschild
h(r) < 0 issue (G118 failure mode) recurs.

The wave base is ready for G120 (rigorous Kerr QNM via Detweiler /
Sasaki-Nakamura on this base).

## Files
- [scripts/G119_kerr_wave_base_geometry.py](../scripts/G119_kerr_wave_base_geometry.py)
- [plots/G119_kerr_wave_base_geometry.png](../plots/G119_kerr_wave_base_geometry.png)
