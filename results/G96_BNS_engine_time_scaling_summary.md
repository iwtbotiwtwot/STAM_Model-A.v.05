# G96 — BNS engine-time scaling test
## Interpretation
The observed GW-to-EM delay is treated as **local engine time near the merger**, not as a long-distance propagation-speed difference. GW and light propagate essentially together after release.
## Locked derivation
```text
A0 = 1/(12π)
A0 = x(1-x)/(1+x), where x = Rs/(2r_threshold)
r_threshold/Rs = 1/(2x)
tau_engine = (5/8)(r_threshold/Rs)^4 (Rs/c)
tau_engine = slope × (M_total/M_sun)
```
- A0 = `0.026525823849`
- x roots = `0.028057277031`, `0.945416899120`
- chosen root = `0.028057277031`
- r_threshold/Rs = `17.820689`
- slope = `0.620970` s/M_sun
## Event comparison
| Event | M_total (M_sun) | Observed delay (s) | STAM predicted (s) | Obs - Pred (s) | Abs % error |
|---|---:|---:|---:|---:|---:|
| GW170817 / GRB 170817A | 2.7000 | 1.7400 | 1.6766 | 0.0634 | 3.643% |

## Falsification handle
Future BNS + EM counterpart events should fall approximately on the linear mass-scaling relation if the STAM doughnut / substance-engine-time mechanism is correct. If events with reliable EM launch-time interpretation do not follow this scaling, this mechanism is wrong or incomplete.
## Files
- CSV: `/mnt/data/G96_run/results/G96_BNS_engine_time_scaling_events.csv`
- Plot: `/mnt/data/G96_run/plots/G96_BNS_engine_time_scaling.png`
