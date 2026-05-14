# G112 - Higher-mode QNM via residual subtraction

Replacement for G110 (two-mode fit failed) and G111 (matrix pencil drifted).
Standard Prony-type residual subtraction: fit n=0 at late times,
subtract template, fit residual at early times for n=1.

## Results
| ell | n | calib | omega_GR | Leaver | omega_exact | shift | status |
|---:|---:|---:|---|---|---|---:|---|
| 2 | 0 | 0.049% | (0.3734906273003711-0.0890064116773112j) | (0.373672-0.088962j) | (0.39215572929660497-0.0831402294272712j) | 5.0958% | OK |
| 2 | 1 | 31.074% | (0.3564542616232808-0.13695750000000004j) | (0.346711-0.273915j) | (0.4084465298737619-0.13695750000000004j) | 13.6155% | POOR_CALIB |
| 3 | 0 | 0.005% | (0.5994358850001957-0.0927326841469219j) | (0.599443-0.092703j) | (0.6115113940404647-0.08866300821337274j) | 2.1008% | OK |
| 3 | 1 | 21.992% | (0.6041768460522794-0.14064900000000002j) | (0.582644-0.281298j) | (0.6110549470086193-0.14064900000000002j) | 1.1088% | POOR_CALIB |
| 4 | 0 | 0.006% | (0.8091609808885415-0.0942099761364803j) | (0.809178-0.094164j) | (0.8179156026950651-0.09069363252304617j) | 1.1581% | OK |
| 4 | 1 | 17.032% | (0.8199302999685578-0.14216700000000002j) | (0.796632-0.284334j) | (0.8260621745326187-0.14216700000000002j) | 0.7369% | POOR_CALIB |

## Files
- [scripts/G112_residual_subtraction_overtones.py](../scripts/G112_residual_subtraction_overtones.py)
- [plots/G112_residual_subtraction_overtones.png](../plots/G112_residual_subtraction_overtones.png)
