# G108 - Rigorous axial perturbation equation from S[Sigma, g, lambda1, lambda2]

## Result
- omega_GR (calibration): (0.3706322594379458-0.08929019551728928j)  vs Leaver (0.373672-0.088962j)
- calibration error: 0.7960%  (threshold 1.0%)
- omega_exact (analytic correction): (0.38837731849782914-0.08257234484597935j)  shift 4.9770%
- omega_canon_FD (G92 correction):   (0.3624655017657346-0.09014496884234341j)  shift 2.1539%
- |exact - canon_FD| = 2.8231%

## Verdict
V_exact and V_canonical DIFFER materially.  Either
           the analytic computation has a sign/factor bug or the
           G92 finite-difference correction was missing structure.
           Investigate before locking either as the framework
           prediction.

## Action correction comparison
- RMS analytic vs FD difference on (sqrt f)''/sqrt f: 4.6857e-03
- max abs analytic vs FD: 6.0961e-02

## Files
- [scripts/G108_axial_from_constrained_sigma_action.py](../scripts/G108_axial_from_constrained_sigma_action.py)
- [plots/G108_axial_from_constrained_sigma_action.png](../plots/G108_axial_from_constrained_sigma_action.png)
