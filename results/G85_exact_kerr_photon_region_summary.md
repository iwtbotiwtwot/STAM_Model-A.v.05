# G85 -- exact Kerr spheroidal photon-region surface

Exact derivation from spherical-null-geodesic conditions R(r) = 0, dR/dr = 0, replacing G84's sin^2(theta) interpolation between r_ph_prograde(a) and r_polar(a).

## Formulas
```
lambda(r_p) = -(r_p^3 - 3 r_p^2 + a^2 r_p + a^2) / [a (r_p - 1)]
eta(r_p)    = r_p^2 [ 4 a^2 Delta - (r_p^2 - 3 r_p + 2 a^2)^2 ]
              / [ a^2 (r_p - 1)^2 ]
Theta(theta)/E^2 = eta - cos^2(theta) [ lambda^2/sin^2(theta) - a^2 ]
```

r_pr_exact(theta; a) is the unique r_p in [r_ph_prograde, r_polar] such that Theta(theta; lambda(r_p), eta(r_p)) = 0.

## Comparison with G84

- max abs |r_pr_exact - r_pr_G84| = 5.3937e-01 M
- max rel |r_pr_exact - r_pr_G84| / r_pr_exact = 39.3454%
- max abs |Sigma_ph_exact - Sigma_ph_G84| = 4.0314e-01
- max rel = 14.7667%

## Verdict
G84 ansatz deviates from exact by up to 39.35%.  Significant; the exact G85 r_pr(theta; a) should replace the sin^2 interpolation in the W_K formula for any precision work.

## Files
- [scripts/G85_exact_kerr_photon_region.py](../scripts/G85_exact_kerr_photon_region.py)
- [plots/G85_exact_kerr_photon_region.png](../plots/G85_exact_kerr_photon_region.png)
