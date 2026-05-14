# STAM Model-A F1: Penrose-Hawking Singularity Theorem Evasion

## The attack

Penrose's 1965 singularity theorem states: in any spacetime satisfying
1. The null energy condition: `T_μν k^μ k^ν ≥ 0` for all null `k^μ`,
2. A non-compact Cauchy surface,
3. **The existence of a closed (compact, no-boundary) *strictly* trapped surface**, meaning a closed 2-surface where both null expansions θ_+ and θ_- are strictly negative,
the spacetime is null-geodesically incomplete — i.e. it contains a singularity.

An opponent will say: bold STAM denies that black holes have singularities. But Penrose proved a singularity must form when matter collapses inside a horizon. Therefore bold STAM contradicts a proven theorem of mathematical physics.

## Bold STAM's defense

Bold STAM does NOT deny the theorem. **The theorem's premise (3) is unmet** in bold STAM. The theorem applies to spacetimes that contain a strict closed trapped surface. Bold-STAM spacetimes do not.

More precisely: bold STAM has a fundamentally different geometric class than a GR Schwarzschild black hole. In Schwarzschild, the interior region `r < Rs` (i.e. `A > 1`) contains spheres on which both null expansions are strictly negative — these are the trapped surfaces that the theorem requires. In bold STAM, **the manifold does not extend to A > 1 at all** (the universe ends at the boundary). There is therefore no interior region in which strict trapped surfaces could form.

## The math: outgoing null expansion θ_+

For a sphere of areal radius `r = Rs/A` in a spherically symmetric metric `ds² = -h(A)c²dt² + dr²/k(A) + r²dΩ²`, the outgoing null expansion (in Schwarzschild-time parameterization) is:
```text
θ_+(r) = (2/r) × sqrt(h(A) k(A))
```

**Schwarzschild** (h = k = 1-A):
```text
θ_+ = (2/r)(1-A)
    > 0  for A < 1   (untrapped, exterior)
    = 0  at A = 1    (marginally trapped, horizon)
    < 0  for A > 1   (strictly trapped, interior)
```

**Bold STAM** (h = 1-A, k = (1-A)(1-A²)²):
```text
θ_+ = (2/r)(1-A)(1-A²)
    > 0  for A < 1   (untrapped)
    = 0  at A = 1    (marginally trapped, with double zero — even more degenerate)
    [no manifold for A > 1]
```

|      A |   r_over_Rs |   theta_plus_Schwarzschild_per_Rs |   theta_plus_bold_STAM_per_Rs |   ratio_bold_over_Schw |
|-------:|------------:|----------------------------------:|------------------------------:|-----------------------:|
| 0.1    |    10       |                        0.18       |                    0.1782     |             0.99       |
| 0.3    |     3.33333 |                        0.42       |                    0.3822     |             0.91       |
| 0.5    |     2       |                        0.5        |                    0.375      |             0.75       |
| 0.7    |     1.42857 |                        0.42       |                    0.2142     |             0.51       |
| 0.9    |     1.11111 |                        0.18       |                    0.0342     |             0.19       |
| 0.99   |     1.0101  |                        0.0198     |                    0.00039402 |             0.0199     |
| 0.999  |     1.001   |                        0.001998   |                    3.994e-06  |             0.001999   |
| 0.9999 |     1.0001  |                        0.00019998 |                    3.9994e-08 |             0.00019999 |


**Verdict on premise (3):** in bold STAM, no closed surface has both null expansions strictly negative. The strict-trapped condition is never satisfied. The theorem's premise is unmet, so its conclusion (incomplete null geodesic) is not forced.

## Geodesic completeness — automatic consequence

A direct check: compute the affine parameter `λ` along an ingoing radial null geodesic from a starting position to the horizon. In Schwarzschild this is finite (geodesic incomplete; reaches r=Rs in bounded affine parameter, then continues into singularity at r=0). In bold STAM the integral diverges logarithmically at A=1 (geodesic asymptotes to boundary, never reaches it in finite affine parameter).

|   A_start |   r_start_over_Rs |   lambda_to_horizon_Schwarzschild_Rs_over_c |   lambda_to_horizon_bold_STAM_Rs_over_c |   lambda_ratio_bold_over_Schw |
|----------:|------------------:|--------------------------------------------:|----------------------------------------:|------------------------------:|
|    0.5    |           2       |                                  1          |                                36.8306  |                       36.8306 |
|    0.9    |           1.11111 |                                  0.111111   |                                10.7918  |                       97.1263 |
|    0.99   |           1.0101  |                                  0.010101   |                                 4.81537 |                      476.722  |
|    0.999  |           1.001   |                                  0.001001   |                                 3.45771 |                     3454.25   |
|    0.9999 |           1.0001  |                                  0.00010001 |                                 2.30273 |                    23025      |


Bold-STAM null geodesics are complete in both directions: they approach but never reach the boundary. This is consistent with — and a consequence of — the absence of trapped surfaces in the strict sense.

## Where bold STAM and Schwarzschild differ as geometric objects

- **Schwarzschild manifold**: extends to A → ∞ (r → 0). Trapped region exists for A > 1. Singularity at A = ∞ (r = 0). Penrose's theorem applies.
- **Bold-STAM manifold**: A < 1 only. No interior. No trapped region. No singularity. Penrose's theorem premise unmet — theorem does not apply.

Bold STAM is not denying Penrose. Bold STAM is constructing a *different geometric object* that the theorem does not address. The theorem says 'IF a strict trapped surface exists, THEN incomplete geodesics exist.' Bold STAM says 'no such surface exists in our manifold.' Both statements are simultaneously true.

## What this gives bold STAM

- A clean answer to 'how do you avoid singularities Penrose proved must form?': by having a manifold class outside the theorem's domain. The theorem applies to spacetimes with strictly trapped surfaces; bold-STAM spacetimes don't have them.
- A *positive* statement of the bubble picture: gravitational collapse in bold STAM asymptotically approaches a marginally-trapped boundary (θ_+ = 0) without ever forming a strictly trapped region (θ_+ < 0). The bubble surface IS the marginally trapped boundary; matter accumulates on it, never enters a (nonexistent) interior.
- Geodesic completeness as a derived consequence, not a postulate.

## What this fatality establishes for the framework

**Fatality 1: SURVIVED.** Bold STAM is consistent with Penrose's 1965 theorem. The theorem applies to spacetimes containing strictly trapped surfaces; bold-STAM spacetimes don't contain them. There is no contradiction with the theorem; bold STAM is in a different geometric class.

Caveat (honest): we have not yet shown that *gravitational collapse dynamics* in bold STAM smoothly produce a manifold of this class without ever transiently forming a strictly trapped region. That requires dynamical-A propagation work (open). The static-state argument here shows the equilibrium configuration evades the theorem; the collapse-process argument is a separate piece of work.

## Generated plots

- `plots/F1_theta_comparison.png`
- `plots/F1_affine_parameter_divergence.png`
- `plots/F1_trapped_region_topology.png`
