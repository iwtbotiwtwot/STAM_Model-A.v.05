# G29: V_3 modified Friedmann redo — derived H(z) vs chronometers

**Date:** 2026-05-11

## Goal

Compute STAM's derived (not ansatz) H(z) from V_3 modified Friedmann, integrate numerically, and compare to cosmic chronometer data. This is the redo of script 36 (V_1 era) with the G8-selected V_3 = alpha/A + beta/(1-A) with minimum at A_0 by construction. Direct backtrack against G28: does V_3 dynamics fit chronometers where the closed-form ansatz failed?

## V_3 calibration

```
V_3(A) = alpha/A + beta/(1-A)
alpha/beta = (A_0/(1-A_0))^2  =>  V_3 has minimum at A = A_0
A_0 = 1/(12*pi) = 0.026526
beta_tilde = beta/rho_crit = Omega_DE * (1-A_0)^2 = 0.649054
alpha_tilde = beta_tilde * (A_0/(1-A_0))^2 = 0.00048191
V_3(A_0)/rho_crit = 0.6849  (matches Omega_DE target 0.6849)
V_3''(A_0) = 53.0480  (omega_A = 7.2834 in H_0 units)
```

## Scenarios — initial conditions at today, integrated backward

| Scenario | A(today) | Adot(today) | chi^2 | chi^2/N |
|---|---|---|---|---|
| Sc.1 — pure CC | A_0 | 0 | 23.46 | 0.757 |
| Sc.2 — A displaced high | 0.10 | 0 | 32.06 | 1.034 |
| Sc.3 — at min, mild kinetic | A_0 | -0.1 | 24.41 | 0.787 |
| Sc.4 — at min, strong kinetic | A_0 | -0.5 | 151.56 | 4.889 |

## Reference models (recomputed)

| Model | H_0 | chi^2 | chi^2/N |
|---|---|---|---|
| STAM closed-form ansatz | 73.04 | 54.97 | 1.773 |
| LCDM (Om=0.315) | 73.04 | 23.42 | 0.756 |
| LCDM (Om=0.315) | 67.4 Planck | 14.88 | 0.480 |
| EdS | 73.04 | 494.89 | 15.964 |

## Plot

![V_3 modified Friedmann vs chronometers](G29_v3_modified_friedmann.png)

Left: H(z) curves for V_3 scenarios + reference models vs chronometer data.
Right: A(z) trajectory in each scenario (how A evolves over cosmic history).

## Reading

V_3 with A sitting at the minimum (Scenario 1) is degenerate with LCDM by construction — the calibration sets V_3(A_0) = Omega_DE_target * rho_crit, which is exactly LCDM's cosmological constant. Any chi^2 difference between Scenario 1 and LCDM at H_0=73.04 reflects numerical drift in the integration, not physics.

The scenarios with non-trivial A dynamics (Scenarios 2-4) test whether V_3 admits tracking-like behavior that would shift H(z) away from LCDM. **V_3's tight confinement around A_0** (the omega_A = sqrt(V_3'') is large in H_0 units) drives A back to A_0 rapidly under Hubble friction — there's no slow-roll tracking comparable to V_1.

**Backtrack on G28's closed-form result:**

The closed-form ansatz H(z) = H_0(1+z)^2/(1 + z + 0.5*z^2) does not match V_3's derived dynamics. V_3 dynamics naturally reduce to LCDM at H_0=73 (because A is confined at A_0 and V_3(A_0) acts as cosmological constant). This suggests the closed-form was a *conjectured* H(z) form not actually derived from V_3 — and that the framework's 'intrinsic H(z)' under V_3 IS essentially LCDM at H_0=73.

**This shifts the narrative on the bridge term.**

If V_3's intrinsic H(z) = LCDM-at-73, then the distance reframe ('STAM intrinsic is flatter') requires the closed-form ansatz to be retained as a *separate* piece of physics — not derived from V_3 dynamics. The closed-form had been informally connected to the bridge-term / photon-A-traversal story; under V_3, those threads need to be re-examined. The bridge term may sit on a different layer than V_3 modified Friedmann.

## What this implies for next work

- V_3 modified Friedmann gives essentially LCDM H(z) when calibrated to Omega_DE today. The framework's cosmic expansion dynamics are not where the distance reframe lives.
- The closed-form H(z) was probably an ad-hoc shape, not derived from V_3. Where did it come from? Worth tracing back to its origin (likely from inverting the bz-form distance modulus). If the derivation is unsound, the closed-form should not be load-bearing.
- The photon-A-traversal mechanism (light traversing A_0 in voids accumulates extra delay) is a *separate* effect from V_3 cosmic expansion. It affects observed distance modulus, not H(z). Chronometer H(z) and SN distance modulus probe different layers of the framework.
- The bridge term's identity becomes clearer: it's about photon paths through the cosmic A_0 field, NOT about cosmic expansion dynamics. V_3 Friedmann doesn't constrain it.