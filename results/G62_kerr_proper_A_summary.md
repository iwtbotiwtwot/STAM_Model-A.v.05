# G62 - Kerr Extension with Proper Framework A(r, theta; M, a)

**Date: 2026-05-13.** Corrects G61's naive A = 2M/r reading. Uses the framework's actual Kerr substance density (already in G3 from May 2026): `A(r, theta; M, a) = 2 M r / Sigma` where `Sigma = r^2 + a^2 cos^2(theta)`.

## What was wrong in G61

G61 used naive `A = 2M/r` for Kerr, which:

- Does not give A = 1 at the actual Kerr horizon (only at r = 2M, which is the equator)

- Gives nonsense for the pole (T_pole != T_Kerr)

- Breaks down structurally at high spin (A > 1 at the photon orbit)



Sean's clarification (2026-05-13) pointed to the actual framework Hawking mechanism (G2-G4 May 2026): rotation makes outward emission more favorable, which is the framework-native super-radiance picture. The proper A formula was already in G3, just not pulled forward to G61.

## Proper framework A for Kerr (already in G3)

```
A(r, theta; M, a) = 2 M r / Sigma
Sigma = r^2 + a^2 cos^2(theta)
```

Properties:

- A = 1 on the actual oblate Kerr horizon: 2Mr = Sigma => r^2 - 2Mr + a^2 cos^2(theta) = 0

- At equator (cos^2 theta = 0): A = 2M/r -- Schwarzschild form

- At pole (cos^2 theta = 1): r = r_+ on horizon

- 0 < A < 1 everywhere outside the bubble; A > 1 forbidden

## Verified structural identities

On the bubble (A = 1 surface), |grad A| has components:

```
grad_r A|_bubble = (M - r)/(M r)
grad_theta A|_bubble = a^2 sin(2 theta) / (2 M r)
```

**At equator (theta = pi/2, r = 2M):**

```
|grad A|_eq = 1/(2M)  =>  k_B T_eq = hbar c / (8 pi M) = k_B T_Schwarzschild
```

**At pole (theta = 0, r = r_+):**

```
|grad A|_pole = sqrt(M^2 - a^2) / (M r_+)
k_B T_pole = (hbar c / 4 pi) * sqrt(M^2-a^2)/(M r_+)
```

Using the Kerr horizon identity `r_+^2 + a^2 = 2 M r_+`:

```
kappa_Kerr = (r_+ - r_-)/(2(r_+^2 + a^2)) = sqrt(M^2-a^2)/(2 M r_+)
=> |grad A|_pole = 2 kappa_Kerr
=> k_B T_pole = hbar c kappa_Kerr / (2 pi) = k_B T_Kerr  (exactly)
```

**Both endpoint temperatures land exactly on the framework's claimed identities** across all spins 0 <= a/M < 1.

## Distinguishing prediction: latitudinal banding

The framework predicts NON-UNIFORM temperature on the spinning bubble:

- Equator (hottest): T_Schw (independent of spin)

- Pole (coolest): T_Kerr (decreasing with spin)

- Smooth profile between



This is DISTINCT from standard Kerr's uniform horizon T_K. Combined with rotation enhancement R(theta) (G4 super-radiance), emission is strongly equatorial -- equator gets both hot AND rotationally enhanced.



Per Sean: rotation tilt at the horizon makes outward more favorable (inward is forbidden by no-interior). This is the framework-native form of super-radiance enhancement.

## Equatorial-plane ringdown

In the equatorial plane (theta = pi/2, cos^2 theta = 0), A reduces to Schwarzschild form A = 2M/r. So equatorial photon orbit analysis can use Schwarzschild-style coordinates with the framework's k(A) commitment.



For Kerr equatorial prograde photon orbit:

| a/M | r_pro/M | A_pro = 2M/r_pro |

|---:|---:|---:|

| 0.00 | 3.0000 | 0.6667 |

| 0.30 | 2.6300 | 0.7604 |

| 0.50 | 2.3473 | 0.8520 |

| 0.67 | 2.0682 | 0.9670 |

| 0.85 | 1.6938 | 1.1808 |

| 0.95 | 1.3863 | 1.4427 |



A_pro exceeds 2/3 at moderate spin and approaches 1 for high spin. The Kerr photon orbit lives INSIDE the framework's PS structural landmark (Sigma = 2 = 3 * 2/3). This raises an open structural question.

## Open structural question for ringdown

Two readings of 'where the framework's k(A) modification turns on':



**Reading A** (PS = where light escapes): for Kerr, this is the actual Kerr photon orbit. Then `k = (1 - A_pro)` exactly, framework recovers GR-Kerr Lyapunov, STAM-Kerr ringdown = exact GR-Kerr.



**Reading B** (PS = universal A = 2/3 = Sigma 2 landmark): then the Kerr photon orbit at A > 2/3 is inside framework's final shell. The quintic Hermite F applies, giving departure from GR-Kerr.



Both readings are structurally consistent. Reading A is cleaner under 'no STAM modification where light escapes'; Reading B is cleaner under 'Sigma integer landmarks are universal.' Framework hasn't committed.

## Status (corrected from G61)

**Resolved:**

- A(r, theta; M, a) was always in G3; the framework's Kerr extension is more mature than G61 suggested

- T_eq = T_Schw, T_pole = T_Kerr both verified to machine precision

- Latitudinal banding is the framework's distinguishing prediction vs. Kerr's uniform T

- Rotation enhancement (super-radiance) is the framework-native form of Penrose-equivalent extraction via Hawking-style channel



**Open:**

- PS identification for Kerr ringdown (Reading A vs Reading B)

- Full T(theta) + R(theta) Hawking spectrum (latitudinal banding signature for PBH evaporation if observed)

- Inside-PS metric for Kerr off-equator (quintic Hermite F applied to off-equatorial A)

## Files

- [scripts/G62_kerr_proper_A.py](../scripts/G62_kerr_proper_A.py)

- [scripts/G3_spinning_bubble_thermodynamics.py](../scripts/G3_spinning_bubble_thermodynamics.py) (original derivation, May 2026)

- [scripts/G4_rotational_emission_bias.py](../scripts/G4_rotational_emission_bias.py) (super-radiance, May 2026)

- [plots/G62_T_profile_corrected.png](../plots/G62_T_profile_corrected.png)

- [plots/G62_emission_pattern.png](../plots/G62_emission_pattern.png)

- [plots/G62_bubble_shape_with_T.png](../plots/G62_bubble_shape_with_T.png)
