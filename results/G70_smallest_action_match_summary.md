# G70 — Smallest covariant action match

**Date: 2026-05-13.**  Step 2 of the Lagrangian-for-A track (Open Problem #1).

## Result: Brans-Dicke (f, V) is INSUFFICIENT.  Smallest match is scalar-tensor.

Smallest covariant level-set / anisotropic-shell action reproducing the framework's exact T^eff_munu on shell:

```
S = (1 / 16 pi G) * integral d^4 x sqrt(-g) [ f(A) R - Z(A) (grad A)^2 - 2 V(A) ]
```

with three free functions of A determined uniquely (up to one integration constant on f) by the EOM:

### f(A) — non-minimal coupling

```
d(ln f)/dA = -2 y^3 (45 A^2 - 33 A - 13) / [ A F(y) ]
         y = 3 A - 2
       F(y) = 1 - 5 y^4 + 4 y^5
            = (1-y)^2 (4 y^3 + 3 y^2 + 2 y + 1)
            = 9 (1-A)^2 (108 A^3 - 189 A^2 + 114 A - 23)
```

Integration constant fixes f(A=2/3) = 1 (GR matching at PS).

### Z(A) — kinetic norm

```
Z(A) / f(A) = 2*(3*A - 2)**2*(393660*A**8 - 2318220*A**7 + 5730669*A**6 - 7684146*A**5 + 5982660*A**4 - 2636280*A**3 + 551025*A**2 - 8580*A - 10790)/(81*A**2*(A - 1)**4*(108*A**3 - 189*A**2 + 114*A - 23)**2)
```

### V(A) — potential

```
V(A) / f(A) = A**2*(3*A - 2)**2*(54675*A**7 - 311283*A**6 + 709074*A**5 - 831546*A**4 + 531819*A**3 - 177021*A**2 + 24258*A + 26)/(36*M**2*(A - 1)*(108*A**3 - 189*A**2 + 114*A - 23))
```

## Counting check

- 3 independent T^eff components on static spherical bg: rho, p_r, p_t

- 1 Bianchi (conservation): leaves 2 free of r

- 3 action functions f(A), Z(A), V(A) -> evaluated on background, 3 functions of r

- 1 background constraint: A = 2M/r determines A(r)

- Net: 3 unknowns matching 3 EOM components -> overdetermined system that closes uniquely.  Consistency is the non-trivial content.

## Why Brans-Dicke failed

With only f(A), V(A), the (tt - thth) and (rr - thth) difference equations give two ODEs for f(A) alone.  These two ODEs are *not* compatible for the committed quintic metric: their residual at A = 0.85 (M = 1) is ~1.4 -- of the same magnitude as G_munu itself.  Adding the Z(A) kinetic term inserts an extra knob that appears only in the (rr - thth) channel (because Z d_mu A d_nu A is purely radial on a static spherical bg), absorbing the inconsistency.

## What this gives the next step

- **Open Problem #1 (Lagrangian for A):**  the smallest covariant action with A as a level-set scalar is identified.  f(A), Z(A), V(A) all closed-form in y = 3A - 2 and F(y).

- **Open Problem #6 (ghost-freedom):**  the kinetic operator's sign is now Z(A).  Sign of Z(A)/f(A) across the final shell determines whether the A-mode is healthy (Z > 0) or ghostly (Z < 0).  G71 should plot Z(A)/f(A) numerically across the shell and report the sign behavior.

## Files

- [scripts/G70_smallest_action_match.py](../scripts/G70_smallest_action_match.py)
