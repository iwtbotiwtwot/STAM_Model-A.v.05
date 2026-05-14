# G71 - Ghost-freedom check on smallest covariant action

**Date: 2026-05-13.**  Step 3 of the Lagrangian-for-A track (Open Problem #6).  Linearizes the smallest covariant action from G70 around the committed-metric / A = 2M/r background and checks the sign of the scalar-mode kinetic operator.

## Setup

Smallest covariant action (G70):

```
S = (1 / 16 pi G) integral sqrt(-g) [ f(A) R - Z(A) (grad A)^2 - 2 V(A) ] d^4 x
```

Standard scalar-tensor result: linearizing around the static spherical background, the propagating scalar mode has effective Brans-Dicke parameter

```
omega_BD(A) = Z(A) f(A) / (f'(A))^2
```

Ghost-freedom: 2 omega_BD + 3 > 0, equivalently

```
G_ghost(A) = 3 (f'(A))^2 + 2 Z(A) f(A)  >  0
```

Graviton positivity: f(A) > 0.

## Sign reduction

Using f' = f * d(ln f)/dA and Z = f * (Z/f), the overall scale of f cancels:

```
G_ghost / f^2 = (4 y^2 / [A^2 F(y)^2]) * S(A)
```

with

```
S(A) = 3 y^4 (45 A^2 - 33 A - 13)^2 + P_8(A)
```

The prefactor 4 y^2 / [A^2 F(y)^2] is non-negative inside the shell, so **sign(G_ghost) = sign(S(A))**.

## Result

**Sign of G_ghost is NOT uniform across the shell.** Zero crossing(s) at A in [0.952938].

Inner-shell region (A just inside PS, up to A = 0.952938) has G_ghost < 0 -- ghost mode in the smallest scalar-tensor embedding.

Outer-shell region (A > 0.952938, up through horizon) has G_ghost > 0 -- healthy.



Zero crossings of S(A) inside (2/3, 1):

- A = 0.9529681926  (y = 0.8589045779)



## Files

- [scripts/G71_ghost_freedom_check.py](../scripts/G71_ghost_freedom_check.py)

- [plots/G71_ghost_freedom_check.png](../plots/G71_ghost_freedom_check.png)
