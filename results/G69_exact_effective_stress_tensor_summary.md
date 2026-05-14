# G69 — Exact Effective Stress Tensor (closed form)

**Date: 2026-05-13.**  Step 1 of the ghost-freedom / Lagrangian-for-A track (Open Problems #6 and #1).  Computes the exact symbolic Einstein tensor G_munu for the framework's committed strong-field metric and defines T^eff_munu = (c^4 / (8 pi G)) G_munu.

## Metric (spinless)

```
ds^2 = -(1 - A) c^2 dt^2 + dr^2 / k(A) + r^2 dOmega^2
A(r) = 2 G M / (c^2 r)
k(A) = (1 - A) F(y),   y = 3 A - 2
F(y) = 1                       for A <= 2/3   (outside PS)
F(y) = 1 - 5 y^4 + 4 y^5       for 2/3 < A < 1 (inside PS, quintic Hermite)
```

## General formula  (any h(r), k(r))

Mixed Einstein-tensor components for ds^2 = -h dt^2 + dr^2/k + r^2 dOmega^2:

```
G^t_t       = (r*Derivative(k(r), r) + k(r) - 1)/r**2
G^r_r       = (r*k(r)*Derivative(h(r), r) + (k(r) - 1)*h(r))/(r**2*h(r))
G^theta_th  = k(r)*Derivative(h(r), (r, 2))/(2*h(r)) + Derivative(h(r), r)*Derivative(k(r), r)/(4*h(r)) - k(r)*Derivative(h(r), r)**2/(4*h(r)**2) + Derivative(k(r), r)/(2*r) + k(r)*Derivative(h(r), r)/(2*r*h(r))
```

Equivalently, in fluid form (8 pi rho = -G^t_t etc., G = c = 1):

```
8 pi rho   = (-r*Derivative(k(r), r) - k(r) + 1)/r**2
8 pi p_r   = (r*k(r)*Derivative(h(r), r) + (k(r) - 1)*h(r))/(r**2*h(r))
8 pi p_t   = k(r)*Derivative(h(r), (r, 2))/(2*h(r)) + Derivative(h(r), r)*Derivative(k(r), r)/(4*h(r)) - k(r)*Derivative(h(r), r)**2/(4*h(r)**2) + Derivative(k(r), r)/(2*r) + k(r)*Derivative(h(r), r)/(2*r*h(r))
```

## Committed metric (inside PS, A = 2M/r, quintic Hermite F)

Mixed Einstein-tensor components (M = G = c = 1):

```
G^t_t       = -16*(-3*M + r)**3*(720*M**3 - 648*M**2*r + 117*M*r**2 + 13*r**3)/r**8

G^r_r       = -16*(-24*M + 13*r)*(-3*M + r)**4/r**7

G^theta_th  = -1440*M*(-3*M + r)**3*(-2*M + r)*(-M + r)/r**8
```

Effective fluid decomposition  T^eff^mu_nu = (1/8pi) G^mu_nu :

```
rho_eff = 2*(-19440*M**6 + 36936*M**5*r - 27135*M**4*r**2 + 9360*M**3*r**3 - 1350*M**2*r**4 + 13*r**6)/(pi*r**8)

p_r_eff = 2*(1944*M**5 - 3645*M**4*r + 2700*M**3*r**2 - 990*M**2*r**3 + 180*M*r**4 - 13*r**5)/(pi*r**7)

p_t_eff = 180*M*(54*M**5 - 135*M**4*r + 126*M**3*r**2 - 56*M**2*r**3 + 12*M*r**4 - r**5)/(pi*r**8)
```

## Energy-condition combinations (exact)

```
rho + p_r            = 360*M*(-3*M + r)**3*(-2*M + r)**2/(pi*r**8)
rho + p_t            = 2*(-3*M + r)**4*(-180*M**2 + 66*M*r + 13*r**2)/(pi*r**8)
rho + p_r + 2 p_t    = -360*M**2*(-3*M + r)**3*(-2*M + r)/(pi*r**8)
trace T              = -4*(-3*M + r)**3*(540*M**3 - 558*M**2*r + 117*M*r**2 + 13*r**3)/(pi*r**8)
```

## Bianchi / TOV verification

The TOV equation derived from ∇_µ G^µ_r = 0 is:

```
d p_r / dr + (h'/(2 h)) (rho + p_r) + (2/r) (p_r - p_t) = 0
```

Sympy evaluates the left-hand side to: `0`



Vanishing => Bianchi identity holds symbolically; conservation is automatic.



## Physical-unit form

```
T^eff_munu = (c^4 / (8 pi G)) G_munu
```

Each component above (in G=c=1 units) carries an implicit c^4/(8 pi G) prefactor for SI conversion.  Components in M = 1 units carry an extra M^{-2} which becomes (c^4 / (G M))^2 / c^4 = c^4 / (G^2 M^2) in SI.



## What this step gives the next steps

- Open Problem #6 (ghost-freedom perturbative check):  the closed-form T^eff_munu is the right-hand side any candidate action S[A, g] must reproduce on shell.  Linearizing around this background and inspecting the quadratic kinetic operator for scalar / tensor modes is the next step.

- Open Problem #1 (Lagrangian for A):  with T^eff_munu in closed form, we can ask which scalar-tensor / non-metric / f(R) actions reproduce exactly these (rho, p_r, p_t) profiles on the spherical background.  G25–G27 attempted standard scalar-tensor ansaetze and found ghost regions plus R_s-dependent V; with the corrected metric we can now do the same inversion under the committed k(A).



## Files

- [scripts/G69_exact_effective_stress_tensor.py](../scripts/G69_exact_effective_stress_tensor.py)
