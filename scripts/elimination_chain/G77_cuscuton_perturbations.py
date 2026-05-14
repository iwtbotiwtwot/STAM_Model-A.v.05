#!/usr/bin/env python3
"""
G77_cuscuton_perturbations.py

Linear perturbations of the cuscuton-embedded committed metric.

Question: does the cuscuton-style square-root kinetic L_cusc = -sigma sqrt((grad A)^2)
actually eliminate the scalar ghost in the scalar (l = 0) perturbation
sector, on a SPACELIKE-gradient background (static A(r))?

The original Afshordi-Chung cuscuton non-propagation result is for a
TIMELIKE gradient (clock field).  For spacelike gradient, the situation
can differ -- the kinetic Lagrangian's degeneracy is along the gradient
direction (radial here), but TIME and ANGULAR perturbations may still
carry kinetic structure.

This script computes the quadratic action S^(2)[delta A] from L_cusc
explicitly:

  A = A_bg(r) + delta A(t, r, theta, phi)
  (grad A)^2 = X + 2 Y + (grad delta A)^2,
    where X = (grad A_bg)^2 = k (A_bg')^2,
          Y = g^{mu nu} d_mu A_bg d_nu delta A
            = k A_bg' d_r delta A   for static A_bg(r)
  s = sqrt((grad A)^2)

  s^(2) (quadratic part of s in delta A)
       = (1/(2 s_bg))[(grad delta A)^2 - Y^2 / s_bg^2]

The bracket evaluates to:
  (grad delta A)^2 - Y^2 / s_bg^2
    = -(d_t delta A)^2/h + k(d_r delta A)^2 + (angular)/r^2
      - k(d_r delta A)^2
    = -(d_t delta A)^2/h + (angular)/r^2

So the RADIAL gradient piece (d_r delta A)^2 cancels -- this is the
cuscuton degeneracy along the gradient direction.  But TIME and
ANGULAR pieces survive.

Therefore L_cusc^(2) gives:
  L^(2) = -sigma_bg * s^(2) (+ sigma'(A) delta A * s^(1) terms)
        = (sigma_bg / (2 s_bg h)) (d_t delta A)^2
          - (sigma_bg / (2 s_bg r^2)) (angular)
        + (other terms involving delta A * d delta A)

Coefficient of (d_t delta A)^2: sigma_bg / (2 s_bg h).
For no ghost in the temporal direction, this MUST be positive
(matching the canonical sign +1/(2h) for standard scalar kinetic).

Since s_bg > 0 and h > 0, sign(K_t) = sign(sigma_bg).
For our committed metric, sigma_bg < 0 (verified in G76).
THEREFORE: K_t < 0, the temporal direction carries a GHOST.

This script verifies this numerically and quantifies the situation.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

# Symbolic
A = sp.symbols('A', positive=True)
M_sym = sp.symbols('M', positive=True)
y = 3*A - 2
F = 1 - 5*y**4 + 4*y**5
Q3 = 108*A**3 - 189*A**2 + 114*A - 23
P8 = (393660*A**8 - 2318220*A**7 + 5730669*A**6 - 7684146*A**5
      + 5982660*A**4 - 2636280*A**3 + 551025*A**2 - 8580*A - 10790)

# G76: sigma(A)/f(A) = y^2 P_8 sqrt(1-A) / (M F^(3/2))
sigma_over_f = y**2 * P8 * sp.sqrt(1 - A) / (M_sym * F**sp.Rational(3, 2))

# s(A) = sqrt(k (A')^2) = sqrt((1-A) F) * |A'| with A' = -A^2/(2M):
# s = (A^2/(2M)) * sqrt((1-A) F)
s_expr = (A**2 / (2 * M_sym)) * sp.sqrt((1 - A) * F)

# Metric h, k
h_expr = 1 - A
k_expr = (1 - A) * F

# From G70:
dlnf_dA = -2 * y**3 * (45*A**2 - 33*A - 13) / (A * F)
Z_over_f = 2 * y**2 * P8 / (81 * A**2 * (1 - A)**4 * Q3**2)

# Lambdify
sigma_over_f_M1 = sigma_over_f.subs(M_sym, 1)
s_M1 = s_expr.subs(M_sym, 1)
h_fn = sp.lambdify(A, h_expr, 'numpy')
k_fn = sp.lambdify(A, k_expr, 'numpy')
sigma_over_f_fn = sp.lambdify(A, sigma_over_f_M1, 'numpy')
s_fn = sp.lambdify(A, s_M1, 'numpy')
dlnf_dA_fn = sp.lambdify(A, dlnf_dA, 'numpy')
Z_over_f_fn = sp.lambdify(A, Z_over_f, 'numpy')


# Integrate f(A) from f(2/3) = 1
def rhs(A_val, lnf):
    return dlnf_dA_fn(A_val)


A0_int = 2/3 + 1e-6
A_end_int = 1 - 1e-4
sol = solve_ivp(rhs, (A0_int, A_end_int), [0.0], dense_output=True,
                rtol=1e-10, atol=1e-13, max_step=1e-4)

A_grid = np.linspace(A0_int, A_end_int, 2000)
lnf = sol.sol(A_grid).flatten()
f_vals = np.exp(lnf)
sigma_vals = f_vals * sigma_over_f_fn(A_grid)
s_vals = s_fn(A_grid)
h_vals = h_fn(A_grid)


# Coefficient of (d_t delta A)^2 in L_cusc^(2):
# K_t = sigma_bg / (2 * s_bg * h)
K_t = sigma_vals / (2 * s_vals * h_vals)

# Coefficient of (1/r^2)(d_theta delta A)^2 in L_cusc^(2):
# K_a = -sigma_bg / (2 * s_bg)   (the 1/r^2 is the angular metric piece)
K_a = -sigma_vals / (2 * s_vals)

# For comparison: standard canonical scalar's K_t would be +1/(2h) (positive).
K_t_canonical = 1 / (2 * h_vals)

print("=" * 80)
print("LINEAR PERTURBATIONS OF CUSCUTON-EMBEDDED METRIC")
print("=" * 80)
print()
print("Quadratic action coefficients for delta A perturbations:")
print()
print(f"  {'A':>8}{'r':>8}{'sigma':>14}{'K_t (d_t delta A)^2':>24}{'K_a (d_theta dA)^2':>22}{'no-ghost?':>14}")
print("-" * 90)
for a in [0.67, 0.70, 0.75, 0.80, 0.85, 0.90, 0.953, 0.99]:
    idx = np.argmin(np.abs(A_grid - a))
    sig = sigma_vals[idx]
    kt = K_t[idx]
    ka = K_a[idx]
    no_ghost = (kt > 0)
    flag = "OK" if no_ghost else "GHOST"
    print(f"  {a:>8.3f}{2/a:>8.3f}{sig:>14.4e}{kt:>24.4e}{ka:>22.4e}{flag:>14}")
print()

# Sign verdict
neg_Kt = np.sum(K_t < 0)
print(f"K_t < 0 (TEMPORAL GHOST) at {neg_Kt}/{len(K_t)} grid points ({100*neg_Kt/len(K_t):.2f}%)")
neg_sigma = np.sum(sigma_vals < 0)
print(f"sigma < 0 at {neg_sigma}/{len(sigma_vals)} grid points ({100*neg_sigma/len(sigma_vals):.2f}%)")
print()

# Comparison to G71's G_ghost (single-scalar)
Z_vals = f_vals * Z_over_f_fn(A_grid)
G_ghost_single = 3 * (f_vals * dlnf_dA_fn(A_grid))**2 + 2 * Z_vals * f_vals
neg_G = np.sum(G_ghost_single < 0)
print(f"G71 G_ghost < 0 (G70 scalar-tensor radial ghost) at "
      f"{neg_G}/{len(G_ghost_single)} ({100*neg_G/len(G_ghost_single):.2f}%)")
print()

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
ax.plot(A_grid, K_t, 'tab:red', linewidth=2, label=r'$K_t$ (cuscuton $(\partial_t \delta A)^2$ coef)')
ax.plot(A_grid, K_t_canonical, 'tab:green', linewidth=1.5, linestyle='--',
        label=r'$+1/(2h)$ (canonical, healthy)')
ax.axhline(0, color='k', linewidth=1, linestyle='--')
ax.axvline(2/3, color='gray', linestyle='--', alpha=0.5)
ax.axvline(1.0, color='black', linestyle=':', alpha=0.5)
ax.set_xlabel('A')
ax.set_ylabel(r'Coefficient of $(\partial_t \delta A)^2$')
ax.set_title('Cuscuton temporal kinetic coefficient\n(negative = temporal ghost)')
ax.set_yscale('symlog', linthresh=1e-2)
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1]
ax.plot(A_grid, G_ghost_single, 'tab:red', linewidth=2,
        label=r'G70 scalar-tensor $G_{\rm ghost}$ (radial)')
ax.plot(A_grid, K_t * 2 * h_vals * s_vals, 'tab:purple', linewidth=2,
        linestyle='--',
        label=r'Cuscuton $\sigma_{bg}$ = K_t·2sh (temporal)')
ax.axhline(0, color='k', linewidth=1, linestyle='--')
ax.axvline(2/3, color='gray', linestyle='--', alpha=0.5)
ax.axvline(1.0, color='black', linestyle=':', alpha=0.5)
ax.set_xlabel('A')
ax.set_ylabel('ghost coefficient')
ax.set_title('Ghost has moved from radial -> temporal direction')
ax.set_yscale('symlog', linthresh=1e-2)
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
out = PLOTS / "G77_cuscuton_perturbations.png"
plt.savefig(out, dpi=200)
plt.close()
print(f"Plot saved: {out}")
print()

# Verdict
print("=" * 80)
print("VERDICT")
print("=" * 80)
print()
print("The cuscuton-style embedding S = (1/16piG) integral sqrt(-g) [f R - 2 sigma s - 2 V]")
print("does NOT eliminate the scalar ghost for our committed metric.")
print()
print("Reason: the original Afshordi-Chung cuscuton non-propagation property")
print("applies to TIMELIKE gradient (clock field).  Our background has SPACELIKE")
print("gradient (static A(r)).  In the spacelike case:")
print()
print("  - The RADIAL direction (along grad A_bg) is degenerate: (d_r delta A)^2")
print("    cancels in L^(2).  No radial scalar propagation.  This DOES remove")
print("    the radial scalar ghost that G71 identified.")
print()
print("  - But the TEMPORAL direction (d_t delta A)^2 retains a kinetic coefficient")
print("    K_t = sigma_bg / (2 s_bg h).  For our framework, sigma_bg < 0 throughout")
print("    the shell, so K_t < 0 -> TEMPORAL GHOST.")
print()
print("So the ghost has MOVED from the radial direction (G70/G71 scalar-tensor)")
print("to the temporal direction (G76 cuscuton).  It hasn't been eliminated.")
print()
print("To genuinely eliminate the scalar DOF in a continuum-field Lagrangian,")
print("we'd need:")
print("  - TIMELIKE-gradient clock field + A as derived (mimetic-style)")
print("  - Full Lagrange-multiplier constraint pinning delta A algebraically")
print("    to delta g (not just (grad A)^2 = W(A) which leaves time-dependence")
print("    of delta A free)")
print("  - DHOST with degenerate higher-derivative structure")
print("  - Or: accept the substance-ontology reformulation (Open Problem #6's")
print("    primary closure) as the foundational answer.")
print()
print("Under substance ontology, A is NOT a continuum field with Cauchy data;")
print("delta A perturbations are not the right kinematic objects to analyze.")
print("The Lagrangian-level ghost is a continuum-field artifact of trying to")
print("force substance-ontology dynamics into a Lagrangian field-theory mold.")
