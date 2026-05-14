#!/usr/bin/env python3
"""
G72a_static_aligned_aether_check.py

Step A of the Einstein-aether track.  Checks whether the framework's
committed metric is in the static-aligned-aether solution class.

Static-aligned Einstein-aether (u^mu aligned with timelike Killing vector)
admits only Schwarzschild-(de Sitter) static spherical solutions
(Eling-Jacobson 2006).  In that solution class, the Einstein tensor
satisfies

    G^t_t = G^r_r = G^theta_theta   (proportional to delta^mu_nu)

i.e., G_munu is a multiple of g_munu (cosmological-constant-like).

If our framework's G_munu does NOT satisfy this, the framework cannot
be embedded in static-aligned Einstein-aether (with or without V(A)
potential).  Then we must escalate to tilted aether (u^r != 0).
"""

from __future__ import annotations

import sys
import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

r, M = sp.symbols('r M', positive=True)

# Mixed Einstein-tensor components from G69 (M = G = c = 1 units)
Gtt = -16 * (r - 3*M)**3 * (720*M**3 - 648*M**2*r + 117*M*r**2 + 13*r**3) / r**8
Grr = -16 * (13*r - 24*M) * (r - 3*M)**4 / r**7
Gthth = -1440 * M * (r - 3*M)**3 * (r - 2*M) * (r - M) / r**8

print("Schwarzschild-de Sitter check: G^t_t = G^r_r = G^th_th ?")
print()
print("G^t_t - G^r_r =")
d1 = sp.simplify(Gtt - Grr)
print(f"  {sp.factor(d1)}")
print()
print("G^t_t - G^th_th =")
d2 = sp.simplify(Gtt - Gthth)
print(f"  {sp.factor(d2)}")
print()
print("G^r_r - G^th_th =")
d3 = sp.simplify(Grr - Gthth)
print(f"  {sp.factor(d3)}")
print()

# Numerical samples inside the shell
print("Numerical sample at r = 2.5 M (inside the final shell):")
sub = {M: 1, r: sp.Rational(5, 2)}
print(f"  G^t_t   = {float(Gtt.subs(sub)):.6f}")
print(f"  G^r_r   = {float(Grr.subs(sub)):.6f}")
print(f"  G^th_th = {float(Gthth.subs(sub)):.6f}")
print()

if d1 == 0 and d2 == 0:
    print("VERDICT: G^mu_nu IS proportional to delta^mu_nu -> Schwarzschild-de Sitter form.")
    print("Static-aligned aether + V(A) cosmological term could match.")
else:
    print("VERDICT: G^mu_nu is NOT proportional to delta^mu_nu.")
    print("Framework's metric is NOT in the Schwarzschild-(de Sitter) class.")
    print("Therefore static-aligned Einstein-aether + any V(A) CANNOT reproduce it.")
    print("Must escalate to tilted aether (u^r != 0).")
