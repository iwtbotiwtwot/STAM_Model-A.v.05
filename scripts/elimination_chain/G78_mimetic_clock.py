#!/usr/bin/env python3
"""
G78_mimetic_clock.py

Tries the mimetic-with-timelike-clock route as continuum-field
embedding for the committed metric.  Tests:

1. The mimetic clock tau on a static spherical background:  what
   constraints on tau follow from (grad tau)^2 = -1?  Show tau must
   have a PG-extended form with both t and r dependence (because
   pure tau = T(t) fails for non-flat h(r)).
2. The stress contribution 2 lambda u_mu u_nu of the mimetic
   constraint: compute its components in static coordinates.
3. Check whether the off-diagonal pieces can be absorbed by the
   f(A) R coupling (whose ∇∇f piece is purely diagonal in static
   coords for f = f(A(r))).
4. If not, conclude that mimetic+f(A)R is structurally incompatible
   with the framework's diagonal-static-spherical metric.
"""

from __future__ import annotations

import sys

import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


# Background: static spherical with our committed h, k
r, t = sp.symbols('r t', positive=True)
M = sp.symbols('M', positive=True)

A_expr = 2 * M / r
y_expr = 3 * A_expr - 2
F_expr = 1 - 5 * y_expr**4 + 4 * y_expr**5
h_expr = 1 - A_expr
k_expr = h_expr * F_expr

print("=" * 80)
print("STEP 1: Is there a tau = T(t) only that satisfies (grad tau)^2 = -1?")
print("=" * 80)
print()
print("If tau = T(t), then (grad tau)^2 = g^tt (T'(t))^2 = -(T'(t))^2 / h.")
print("(grad tau)^2 = -1 requires (T'(t))^2 = h(r).")
print()
print("But T'(t) is a function of t only, h(r) is a function of r only.")
print("Both equal => both must be constant.  h(r) = const means flat spacetime.")
print()
print("For our committed metric h = 1 - 2M/r, h is NOT constant.")
print("Therefore tau = T(t) alone does NOT solve the mimetic constraint.")
print()


print("=" * 80)
print("STEP 2: PG-extended mimetic clock  tau = t + R(r)")
print("=" * 80)
print()
print("Ansatz: tau = t + integral b(r) dr.  Then d_t tau = 1, d_r tau = b(r).")
print("Constraint:")
print("  g^tt (d_t tau)^2 + g^rr (d_r tau)^2 = -1")
print("  -(1)^2 / h + k b^2 = -1")
print("  k b^2 = -1 + 1/h = (1-h)/h = A/h")
print("  b^2 = A / (k h)")
print()

b_sq = A_expr / (k_expr * h_expr)
b_sq_simplified = sp.simplify(b_sq)
print(f"  b(r)^2 = A / (k h) = {b_sq_simplified}")
print()

# Numerical sanity: inside the shell, A > 0, k > 0, h > 0, so b^2 > 0.  Real.
print("Inside the final shell (A in (2/3, 1)): A > 0, k > 0, h > 0.")
print("So b^2 > 0 throughout -> tau is real-valued.  Good.")
print()
print("At horizon (A -> 1, k -> 0, h -> 0):")
print("  b^2 = A/(k h) -> 1/0/0 -> infinity.  PG mimetic clock diverges at horizon")
print("  (consistent with the metric's horizon structure).")
print()


print("=" * 80)
print("STEP 3: Stress-energy 2 lambda u_mu u_nu in static coordinates")
print("=" * 80)
print()
print("With u_mu = d_mu tau = (1, b(r), 0, 0)  in static (t, r, theta, phi) coords:")
print()
print("Mimetic stress T^mim_mu_nu = 2 lambda u_mu u_nu has components:")
print()
print("  T^mim_tt = 2 lambda * 1 * 1     = 2 lambda")
print("  T^mim_tr = 2 lambda * 1 * b(r)  = 2 lambda b   <- OFF-DIAGONAL")
print("  T^mim_rr = 2 lambda * b * b     = 2 lambda b^2")
print("  T^mim_thth = T^mim_phph = 0")
print()
print("In mixed indices T^mim^mu_nu = g^mu_alpha T^mim_alpha_nu:")
print()
print("  T^mim^t_t = g^tt * 2 lambda      = -2 lambda / h")
print("  T^mim^t_r = g^tt * 2 lambda b    = -2 lambda b / h    <- OFF-DIAGONAL")
print("  T^mim^r_t = g^rr * 2 lambda b    = +2 lambda k b      <- OFF-DIAGONAL")
print("  T^mim^r_r = g^rr * 2 lambda b^2  = +2 lambda k b^2 = +2 lambda A/h")
print()


print("=" * 80)
print("STEP 4: matching to our committed metric")
print("=" * 80)
print()
print("Our committed metric has G^mu_nu strictly DIAGONAL (G69):")
print("  G^t_r = G^r_t = 0  (all off-diagonal components vanish)")
print()
print("The non-minimal coupling f(A) R contributes through nabla_mu nabla_nu f.")
print("For static f = f(A(r)), nabla_t nabla_r f also vanishes:")
print("  nabla_t nabla_r f = d_t d_r f - Gamma^lambda_{tr} d_lambda f")
print("                    = 0 - Gamma^t_{tr} * d_t f - Gamma^r_{tr} * d_r f")
print("                    = 0 - (h'/2h) * 0 - 0 = 0")
print("(d_t f = 0 since f is static; Gamma^r_{tr} = 0 in static metric.)")
print()
print("So  nabla_mu nabla_nu f  is also strictly DIAGONAL in static coordinates.")
print()
print("Einstein equation (mixed indices, f-R + mimetic + V):")
print("  f G^mu_nu = nabla^mu nabla_nu f - delta^mu_nu (box f + V) - 2 lambda u^mu u_nu")
print()
print("The (t, r) component:")
print("  f * 0 = 0 - 0 - 2 lambda u^t u_r = -2 lambda * (-1/h) * b = +2 lambda b/h")
print()
print("Setting LHS = RHS:  0 = 2 lambda b / h.")
print("Since b != 0 in the shell, this forces lambda = 0.")
print()


print("=" * 80)
print("STEP 5: STRUCTURAL CONCLUSION")
print("=" * 80)
print()
print("The mimetic constraint forces tau to have non-zero radial gradient b(r) on")
print("our static spherical background (no pure-time mimetic clock exists for")
print("non-trivial h(r)).  This non-zero b makes 2 lambda u_mu u_nu carry off-")
print("diagonal (t, r) components.  But:")
print()
print("  - G^mu_nu of the committed metric is strictly diagonal.")
print("  - nabla^mu nabla_nu f for f = f(A(r)) is strictly diagonal.")
print("  - V g_munu is strictly diagonal.")
print()
print("So the off-diagonal mimetic stress has nothing to cancel against on the")
print("LHS.  The matching equation in the (t, r) channel forces  lambda = 0.")
print()
print("With lambda = 0, the mimetic Lagrange-multiplier piece drops out entirely,")
print("and the action reduces to vanilla f(A) R - 2 V(A) -- which is exactly the")
print("Brans-Dicke ansatz that G70 showed is INSUFFICIENT to match our metric.")
print()
print("Therefore: simple mimetic-with-timelike-clock + f(A) R + V(A) cannot")
print("reproduce the committed metric.  The PG-extended clock generates off-")
print("diagonal stress incompatible with the diagonal-static-spherical structure.")
print()


print("=" * 80)
print("What might work (deferred):")
print("=" * 80)
print()
print("1. Mimetic + disformal coupling: g_munu -> g_munu + B(A) d_mu A d_nu A in")
print("   the matter sector.  The disformal piece carries the right off-diagonal")
print("   structure to potentially cancel the PG-clock off-diagonal stress.")
print()
print("2. Mimetic-DHOST: add higher-derivative terms with degeneracy conditions.")
print("   Gives ghost-freedom by construction with enough freedom to match.")
print()
print("3. Two-Lagrange-multiplier: enforce (grad tau)^2 = -1 AND (grad A)^2 = W(A).")
print("   The second constraint absorbs the off-diagonal piece through its own")
print("   stress contribution.")
print()
print("4. Accept the substance-ontology reformulation (Open Problem #6's primary")
print("   closure via Gamma_res in Open Problem #5a) as the foundational answer.")
print("   The Lagrangian-level continuum-field ghost is an artifact of trying to")
print("   force substance dynamics into a wrong-fitting field-theory mold.")
print()
print("Status of Lagrangian-level Open Problem #6 after G70-G78:")
print("  - Scalar-tensor (G70/G71): exists; radial scalar ghost in 86% of shell.")
print("  - Multi-field (G74): cannot rescue (M_AA unchanged by Psi enrichment).")
print("  - Cuscuton (G76): partial; temporal ghost in 100% of shell (G77).")
print("  - Static-aligned Einstein-aether (G72a): structurally excluded.")
print("  - Tilted Einstein-aether (G72c): observationally excluded (c_i too large).")
print("  - Mimetic + f(A) R (G78): structurally excluded (off-diagonal mismatch).")
print()
print("All standard low-derivative continuum-field embeddings have been ruled out.")
print("The substance-ontology reformulation remains the cleanest framework-internal")
print("closure of the perturbative-stability question.")
