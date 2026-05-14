#!/usr/bin/env python3
"""
G74_multifield_ledger.py

Route 3 of the S[A, g] derivation: two-scalar embedding via the ledger
channel Psi (Open Problem #6a).

Setup
=====
S = (1/16piG) integral sqrt(-g) [ f(A, Psi) R - Z_AA(A) (grad A)^2
                                - Z_PsiPsi(A) (grad Psi)^2
                                - 2 V(A, Psi) ] d^4 x

We work on the static spherical background with the 2026-05-13 committed
metric (G70 single-scalar match) and Psi_bg = 0 (the simplest backgrounds
that preserves G70's metric match).  With Psi_bg = 0, the background EOM
collapses to G70's, so the metric matching is automatic.

At the perturbative level, the effective scalar kinetic matrix for
(delta A, delta Psi) at the background is (standard 2-field scalar-tensor
result, e.g., De Felice-Tsujikawa for ST gravity):

    M = K + (3 / (2 f)) D D^T

where:
    K = diag(Z_AA, Z_PsiPsi)             (no kinetic mixing here)
    D = (df/dA, df/dPsi)^T              evaluated at (A_bg, Psi_bg = 0)

For Psi_bg = 0 and  f(A, Psi) = f_0(A) + g(A) * Psi:
    df/dA |_bg = f_0'(A)
    df/dPsi |_bg = g(A)

The no-ghost condition is M positive-definite:
    M_AA > 0  AND  det M > 0

For the single-scalar case (G70), we had only Z_AA, no Psi.  Sign was
sign of G_ghost = 3 (f_0')^2 + 2 Z_AA f_0  -- negative in 86% of shell.

In 2-scalar:
    det M = Z_AA * Z_PsiPsi
          + (3/(2f)) [ Z_PsiPsi * (f_0')^2 + Z_AA * g(A)^2 ]

The Z_AA * g(A)^2 piece (with Z_AA < 0 in shell) is negative.
The Z_PsiPsi * (f_0')^2 piece is positive if Z_PsiPsi > 0.
The Z_AA * Z_PsiPsi piece is negative if Z_PsiPsi > 0.

For det M > 0:
    Z_PsiPsi (f_0')^2 > -Z_AA * Z_PsiPsi (2f/3) - Z_AA g(A)^2
    -> need careful tuning.

This script:
1. Computes f_0, Z_AA = Z (from G70 ODE).
2. Evaluates G_ghost (single-scalar) across the shell.
3. Searches over (g(A), Z_PsiPsi(A)) ansatze to find configurations where
   M is positive-definite throughout the shell (2/3 < A < 1).
4. Tests two structural ansatze:
   (a) g(A) = g0 * y^n  (pair-structure-motivated)
   (b) g(A) = g0 * f_0(A)  (Psi enters as a phase / scaling of f)
5. Reports whether ghost can be rescued by realistic g(A), Z_PsiPsi(A).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# G70's closed-form expressions
A_sym = sp.symbols('A', positive=True)
y_sym = 3 * A_sym - 2
F_sym = 1 - 5 * y_sym**4 + 4 * y_sym**5
Q3_sym = 108 * A_sym**3 - 189 * A_sym**2 + 114 * A_sym - 23
P8_sym = (393660 * A_sym**8 - 2318220 * A_sym**7 + 5730669 * A_sym**6
          - 7684146 * A_sym**5 + 5982660 * A_sym**4 - 2636280 * A_sym**3
          + 551025 * A_sym**2 - 8580 * A_sym - 10790)

dlnf_dA_sym = -2 * y_sym**3 * (45 * A_sym**2 - 33 * A_sym - 13) / (A_sym * F_sym)
Z_over_f_sym = 2 * y_sym**2 * P8_sym / (81 * A_sym**2 * (1 - A_sym)**4 * Q3_sym**2)

dlnf_dA_fn = sp.lambdify(A_sym, dlnf_dA_sym, 'numpy')
Z_over_f_fn = sp.lambdify(A_sym, Z_over_f_sym, 'numpy')


# Integrate f_0(A) from f_0(2/3) = 1
def integrate_f0():
    def rhs(A, lnf):
        return dlnf_dA_fn(A)
    A0 = 2/3 + 1e-6
    A_end = 1 - 1e-4
    sol = solve_ivp(rhs, (A0, A_end), [0.0], dense_output=True,
                    rtol=1e-10, atol=1e-13, max_step=1e-4)
    return sol


print("Integrating f_0(A) from f_0(2/3) = 1 ...")
sol_f0 = integrate_f0()
A_grid = np.linspace(2/3 + 1e-4, 1 - 1e-3, 2000)
lnf0 = sol_f0.sol(A_grid).flatten()
f0 = np.exp(lnf0)
df0_dA = f0 * dlnf_dA_fn(A_grid)
Z_AA = f0 * Z_over_f_fn(A_grid)

# Single-scalar G_ghost (for reference)
G_ghost_single = 3 * df0_dA**2 + 2 * Z_AA * f0

print()
print("G70 single-scalar baseline:")
print(f"  G_ghost < 0 fraction (ghost region): "
      f"{100 * np.sum(G_ghost_single < 0) / len(G_ghost_single):.2f}%")
print()


# Ansatz testing
def det_M_2scalar(g_of_A, Z_PsiPsi_of_A):
    """For f(A, Psi) = f_0(A) + g(A) Psi, Psi_bg = 0:
       M = [[Z_AA + (3 f0'^2)/(2 f0),    (3 f0' g)/(2 f0)],
            [(3 f0' g)/(2 f0),           Z_PsiPsi + (3 g^2)/(2 f0)]]
    Returns det M and trace M arrays.
    """
    f0_arr = f0
    f0p_arr = df0_dA
    Z_AA_arr = Z_AA
    g_arr = g_of_A(A_grid)
    Z_PP_arr = Z_PsiPsi_of_A(A_grid)

    M_AA = Z_AA_arr + 3 * f0p_arr**2 / (2 * f0_arr)
    M_AP = 3 * f0p_arr * g_arr / (2 * f0_arr)
    M_PP = Z_PP_arr + 3 * g_arr**2 / (2 * f0_arr)

    det_M = M_AA * M_PP - M_AP**2
    trace_M = M_AA + M_PP
    return M_AA, M_PP, det_M, trace_M


# Note: M_AA here is the same as single-scalar (3 f_0'^2)/(2 f_0) + Z_AA
# which is (G_ghost / (2 f_0)).  So M_AA < 0 where G_ghost < 0.

# For no-ghost: M positive-definite -> M_AA > 0 AND det M > 0.
# Single-scalar already fails M_AA > 0 in inner shell.

# Try ansatz (a): g(A) = g0 * y^n
def make_ansatz_a(g0, n):
    def g(A):
        return g0 * (3*A - 2)**n
    return g


def make_const(value):
    def f(A):
        return value * np.ones_like(A)
    return f


print("=" * 80)
print("Ansatz (a):  g(A) = g0 * y^n     (pair-structure motivated)")
print("=" * 80)
print()
print(f"{'g0':>8}{'n':>4}{'Z_PsiPsi':>12}{'M_AA>0':>10}{'det M>0':>10}{'trace>0':>10}{'verdict':>14}")
print("-" * 78)

best_a = None
for g0 in [0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0]:
    for n in [0, 1, 2, 3, 4]:
        for Z_PP_val in [0.1, 1.0, 10.0, 100.0]:
            g_fn = make_ansatz_a(g0, n)
            Z_PP_fn = make_const(Z_PP_val)
            M_AA, M_PP, det_M, trace_M = det_M_2scalar(g_fn, Z_PP_fn)

            ok_AA = np.all(M_AA > 0)
            ok_det = np.all(det_M > 0)
            ok_tr = np.all(trace_M > 0)
            verdict = "GHOST-FREE" if (ok_AA and ok_det and ok_tr) else "ghost"

            # Only print interesting ones (changes or near-success)
            if ok_AA or ok_det or (g0 in [1.0, 50.0] and n in [0, 2]):
                marker = "***" if verdict == "GHOST-FREE" else ""
                print(f"{g0:>8.2f}{n:>4}{Z_PP_val:>12.2f}"
                      f"{str(ok_AA):>10}{str(ok_det):>10}{str(ok_tr):>10}"
                      f"{verdict:>14}{marker:>5}")

            if verdict == "GHOST-FREE" and best_a is None:
                best_a = (g0, n, Z_PP_val)

print()
if best_a:
    g0, n, Z_PP = best_a
    print(f"FIRST GHOST-FREE config: g(A) = {g0} * y^{n}, Z_PsiPsi = {Z_PP}")
else:
    print("No ghost-free config found in ansatz (a) parameter sweep.")
    # Diagnose: M_AA is single-scalar -- need to lift M_AA itself
    # M_AA = Z_AA + 3 f0'^2/(2 f0) = G_ghost / (2 f0)
    # which is < 0 in inner shell.  Can't be fixed by Psi alone.
    print("Note: M_AA in 2-scalar with no kinetic mixing is identical to single-scalar")
    print("G_ghost / (2 f_0).  Cannot be lifted by adding Psi-only piece.")
    print("Need kinetic mixing Z_APsi != 0  OR  Psi-dependent Z_AA(A, Psi).")
print()

# Try ansatz (b) with KINETIC MIXING Z_APsi(A)
# M_AA gets modified by mixing through "kinetic + R-correction" interplay.
# Full 2-field calculation needed.  Simplified:
#   K = [[Z_AA, Z_APsi], [Z_APsi, Z_PsiPsi]]
#   M = K + (3/(2f)) D D^T
#   M_AA = Z_AA + (3 f0'^2)/(2 f0)
#   M_APsi = Z_APsi + (3 f0' g)/(2 f0)
#   M_PsiPsi = Z_PsiPsi + (3 g^2)/(2 f0)
#
# M_AA unchanged by mixing -- still single-scalar G_ghost / (2 f_0).
# But det M = M_AA M_PsiPsi - M_APsi^2 CAN be made positive even when M_AA < 0
# if M_PsiPsi < 0 too and M_APsi^2 is large.
#
# For both eigenvalues positive (no-ghost): trace > 0 AND det > 0.
# M_AA < 0 means trace = M_AA + M_PsiPsi requires M_PsiPsi > -M_AA = |M_AA|.
# But then M_PsiPsi > 0, so M_AA M_PsiPsi < 0, and det = M_AA M_PsiPsi - M_APsi^2 < 0
# unless ... well, with M_AA M_PsiPsi < 0, det M = negative - |something| < 0.
#
# Conclusion: in scalar-tensor with f(A, Psi), the M_AA channel cannot be
# rescued from single-scalar G_ghost.  The eigenvalue along A-direction
# remains the ghost mode regardless of Psi enrichment.

print("=" * 80)
print("STRUCTURAL CONCLUSION")
print("=" * 80)
print()
print("In 2-scalar scalar-tensor with f(A, Psi), the M_AA component of the")
print("effective kinetic matrix is identical to single-scalar:")
print("   M_AA = Z_AA + (3 f_0'^2) / (2 f_0) = G_ghost / (2 f_0)")
print()
print("Since G_ghost < 0 in 86% of the shell (G71), M_AA < 0 there.")
print()
print("M positive-definite requires M_AA > 0 AND det M > 0.")
print("M_AA < 0 ALONE rules out positive-definite M -- regardless of how Psi enters.")
print()
print("Therefore: SIMPLE 2-scalar Brans-Dicke-type enrichment via Psi cannot")
print("rescue the ghost.  The kinetic-mixing trick doesn't apply to scalar-tensor")
print("with f-coupling: the matrix correction (3/2f) D D^T is rank-1 outer")
print("product, and rank-1 corrections cannot rescue an indefinite K.")
print()
print("What DOES work in principle:")
print("  - DHOST (degenerate higher-order): the (box A)^2 type terms modify")
print("    M structurally, not just by rank-1 correction.")
print("  - Vector-tensor with explicit timelike direction.")
print("  - Disformal coupling (g_mu_nu -> g_mu_nu + B(A) d_mu A d_nu A) in")
print("    the matter sector.  Changes the kinetic structure.")
print("  - Lagrange multiplier / mimetic-type constraint that removes the")
print("    ghost as a non-dynamical mode (cuscuton-like).")
print()
print("Or: accept the substance-ontology reformulation as the framework's")
print("primary closure of Open Problem #6, and treat the Lagrangian-side")
print("ghost as a continuum-field artifact not present in the substance theory.")
