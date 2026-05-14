#!/usr/bin/env python3
"""
G71_ghost_freedom_check.py

Step 3 of the Lagrangian-for-A track (Open Problem #6).

Takes the smallest covariant action from G70:

    S = (1/16 pi G) integral sqrt(-g) [ f(A) R - Z(A) (grad A)^2 - 2 V(A) ] d^4 x

and linearizes around the committed-metric / A = 2M/r background.  The
standard scalar-tensor result is that the propagating scalar mode has
effective kinetic norm proportional to

    G_ghost(A) = 3 [ f'(A) ]^2  +  2 Z(A) f(A)

(equivalently, the effective Brans-Dicke parameter omega_BD = Z f / (f')^2
must satisfy 2 omega_BD + 3 > 0).  Ghost-freedom requires G_ghost > 0
throughout the support of the modification (the final shell A in (2/3, 1)).
Graviton positivity additionally requires f(A) > 0.

From G70, using f' = f * d(ln f)/dA and Z = f * (Z/f), the overall scale of
f cancels and:

    G_ghost(A) / f(A)^2  =  3 [d(ln f)/dA]^2  +  2 (Z/f)
                         = (4 y^2 / [A^2 F(y)^2]) * S(A)

with the *signed* polynomial:

    S(A)  =  3 y^4 (45 A^2 - 33 A - 13)^2  +  P_8(A)

where y = 3A - 2 and P_8(A) is the degree-8 polynomial determined in G70.
The sign of G_ghost(A) is the sign of S(A) (everything else in the
prefactor is non-negative inside the shell).

This script:
  1. Constructs S(A) symbolically.
  2. Evaluates S(A) on a fine grid over (2/3, 1).
  3. Locates zero crossings of S(A).
  4. Also numerically integrates f(A) from f(2/3) = 1 and verifies
     G_ghost(A) directly, cross-checking the structural shortcut.
  5. Reports verdict on Open Problem #6 (perturbative-level ghost
     freedom of the smallest scalar-tensor embedding of the committed
     metric).
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


# ---------------------------------------------------------------------------
# Symbolic ghost polynomial S(A)
# ---------------------------------------------------------------------------

def build_S_of_A():
    A = sp.symbols('A', positive=True)
    y = 3 * A - 2

    # P_8 from G70
    P8 = (393660 * A**8
          - 2318220 * A**7
          + 5730669 * A**6
          - 7684146 * A**5
          + 5982660 * A**4
          - 2636280 * A**3
          + 551025 * A**2
          - 8580 * A
          - 10790)

    S = 3 * y**4 * (45 * A**2 - 33 * A - 13)**2 + P8
    S_expanded = sp.expand(S)
    return A, S_expanded, P8


# ---------------------------------------------------------------------------
# Numerical evaluation
# ---------------------------------------------------------------------------

def evaluate_ghost_polynomial():
    A_sym, S_expr, P8_expr = build_S_of_A()

    print("=" * 80)
    print("STEP 1: Construct S(A) symbolically")
    print("=" * 80)
    print()
    print("S(A) = 3 y^4 (45 A^2 - 33 A - 13)^2 + P_8(A),  y = 3 A - 2")
    print()
    print("Expanded S(A):")
    print(f"  {S_expr}")
    print()

    # Check degree and leading term
    poly = sp.Poly(S_expr, A_sym)
    print(f"  degree(S) = {poly.degree()}")
    print(f"  leading coefficient = {poly.LC()}")
    print(f"  constant term (S(0)) = {poly.eval(0)}")
    print()

    # Sample values at endpoints and interior
    print("=" * 80)
    print("STEP 2: Sample S(A) at endpoints and selected interior points")
    print("=" * 80)
    print()
    print(f"{'A':>10}{'y = 3A-2':>14}{'S(A)':>20}{'sign':>10}")
    print("-" * 54)
    sample_A = [
        sp.Rational(2, 3),       # PS
        sp.Rational(0.67),       # just inside PS
        sp.Rational(0.7),
        sp.Rational(0.75),
        sp.Rational(0.8),
        sp.Rational(0.85),
        sp.Rational(0.9),
        sp.Rational(0.95),
        sp.Rational(0.99),
        sp.Rational(1),           # horizon
    ]
    for a in sample_A:
        y_val = 3 * a - 2
        S_val = float(sp.N(S_expr.subs(A_sym, a)))
        sign = '+' if S_val > 0 else ('-' if S_val < 0 else '0')
        print(f"{float(a):>10.4f}{float(y_val):>14.4f}{S_val:>20.4f}{sign:>10}")
    print()

    # Find real roots inside (2/3, 1)
    print("=" * 80)
    print("STEP 3: Find zero crossings of S(A) in (2/3, 1)")
    print("=" * 80)
    print()

    # Coefficients of S(A) as a polynomial in A
    coeffs = [float(c) for c in poly.all_coeffs()]
    roots = np.roots(coeffs)
    real_roots = [r.real for r in roots if abs(r.imag) < 1e-10]
    interior = sorted([r for r in real_roots if 2/3 - 1e-12 < r < 1 + 1e-12])
    print(f"Real roots of S(A) in the closed interval [2/3, 1]:")
    if not interior:
        print("  None.")
    else:
        for r in interior:
            print(f"  A = {r:.10f}")
    print()
    all_roots_str = ", ".join(f"{r:.6f}" for r in sorted(real_roots))
    print(f"All real roots: [{all_roots_str}]")
    print()

    # Verdict on sign behavior
    print("=" * 80)
    print("STEP 4: Sign verdict across (2/3, 1)")
    print("=" * 80)
    print()
    A_grid = np.linspace(2/3 + 1e-6, 1 - 1e-6, 100000)
    S_grid = np.array([float(sp.N(S_expr.subs(A_sym, sp.Float(a)))) for a in A_grid[::1000]])
    # Use numpy for fast evaluation instead
    S_fn = sp.lambdify(A_sym, S_expr, 'numpy')
    S_grid = S_fn(A_grid)

    neg_count = np.sum(S_grid < 0)
    pos_count = np.sum(S_grid > 0)
    total = len(S_grid)

    print(f"  Sample size: {total} points uniformly across (2/3 + epsilon, 1 - epsilon)")
    print(f"  S(A) < 0 : {neg_count} points  ({100*neg_count/total:.2f}%)")
    print(f"  S(A) > 0 : {pos_count} points  ({100*pos_count/total:.2f}%)")
    print()

    if neg_count == 0:
        print("  G_ghost = (positive prefactor) * S(A) > 0 throughout the shell.")
        print("  -> NO GHOST.  Smallest covariant embedding is ghost-free.")
    elif pos_count == 0:
        print("  G_ghost < 0 throughout the shell.  Smallest covariant embedding")
        print("  carries a ghost everywhere inside PS.")
    else:
        # Find the crossover
        # locate sign change
        sign = np.sign(S_grid)
        change_idx = np.where(np.diff(sign) != 0)[0]
        print(f"  Sign changes at approximately:")
        for idx in change_idx:
            A_cross = (A_grid[idx] + A_grid[idx + 1]) / 2
            print(f"     A = {A_cross:.6f}  (y = {3*A_cross - 2:.6f})")
        print()
        # Sign at PS+, mid-shell, horizon-
        print(f"  S(A=2/3+) ~ {S_grid[0]:.4f}    (just inside PS)")
        print(f"  S(A=5/6)  ~ {S_fn(np.array([5/6]))[0]:.4f}    (mid-shell)")
        print(f"  S(A=1-)   ~ {S_grid[-1]:.4f}   (just below horizon)")
        print()
        if S_grid[0] < 0 and S_grid[-1] > 0:
            print("  G_ghost < 0 in the inner part of the shell (ghost present),")
            print("  G_ghost > 0 in the outer part (healthy).")
        elif S_grid[0] > 0 and S_grid[-1] < 0:
            print("  G_ghost > 0 near PS, ghost in the deeper shell.")
        else:
            print("  Sign pattern is complicated -- multiple crossings.")

    return A_sym, S_expr, S_fn, A_grid, S_grid


# ---------------------------------------------------------------------------
# Cross-check: numerically integrate f(A) and verify G_ghost
# ---------------------------------------------------------------------------

def numerical_crosscheck():
    print("=" * 80)
    print("STEP 5: Numerical cross-check")
    print("=" * 80)
    print()
    print("Integrate f(A) from f(2/3) = 1, compute G_ghost(A) = 3(f')^2 + 2 Z f")
    print("directly, and verify sign(G_ghost) == sign(S).")
    print()

    A_sym = sp.symbols('A', positive=True)
    y = 3 * A_sym - 2
    F_metric = 1 - 5 * y**4 + 4 * y**5
    Q3 = 108 * A_sym**3 - 189 * A_sym**2 + 114 * A_sym - 23
    P8 = (393660 * A_sym**8 - 2318220 * A_sym**7 + 5730669 * A_sym**6
          - 7684146 * A_sym**5 + 5982660 * A_sym**4 - 2636280 * A_sym**3
          + 551025 * A_sym**2 - 8580 * A_sym - 10790)

    dlnf_dA = -2 * y**3 * (45 * A_sym**2 - 33 * A_sym - 13) / (A_sym * F_metric)
    Z_over_f = 2 * y**2 * P8 / (81 * A_sym**2 * (1 - A_sym)**4 * Q3**2)

    dlnf_dA_fn = sp.lambdify(A_sym, dlnf_dA, 'numpy')
    Z_over_f_fn = sp.lambdify(A_sym, Z_over_f, 'numpy')

    def rhs(A, lnf):
        return dlnf_dA_fn(A)

    # Integrate ln f from PS+ to near horizon
    A0 = 2/3 + 1e-6
    A_end = 1 - 1e-4
    sol = solve_ivp(rhs, (A0, A_end), [0.0], dense_output=True,
                    rtol=1e-10, atol=1e-13, max_step=1e-4)

    A_grid = np.linspace(A0, A_end, 5000)
    lnf = sol.sol(A_grid).flatten()
    f = np.exp(lnf)
    fprime = f * dlnf_dA_fn(A_grid)
    Z = f * Z_over_f_fn(A_grid)
    G_ghost = 3 * fprime**2 + 2 * Z * f

    f_min = f.min()
    f_max = f.max()
    print(f"  f(A=2/3+)  = {f[0]:.6f}")
    print(f"  f(A near horizon)  = {f[-1]:.6e}")
    print(f"  f always positive: {bool(np.all(f > 0))} (graviton sign OK)")
    print()
    print(f"  Z(A=2/3+)  = {Z[0]:.6f}     (should be ~0 at PS)")
    print(f"  Z(A=0.85)  = {Z[np.argmin(np.abs(A_grid - 0.85))]:.6e}")
    print(f"  Z(A near horizon) = {Z[-1]:.6e}")
    print()
    print(f"  3 (f')^2 (A=2/3+)  = {3*fprime[0]**2:.6e}")
    print(f"  3 (f')^2 (A=0.85)  = {3*fprime[np.argmin(np.abs(A_grid - 0.85))]**2:.6e}")
    print(f"  3 (f')^2 (A near horizon) = {3*fprime[-1]**2:.6e}")
    print()
    print(f"  G_ghost (A=2/3+)  = {G_ghost[0]:.6e}")
    print(f"  G_ghost (A=0.75)  = {G_ghost[np.argmin(np.abs(A_grid - 0.75))]:.6e}")
    print(f"  G_ghost (A=0.85)  = {G_ghost[np.argmin(np.abs(A_grid - 0.85))]:.6e}")
    print(f"  G_ghost (A=0.95)  = {G_ghost[np.argmin(np.abs(A_grid - 0.95))]:.6e}")
    print(f"  G_ghost (A near horizon) = {G_ghost[-1]:.6e}")
    print()

    neg = np.sum(G_ghost < 0)
    pos = np.sum(G_ghost > 0)
    print(f"  G_ghost < 0: {neg}/{len(G_ghost)}  ({100*neg/len(G_ghost):.2f}%)")
    print(f"  G_ghost > 0: {pos}/{len(G_ghost)}  ({100*pos/len(G_ghost):.2f}%)")
    print()

    # Find crossings
    sign_changes = np.where(np.diff(np.sign(G_ghost)) != 0)[0]
    if len(sign_changes) > 0:
        print(f"  G_ghost zero crossings inside the shell:")
        for idx in sign_changes:
            A_c = (A_grid[idx] + A_grid[idx + 1]) / 2
            print(f"     A ~ {A_c:.6f}")
        print()

    return A_grid, f, fprime, Z, G_ghost


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def make_plots(A_grid, f, fprime, Z, G_ghost, S_fn):
    A_sym = sp.symbols('A', positive=True)
    S_grid = S_fn(A_grid)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    ax = axes[0, 0]
    ax.plot(A_grid, f, 'tab:blue', linewidth=2)
    ax.axvline(2/3, color='gray', linestyle='--', alpha=0.5, label='PS (A=2/3)')
    ax.axvline(1, color='black', linestyle=':', alpha=0.5, label='Horizon (A=1)')
    ax.set_xlabel('A')
    ax.set_ylabel('f(A)')
    ax.set_title('Non-minimal coupling f(A)\n(graviton kinetic norm)')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(A_grid, Z, 'tab:orange', linewidth=2)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.axvline(2/3, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(1, color='black', linestyle=':', alpha=0.5)
    ax.set_xlabel('A')
    ax.set_ylabel('Z(A)')
    ax.set_title('Bare kinetic coefficient Z(A)\n(NOT the no-ghost coefficient)')
    ax.set_yscale('symlog', linthresh=1e-3)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(A_grid, G_ghost, 'tab:red', linewidth=2)
    ax.axhline(0, color='black', linewidth=1, linestyle='--',
               label='G_ghost = 0  (ghost threshold)')
    ax.axvline(2/3, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(1, color='black', linestyle=':', alpha=0.5)
    ax.set_xlabel('A')
    ax.set_ylabel(r'$G_{\rm ghost} = 3(f\')^2 + 2 Z f$')
    ax.set_title('No-ghost coefficient G_ghost(A)\nNegative = ghost; positive = healthy')
    ax.set_yscale('symlog', linthresh=1e-3)
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(A_grid, S_grid, 'tab:purple', linewidth=2,
            label=r'$S(A) = 3 y^4 (45A^2 - 33A - 13)^2 + P_8(A)$')
    ax.axhline(0, color='black', linewidth=1, linestyle='--', label='S = 0')
    ax.axvline(2/3, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(1, color='black', linestyle=':', alpha=0.5)
    ax.set_xlabel('A')
    ax.set_ylabel('S(A)  (sign-determining polynomial)')
    ax.set_title('Sign of G_ghost is sign of S(A) inside the shell')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out = PLOTS / "G71_ghost_freedom_check.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def write_summary(verdict_lines, S_expr, A_root_locations):
    md = []
    md.append("# G71 - Ghost-freedom check on smallest covariant action\n")
    md.append("**Date: 2026-05-13.**  Step 3 of the Lagrangian-for-A track "
              "(Open Problem #6).  Linearizes the smallest covariant action "
              "from G70 around the committed-metric / A = 2M/r background and "
              "checks the sign of the scalar-mode kinetic operator.\n")
    md.append("## Setup\n")
    md.append("Smallest covariant action (G70):\n")
    md.append("```\n"
              "S = (1 / 16 pi G) integral sqrt(-g) [ f(A) R "
              "- Z(A) (grad A)^2 - 2 V(A) ] d^4 x\n"
              "```\n")
    md.append("Standard scalar-tensor result: linearizing around the static "
              "spherical background, the propagating scalar mode has "
              "effective Brans-Dicke parameter\n")
    md.append("```\n"
              "omega_BD(A) = Z(A) f(A) / (f'(A))^2\n"
              "```\n")
    md.append("Ghost-freedom: 2 omega_BD + 3 > 0, equivalently\n")
    md.append("```\n"
              "G_ghost(A) = 3 (f'(A))^2 + 2 Z(A) f(A)  >  0\n"
              "```\n")
    md.append("Graviton positivity: f(A) > 0.\n")
    md.append("## Sign reduction\n")
    md.append("Using f' = f * d(ln f)/dA and Z = f * (Z/f), the overall scale "
              "of f cancels:\n")
    md.append("```\n"
              "G_ghost / f^2 = (4 y^2 / [A^2 F(y)^2]) * S(A)\n"
              "```\n")
    md.append("with\n")
    md.append("```\n"
              "S(A) = 3 y^4 (45 A^2 - 33 A - 13)^2 + P_8(A)\n"
              "```\n")
    md.append(f"The prefactor 4 y^2 / [A^2 F(y)^2] is non-negative inside the shell, so "
              "**sign(G_ghost) = sign(S(A))**.\n")
    md.append("## Result\n")
    for line in verdict_lines:
        md.append(line + "\n")
    md.append("\n")
    md.append("Zero crossings of S(A) inside (2/3, 1):\n")
    if not A_root_locations:
        md.append("- None.\n")
    else:
        for r in A_root_locations:
            md.append(f"- A = {r:.10f}  (y = {3*r - 2:.10f})\n")
    md.append("\n")
    md.append("## Files\n")
    md.append("- [scripts/G71_ghost_freedom_check.py](../scripts/G71_ghost_freedom_check.py)\n")
    md.append("- [plots/G71_ghost_freedom_check.png](../plots/G71_ghost_freedom_check.png)\n")

    out = RESULTS / "G71_ghost_freedom_check_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    print(f"Summary written: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("G71: Ghost-freedom check on smallest covariant action")
    print("=" * 80)
    print()

    A_sym, S_expr, S_fn, A_grid_sym, S_grid_sym = evaluate_ghost_polynomial()
    A_grid, f, fprime, Z, G_ghost = numerical_crosscheck()

    plot_path = make_plots(A_grid, f, fprime, Z, G_ghost, S_fn)
    print(f"  Plot: {plot_path}")
    print()

    # Verdict lines for summary
    neg = np.sum(G_ghost < 0)
    pos = np.sum(G_ghost > 0)
    verdict = []
    if neg == 0:
        verdict.append("**G_ghost > 0 throughout the final shell.**  Smallest covariant "
                       "scalar-tensor embedding of the committed metric is **ghost-free** "
                       "at the perturbative kinetic-operator level.")
    elif pos == 0:
        verdict.append("**G_ghost < 0 throughout the final shell.**  Smallest covariant "
                       "scalar-tensor embedding carries a ghost everywhere inside PS.")
    else:
        sign_changes = np.where(np.diff(np.sign(G_ghost)) != 0)[0]
        cross_A = [(A_grid[i] + A_grid[i + 1]) / 2 for i in sign_changes]
        verdict.append(
            f"**Sign of G_ghost is NOT uniform across the shell.** "
            f"Zero crossing(s) at A in [{', '.join(f'{c:.6f}' for c in cross_A)}].")
        if G_ghost[0] < 0:
            verdict.append(
                f"Inner-shell region (A just inside PS, up to A = {cross_A[0]:.6f}) "
                f"has G_ghost < 0 -- ghost mode in the smallest scalar-tensor embedding.")
            verdict.append(
                f"Outer-shell region (A > {cross_A[-1]:.6f}, up through horizon) "
                f"has G_ghost > 0 -- healthy.")
        else:
            verdict.append(
                f"Inner-shell region has G_ghost > 0 (healthy); deeper region has ghost.")

    # find root locations (in shell)
    poly = sp.Poly(S_expr, A_sym)
    coeffs = [float(c) for c in poly.all_coeffs()]
    roots = np.roots(coeffs)
    real_roots = [r.real for r in roots if abs(r.imag) < 1e-10]
    interior = sorted([r for r in real_roots if 2/3 - 1e-12 < r < 1 + 1e-12])

    write_summary(verdict, S_expr, interior)

    print("=" * 80)
    print("VERDICT")
    print("=" * 80)
    for line in verdict:
        print(line)
    print()


if __name__ == "__main__":
    main()
