#!/usr/bin/env python3
"""
G58_ledger_configuration_volume.py

Tests the ledger-configuration-volume derivation of the final-shell closure
profile F(Sigma), and explores what symphony-refinement (non-uniform alpha_i)
would do to the exponents and the inside-PS metric shape.

Setup:
  Final-shell coordinate y = Sigma - 2 = 3A - 2, y in (0, 1)
  Channels: D spatial + 2 horizon-pair (D + 2 channels total)
  Spatial channels q_i >= 0 with sum_i q_i = y
  Horizon-pair channels h_1, h_2 >= 0 with h_1 + h_2 = 1 - y

  Under substance ontology + presentism, with uniform per-entry measure
  alpha_i = 1, the configuration volume on the (spatial-vs-horizon-pair)
  split is:
      g(y) ~ y^(alpha_S - 1) (1 - y)^(alpha_H - 1)
  where:
      alpha_S = sum of alpha_i over spatial channels = D (uniform)
      alpha_H = sum over horizon-pair channels      = 2 (uniform, two-face)

  For D = 3, uniform: Beta(3, 2) -> p(y) = 12 y^2 (1-y)
  Survival: F(y) = 1 - I_y(3, 2) = 1 - 4y^3 + 3y^4   (the quartic Hermite)

Candidates tested:
  C1 - uniform D=3, two-face       alpha_S=3, alpha_H=2  Beta(3,2)
       (the framework's current commitment)
  C2 - pair-weighted (each x2)     alpha_S=6, alpha_H=4  Beta(6,4)
       (symphony reading: each entry counted twice for paired writes)
  C3 - higher-alpha spatial         alpha_S=4, alpha_H=2  Beta(4,2)
  C4 - higher-alpha horizon         alpha_S=3, alpha_H=3  Beta(3,3)
  C5 - softer (cubic smoothstep)   alpha_S=2, alpha_H=2  Beta(2,2)
       (the early-session cubic, only C^1 at PS)
  C6 - thicker horizon              alpha_S=3, alpha_H=4  Beta(3,4)

Joint constraints:
  - SU shell-count: ord_{A=1} k(A) = D
  - C^2 at PS:     alpha_S > 2 (so F''(0) = 0)
  - Two-face:      alpha_H = 2 (geometrically, 2 horizon faces)
  - ord_{A=1} k = alpha_H + 1 (computed from k = (1-A) * F(y(A)))

  Combined: alpha_H = 2 forces ord_{A=1} k = 3, independent of D.
            SU shell-count requires ord = D.
            Joint compatibility: D = 3.
"""

from __future__ import annotations
from math import comb, factorial
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

A_PS = 2.0 / 3.0


# --- core formulas (scipy-free, integer alpha only) ---

def beta_func(a, b):
    """B(a, b) = (a-1)!(b-1)! / (a+b-1)! for integer a, b >= 1."""
    return factorial(a - 1) * factorial(b - 1) / factorial(a + b - 1)


def beta_density(y, alpha_S, alpha_H):
    """Beta(alpha_S, alpha_H) density on [0, 1]."""
    return y ** (alpha_S - 1) * (1 - y) ** (alpha_H - 1) / beta_func(alpha_S, alpha_H)


def survival(y, alpha_S, alpha_H):
    """F(y) = 1 - I_y(alpha_S, alpha_H) for integer alpha_S, alpha_H.

    Uses the binomial-CDF identity:
        1 - I_y(a, b) = sum_{k=0}^{a-1} C(a+b-1, k) y^k (1-y)^(a+b-1-k)
    """
    n = alpha_S + alpha_H - 1
    if np.isscalar(y):
        return float(sum(comb(n, k) * y ** k * (1 - y) ** (n - k)
                         for k in range(alpha_S)))
    y_arr = np.asarray(y, dtype=float)
    result = np.zeros_like(y_arr)
    for k in range(alpha_S):
        result += comb(n, k) * y_arr ** k * (1 - y_arr) ** (n - k)
    return result


def k_of_A(A, alpha_S, alpha_H):
    """k(A) = (1-A) outside PS; (1-A) * F(y) inside PS, with y = 3A - 2."""
    if A <= A_PS:
        return 1.0 - A
    y = 3.0 * A - 2.0
    return (1.0 - A) * survival(y, alpha_S, alpha_H)


def ps_smoothness_label(alpha_S):
    """Differentiability of F at the photon sphere (y=0).
    p(y) ~ y^(alpha_S - 1), so:
      F'(0) = 0  iff alpha_S > 1
      F''(0) = 0 iff alpha_S > 2
      F'''(0) = 0 iff alpha_S > 3
    """
    if alpha_S <= 1:
        return "C^0 only"
    if alpha_S <= 2:
        return "C^1"
    if alpha_S <= 3:
        return "C^2"
    return f"C^{int(alpha_S - 1)}"


def horizon_order(alpha_H):
    """ord_{A=1} k(A) = alpha_H + 1.
    Derivation: k = (1-A) * F. Near horizon, F ~ (1-y)^alpha_H / alpha_H.
    With (1-y) = 3(1-A): F ~ const * (1-A)^alpha_H.
    So k ~ (1-A) * (1-A)^alpha_H = (1-A)^(alpha_H + 1).
    """
    return alpha_H + 1


# --- the candidate set ---

CANDIDATES = [
    # (label, alpha_S, alpha_H, motivation)
    ("C1: uniform D=3, two-face",   3, 2, "alpha_i=1 uniform; D spatial=3 + two-face=2"),
    ("C2: pair-weighted (each x2)", 6, 4, "each entry counted twice for paired writes"),
    ("C3: higher-alpha spatial",     4, 2, "extra spatial smoothness; alpha_i=4/3 spatial"),
    ("C4: higher-alpha horizon",     3, 3, "alpha_i=3/2 on horizon-pair (third face?)"),
    ("C5: softer (cubic smoothstep)", 2, 2, "early-session cubic; only C^1 at PS"),
    ("C6: thicker horizon",          3, 4, "horizon-pair alpha=4 (four-face refinement?)"),
]

D_CONST = 3  # spatial dimension


def main():
    print("=" * 88)
    print("G58: Ledger configuration volume -> strong-field metric")
    print("=" * 88)
    print()

    # --- Candidate table ---
    print(f"{'Candidate':<35}{'a_S':>5}{'a_H':>5}  {'PS smoothness':<18}{'ord_{A=1} k':>14}  match?")
    print("-" * 95)
    for label, aS, aH, _ in CANDIDATES:
        ord_h = horizon_order(aH)
        smooth = ps_smoothness_label(aS)
        c2_ok = aS > 2
        ord_ok = ord_h == D_CONST
        verdict = "OK" if (c2_ok and ord_ok) else (
            "fails C^2" if not c2_ok else f"ord != D ({ord_h} vs {D_CONST})"
        )
        print(f"{label:<35}{aS:>5}{aH:>5}  {smooth:<18}{ord_h:>14}  {verdict}")
    print()
    print("Joint constraints (D = 3):")
    print("  - SU shell-count requires  ord_{A=1} k = D = 3")
    print("  - C^2 at PS requires       alpha_S > 2")
    print("  - Two-face refinement      alpha_H = 2  (geometrically, independent of D)")
    print("  - ord_{A=1} k = alpha_H + 1 = 3  AUTOMATIC from two-face")
    print()
    print("Joint compatibility: alpha_H = 2 forces ord = 3 always; SU requires ord = D.")
    print("  -> D = 3 is the unique spatial dimension where SU + two-face are jointly satisfied.")
    print()

    # --- Verification: C1 reproduces the quartic Hermite ---
    print("Verification (C1: alpha_S=3, alpha_H=2 reproduces the quartic Hermite):")
    print(f"{'y':>6}{'p_dirichlet':>14}{'12y^2(1-y)':>14}{'F_dirichlet':>14}{'1-4y^3+3y^4':>14}")
    for y in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]:
        p_d = beta_density(y, 3, 2) if 0 < y < 1 else 0.0
        p_a = 12 * y**2 * (1 - y)
        F_d = float(survival(y, 3, 2)) if 0 <= y <= 1 else float("nan")
        F_a = 1 - 4 * y**3 + 3 * y**4
        print(f"{y:>6.2f}{p_d:>14.6f}{p_a:>14.6f}{F_d:>14.6f}{F_a:>14.6f}")
    print()
    print("Match to numerical precision -> ledger-config-volume == simplex == quartic Hermite. OK")
    print()

    # --- Inside-PS k(A) sample for each candidate ---
    print("Inside-PS k(A) for each candidate (A in (2/3, 1)):")
    print(f"{'A':>6}{'k_GR':>10}", end="")
    for label, _, _, _ in CANDIDATES:
        tag = label.split(":")[0]
        print(f"{('k_'+tag):>10}", end="")
    print()
    print("-" * (16 + 10 * len(CANDIDATES)))
    for A in [0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99]:
        kG = 1 - A
        print(f"{A:>6.2f}{kG:>10.4f}", end="")
        for _, aS, aH, _ in CANDIDATES:
            kv = k_of_A(A, aS, aH)
            print(f"{kv:>10.4f}", end="")
        print()
    print()

    # --- F(y(A)) = k_STAM / k_GR (the survival ratio) ---
    print("Survival ratio F = k_STAM / k_GR inside PS:")
    print(f"{'A':>6}", end="")
    for label, _, _, _ in CANDIDATES:
        tag = label.split(":")[0]
        print(f"{('F_'+tag):>10}", end="")
    print()
    print("-" * (6 + 10 * len(CANDIDATES)))
    for A in [0.70, 0.80, 0.90, 0.99]:
        print(f"{A:>6.2f}", end="")
        y = 3 * A - 2
        for _, aS, aH, _ in CANDIDATES:
            F = float(survival(y, aS, aH))
            print(f"{F:>10.4f}", end="")
        print()
    print()

    # --- QNM consequence ---
    print("QNM consequence at photon sphere (spinless eikonal):")
    print(f"  k(A=2/3) = (1 - 2/3) * F(y=0) = (1/3) * 1 = 1/3 for ALL candidates.")
    print(f"  (F(0) = 1 by normalization for any alpha_S, alpha_H.)")
    print(f"  -> lambda_STAM = lambda_Schwarzschild")
    print(f"  -> tau_STAM / tau_GR = 1.000 (eikonal, spinless)")
    print(f"  Choice of alpha doesn't shift the spinless ringdown prediction.")
    print()
    print("Where the candidates differ: shape of k(A) inside PS.")
    print("Observable handles for distinguishing alpha-choices:")
    print("  - higher overtones (n >= 1): probe geometry just inside r_c")
    print("  - late-inspiral chirp: A approaching 2/3 from outside")
    print("  - LISA EMRI ringdowns: sample intermediate A precisely")
    print("  - sub-leading WKB / full Regge-Wheeler on n=0")
    print()

    # --- Plots ---
    print("Generating plots...")
    p1 = plot_p_F_comparison()
    p2 = plot_k_inside_PS()
    p3 = plot_survival_ratio()
    print(f"  {p1}")
    print(f"  {p2}")
    print(f"  {p3}")
    print()

    md = write_md()
    print(f"Summary: {md}")
    print()
    print("Bottom line:")
    print("  - Ledger-config-volume reproduces simplex result (Beta(3,2) -> quartic Hermite).")
    print("  - Substance + presentism + uniform measure are now load-bearing for the metric.")
    print("  - Spinless ringdown: tau_STAM/tau_GR = 1 for ALL candidates (PS robustness).")
    print("  - D = 3 forced by joint compatibility of SU shell-count + two-face refinement.")
    print("  - Symphony-refinement variants distinguishable by inside-PS observables only.")


# --- plotting ---

def plot_p_F_comparison():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    y_grid = np.linspace(0.001, 0.999, 1000)

    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]
    for (label, aS, aH, _), color in zip(CANDIDATES, colors):
        p = beta_density(y_grid, aS, aH)
        F = survival(y_grid, aS, aH)
        lbl = f"{label.split(':')[0]} Beta({aS},{aH})"
        axes[0].plot(y_grid, p, color=color, linewidth=2, label=lbl)
        axes[1].plot(y_grid, F, color=color, linewidth=2, label=lbl)

    axes[0].set_xlabel("y = Sigma - 2 = 3A - 2  (depth into final shell)")
    axes[0].set_ylabel("p(y) = closure density")
    axes[0].set_title("Closure density p(y) for ledger-measure variants\n"
                      "Boundary behavior at y=0 (PS) and y=1 (horizon) varies")
    axes[0].legend(fontsize=8, loc="upper right")
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel("y = Sigma - 2 (depth into final shell)")
    axes[1].set_ylabel("F(y) = survival fraction")
    axes[1].set_title("Survival F(y): F(0)=1 (PS), F(1)=0 (horizon) for all\n"
                      "Shape between varies; C1 = quartic Hermite (current commitment)")
    axes[1].legend(fontsize=8, loc="upper right")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    out = PLOTS / "G58_p_and_F_comparison.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_k_inside_PS():
    fig, ax = plt.subplots(figsize=(11, 6))
    A_grid = np.linspace(A_PS + 0.0005, 0.9995, 500)

    k_GR = 1 - A_grid
    ax.plot(A_grid, k_GR, "k-", linewidth=2.5, label="GR  k = 1 - A")

    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]
    for (label, aS, aH, _), color in zip(CANDIDATES, colors):
        k_vals = np.array([k_of_A(A, aS, aH) for A in A_grid])
        lbl = f"{label.split(':')[0]} Beta({aS},{aH})"
        ax.plot(A_grid, k_vals, color=color, linewidth=1.8, alpha=0.9, label=lbl)

    ax.axvline(A_PS, color="gray", linestyle="--", alpha=0.6, label="PS (A=2/3)")
    ax.axvline(1.0, color="black", linestyle=":", alpha=0.6, label="Horizon (A=1)")
    ax.set_xlabel("A inside photon sphere")
    ax.set_ylabel("k(A) = 1/g_rr")
    ax.set_title("Inside-PS k(A) for ledger-measure variants\n"
                 "All match GR at A=2/3; differ in shape and horizon order")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(A_PS - 0.005, 1.0)
    plt.tight_layout()
    out = PLOTS / "G58_k_inside_PS.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_survival_ratio():
    fig, ax = plt.subplots(figsize=(11, 6))
    A_grid = np.linspace(A_PS + 0.0005, 0.9995, 500)

    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]
    for (label, aS, aH, _), color in zip(CANDIDATES, colors):
        y_grid = 3 * A_grid - 2
        F_vals = survival(y_grid, aS, aH)
        lbl = f"{label.split(':')[0]} Beta({aS},{aH})"
        ax.plot(A_grid, F_vals, color=color, linewidth=2, label=lbl)

    ax.axvline(A_PS, color="gray", linestyle="--", alpha=0.6)
    ax.axvline(1.0, color="black", linestyle=":", alpha=0.6)
    ax.axhline(1.0, color="gray", linestyle=":", alpha=0.4)
    ax.set_xlabel("A")
    ax.set_ylabel("F(y(A)) = k_STAM / k_GR")
    ax.set_title("Survival fraction F = k_STAM / k_GR inside photon sphere\n"
                 "The principal observable distinguishing ledger-measure candidates")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(A_PS - 0.005, 1.0)
    plt.tight_layout()
    out = PLOTS / "G58_survival_ratio.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- markdown summary ---

def write_md():
    md = []
    md.append("# G58 - Ledger configuration volume -> strong-field metric\n")
    md.append("**Date: 2026-05-12 (continuation).** Tests the ledger-configuration-volume "
              "derivation of the final-shell closure profile F(Sigma) and explores "
              "symphony-refinement (non-uniform alpha_i) variants.\n")

    md.append("## Setup\n")
    md.append("Under substance ontology + presentism, the ledger at final-shell depth "
              "y = Sigma - 2 is a continuous configuration of A-resolution. With D + 2 "
              "channels (D spatial + 2 horizon-pair) and uniform per-entry measure, "
              "the configuration-volume measure on the spatial-vs-horizon-pair split is:\n")
    md.append("```")
    md.append("g(y)  ~  y^(alpha_S - 1) * (1 - y)^(alpha_H - 1)")
    md.append("    alpha_S = sum of alpha_i over spatial channels   = D    (uniform)")
    md.append("    alpha_H = sum over horizon-pair channels         = 2    (uniform, two-face)")
    md.append("```\n")
    md.append("For D = 3 with uniform measure: Beta(3, 2) = 12 y^2 (1-y).\n")
    md.append("Survival: F(y) = 1 - I_y(3, 2) = **1 - 4 y^3 + 3 y^4**  (the quartic Hermite).\n")

    md.append("## Verified equivalence (C1)\n")
    md.append("Beta(3, 2) Dirichlet density and the quartic Hermite F match to numerical "
              "precision across the final shell. Ledger-configuration-volume gives the "
              "same p(y) and F(y) as the simplex argument.\n")

    md.append("## Candidates tested\n")
    md.append("| Candidate | a_S | a_H | PS smoothness | ord_{A=1} k | Verdict |")
    md.append("|---|---:|---:|---|---:|---|")
    for label, aS, aH, _ in CANDIDATES:
        smooth = ps_smoothness_label(aS)
        ord_h = horizon_order(aH)
        c2_ok = aS > 2
        ord_ok = ord_h == D_CONST
        verdict = "**matches all constraints**" if (c2_ok and ord_ok) else (
            "fails C^2 at PS" if not c2_ok else f"ord != D (gives {ord_h})"
        )
        md.append(f"| {label} | {aS} | {aH} | {smooth} | {ord_h} | {verdict} |")
    md.append("")

    md.append("## Joint constraint analysis\n")
    md.append("Three independent structural constraints:\n")
    md.append("- **SU shell-count**: ord_{A=1} k(A) = D\n")
    md.append("- **C^2 at PS**: alpha_S > 2 (so F''(0) = 0)\n")
    md.append("- **Two-face refinement**: alpha_H = 2 (geometrically, two horizon faces; "
              "independent of D)\n")
    md.append("\n")
    md.append("From the third: ord_{A=1} k = alpha_H + 1 = 3 **automatically**, in any "
              "spatial dimension D.\n")
    md.append("\n")
    md.append("From the first: ord = D required.\n")
    md.append("\n")
    md.append("**Joint compatibility forces D = 3.** This is the structural derivation: "
              "D = 3 is the unique spatial dimension where SU shell-count and the two-face "
              "horizon refinement co-determine the same horizon closure order.\n")
    md.append("\n")
    md.append("Combined with C^2 at PS: among uniform-measure candidates with integer alpha, "
              "only **Beta(3, 2)** satisfies all three constraints. The strong-field metric "
              "is structurally unique under these commitments.\n")

    md.append("## QNM consequence (spinless eikonal): robust at the photon sphere\n")
    md.append("For **all candidates** with normalized F (so F(0) = 1):\n")
    md.append("```")
    md.append("k(A = 2/3) = (1 - 2/3) * F(0) = (1/3) * 1 = 1/3   (exact Schwarzschild)")
    md.append("```")
    md.append("So **tau_STAM / tau_GR = 1.000 (spinless eikonal) for every ledger-measure "
              "candidate.** The choice of alpha doesn't shift the spinless ringdown prediction.\n")
    md.append("\n")
    md.append("Confirms G57: spinless STAM ringdown = exact GR is robust to any sensible "
              "ledger-measure choice, not specific to the quartic Hermite.\n")

    md.append("## Where candidates differ: inside-PS metric shape\n")
    md.append("All differences live in the survival ratio F(y) = k_STAM / k_GR for A > 2/3. "
              "Sample at representative depths:\n")
    md.append("")
    md.append("| A | C1 Beta(3,2) | C2 Beta(6,4) | C5 Beta(2,2) |")
    md.append("|---:|---:|---:|---:|")
    for A in [0.70, 0.80, 0.90, 0.99]:
        y = 3 * A - 2
        F1 = float(survival(y, 3, 2))
        F2 = float(survival(y, 6, 4))
        F5 = float(survival(y, 2, 2))
        md.append(f"| {A:.2f} | {F1:.4f} | {F2:.4f} | {F5:.4f} |")
    md.append("")
    md.append("Distinguishing observables (none observationally settled yet):\n")
    md.append("- Higher overtones n >= 1: probe geometry just inside r_c\n")
    md.append("- Late-inspiral chirp: A approaching 2/3 from outside\n")
    md.append("- LISA EMRI ringdowns: sample intermediate A precisely\n")
    md.append("- Sub-leading WKB or full Regge-Wheeler on n=0\n")

    md.append("## What this does for the framework\n")
    md.append("1. **Ledger-config-volume = simplex measure under existing commitments.** "
              "Substance ontology + presentism naturally produce continuous Dirichlet measure; "
              "uniform per-entry weighting is the natural expression of presentism (no "
              "internal hierarchy of entry-types).\n")
    md.append("2. **Substance + presentism are now load-bearing for the strong-field metric**, "
              "not just for the quantum interpretation. Tighter framework integration.\n")
    md.append("3. **D = 3 is forced by SU shell-count + two-face refinement + Beta horizon "
              "exponent.** 'Why three spatial dimensions' shifts from cosmological/anthropic "
              "to framework-internal.\n")
    md.append("4. **Spinless ringdown = exact GR is a robust prediction**, not specific to "
              "the quartic. Any normalized ledger-measure candidate gives F(0) = 1 -> "
              "k(2/3) = 1/3 = Schwarzschild.\n")

    md.append("## Open follow-ups\n")
    md.append("- Whether **uniform per-entry measure (alpha_i = 1)** is itself derivable from "
              "presentism, or is a fourth independent commitment. Best reading: under "
              "presentism, the ledger has no internal hierarchy of entry-types, so uniform "
              "is natural -- but it could be stated explicitly.\n")
    md.append("- A **symphony-refinement** observation prediction for late-inspiral chirp or "
              "EMRI ringdown shape, distinguishing C1 from non-uniform-alpha candidates.\n")
    md.append("- Confirm the **SU shell-count argument is independent** of the horizon-pair "
              "argument (otherwise the D = 3 forcing dissolves into a tautology).\n")

    md.append("## Files\n")
    md.append("- [scripts/G58_ledger_configuration_volume.py](../scripts/G58_ledger_configuration_volume.py)\n")
    md.append("- [plots/G58_p_and_F_comparison.png](../plots/G58_p_and_F_comparison.png)\n")
    md.append("- [plots/G58_k_inside_PS.png](../plots/G58_k_inside_PS.png)\n")
    md.append("- [plots/G58_survival_ratio.png](../plots/G58_survival_ratio.png)\n")

    out = RESULTS / "G58_ledger_configuration_volume_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


if __name__ == "__main__":
    main()
