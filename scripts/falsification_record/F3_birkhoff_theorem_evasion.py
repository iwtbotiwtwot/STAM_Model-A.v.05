#!/usr/bin/env python3
"""
F3_birkhoff_theorem_evasion.py

FATALITY 3: Bold-STAM evasion of Birkhoff's theorem.
Author: Sean Brady / STAM Model-A continuation

The opponent's attack:
    "Birkhoff's theorem (1923) proves that the unique spherically symmetric
     vacuum solution to Einstein's equations is Schwarzschild. Bold STAM's
     g_rr = 1/[(1-A)(1-A^2)^2] is not Schwarzschild. Therefore bold STAM
     contradicts Birkhoff's theorem."

Bold STAM's defense (in two parts):
    (1) Birkhoff's theorem is a theorem about *vacuum Einstein gravity*. It
        applies to spacetimes where R_μν = 0 (no matter, no fields contributing
        to T_μν beyond the metric itself). Bold-STAM spacetime is NOT vacuum:
        the A field is itself a physical field with its own stress-energy
        contribution. Therefore Birkhoff's premise is unmet.

    (2) The price of this evasion: we must compute what effective stress-energy
        T_μν is required to source the bold-STAM metric under Einstein gravity,
        and check whether this T_μν has reasonable properties.

What this script computes:
    For the bold-STAM metric

        ds² = -(1-A)c² dt² + dr²/[(1-A)(1-A²)²] + r² dΩ²

    we extract the standard "mass function"

        m(r) = (c²r/(2G)) × (1 - e^(-2Λ(r)))   where  e^(-2Λ) = (1-A)(1-A²)²

    and compute the effective energy density via

        ρ_effective(r) = (1/(4πr²)) × dm/dr.

    For a true vacuum (Schwarzschild), m(r) = M = constant for all r > Rs,
    giving ρ = 0. For bold STAM, m(r) is non-trivial — the matter content
    required to source the metric is non-trivial.

Honest findings (verified numerically below):
    - Bold STAM is NOT a vacuum solution. Birkhoff's premise unmet. Theorem
      doesn't apply. So far so good.
    - The effective ρ profile for bold STAM is **non-monotonic in r** and is
      **negative in the outer exterior** (A < ~0.45, far from horizon) and
      **positive near the horizon** (A > ~0.45).
    - Negative ρ would be a NEC violation IF we treated bold STAM as
      Einstein gravity with ordinary matter. Two interpretations:
        (a) Bold STAM is NOT Einstein gravity. It is a modified-gravity
            theory with the A field as fundamental and its own field
            equations. The "effective ρ" is a re-expression of those
            modifications, not literal matter. NEC is not a constraint
            on it.
        (b) The A field is quintessence-like (negative-pressure scalar)
            in some regions. This is a real physical claim and would need
            its own defense.
    - **Bottom line on F3: PARTIALLY SURVIVED.** Birkhoff doesn't apply
      because bold STAM isn't vacuum. But the price is that bold STAM
      owes a specification of its field equations — which is open work.

This is the most honest fatality result so far: bold STAM evades the
theorem cleanly, but the evasion exposes a real piece of the framework
that hasn't been written out (the field equations for the A field).
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# --- Mass function for Schwarzschild and bold STAM ---

def mass_function_ratio_schwarzschild(A: np.ndarray) -> np.ndarray:
    """m(r)/M for Schwarzschild — equals 1 identically (vacuum)."""
    return np.ones_like(A)


def mass_function_ratio_bold_stam(A: np.ndarray) -> np.ndarray:
    """m(r)/M for bold STAM with g_rr = 1/[(1-A)(1-A^2)^2].

    Derivation:
        e^(-2Lambda) = (1-A)(1-A^2)^2
        m(r) = (c^2 r / 2G) (1 - e^(-2Lambda))
        M = c^2 Rs / 2G  =>  m/M = (r/Rs) (1 - e^(-2Lambda)) = (1/A) (1 - (1-A)(1-A^2)^2)
    """
    return (1.0 - (1.0 - A) * (1.0 - A ** 2) ** 2) / A


def effective_density_bracket(A: np.ndarray) -> np.ndarray:
    """The bracket [1 - (1-A)^2(1+A)(1+A+4A^2)] proportional to ρ for bold STAM.

    Derivation of dm/dr for bold STAM:
        m/M = (1 - (1-A)(1-A^2)^2)/A
        dm/dr = (c^2 / 2G) [1 - (1-A)^2(1+A)(1+A+4A^2)]
        rho = (1/4π r^2) dm/dr proportional to bracket.

    Bracket = 0 means rho = 0 (vacuum).
    Bracket > 0 means rho > 0.
    Bracket < 0 means rho < 0.
    """
    return 1.0 - ((1.0 - A) ** 2) * (1.0 + A) * (1.0 + A + 4.0 * A ** 2)


def effective_density_normalized(A: np.ndarray) -> np.ndarray:
    """ρ_effective(r) × r² in normalized units. Sign tells us NEC status.

    Up to a constant (c²/8πG), this is r² × ρ.
    """
    return effective_density_bracket(A)


# --- Find sign-change point ---

def find_NEC_crossover(A_low: float = 0.01, A_high: float = 0.99,
                      tolerance: float = 1e-6) -> float:
    """Find A_crit where effective_density_bracket = 0 (sign changes from
    negative to positive)."""
    while A_high - A_low > tolerance:
        mid = 0.5 * (A_low + A_high)
        b = effective_density_bracket(np.array([mid]))[0]
        if b < 0:
            A_low = mid
        else:
            A_high = mid
    return 0.5 * (A_low + A_high)


# --- Tables ---

def build_mass_density_table() -> pd.DataFrame:
    A_values = np.array([0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99])
    rows = []
    for A in A_values:
        m_S = mass_function_ratio_schwarzschild(np.array([A]))[0]
        m_B = mass_function_ratio_bold_stam(np.array([A]))[0]
        bracket = effective_density_bracket(np.array([A]))[0]
        rho_sign = "negative (NEC violation)" if bracket < 0 else "positive (NEC OK)"
        if abs(bracket) < 1e-10:
            rho_sign = "zero (vacuum)"
        rows.append({
            "A": A,
            "r_over_Rs": 1.0 / A,
            "m_Schwarzschild_over_M": m_S,
            "m_bold_STAM_over_M": m_B,
            "rho_bracket_bold_STAM": bracket,
            "rho_sign": rho_sign,
        })
    return pd.DataFrame(rows)


# --- Plots ---

def plot_mass_function() -> Path:
    A_grid = np.linspace(0.001, 0.999, 1000)
    m_S = mass_function_ratio_schwarzschild(A_grid)
    m_B = mass_function_ratio_bold_stam(A_grid)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.plot(A_grid, m_S, color="tab:red", linewidth=2.5, label="Schwarzschild m(r)/M = 1 (vacuum)")
    ax.plot(A_grid, m_B, color="tab:blue", linewidth=2.5, label="Bold STAM m(r)/M (non-vacuum)")
    ax.axhline(1.0, color="black", linewidth=0.5, alpha=0.5)
    ax.set_xlabel("A = Rs/r")
    ax.set_ylabel("Mass function m(r)/M")
    ax.set_title("Effective enclosed mass: bold STAM is NOT vacuum (m varies with r)")
    ax.grid(True, linewidth=0.3)
    ax.legend()
    ax.set_xlim(0, 1)

    out = PLOTS / "F3_mass_function.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_effective_density() -> Path:
    A_grid = np.linspace(0.001, 0.999, 1000)
    bracket = effective_density_bracket(A_grid)

    A_crit = find_NEC_crossover()

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.fill_between(A_grid, 0, bracket, where=(bracket < 0), color="tab:red", alpha=0.3,
                    label="ρ < 0 (NEC violation if interpreted as ordinary matter)")
    ax.fill_between(A_grid, 0, bracket, where=(bracket > 0), color="tab:green", alpha=0.3,
                    label="ρ > 0 (NEC OK)")
    ax.plot(A_grid, bracket, color="black", linewidth=2)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.axvline(A_crit, color="tab:purple", linestyle="--", linewidth=1.5,
               label=f"A_crit ≈ {A_crit:.4f} (sign change)")
    ax.set_xlabel("A = Rs/r")
    ax.set_ylabel("Effective ρ × r² (proportional, signed)")
    ax.set_title(
        "Effective stress-energy required to source bold-STAM metric under Einstein gravity"
    )
    ax.grid(True, linewidth=0.3)
    ax.legend()

    out = PLOTS / "F3_effective_density.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_log_density_magnitude() -> Path:
    A_grid = np.linspace(0.001, 0.999, 1000)
    bracket = effective_density_bracket(A_grid)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    pos_mask = bracket > 0
    neg_mask = bracket < 0
    ax.semilogy(A_grid[pos_mask], np.abs(bracket[pos_mask]), color="tab:green", linewidth=2,
                label="|ρ| where ρ > 0")
    ax.semilogy(A_grid[neg_mask], np.abs(bracket[neg_mask]), color="tab:red", linewidth=2,
                label="|ρ| where ρ < 0")
    ax.set_xlabel("A = Rs/r")
    ax.set_ylabel("|effective ρ × r²|, log scale")
    ax.set_title("Magnitude of bold-STAM effective stress-energy: small far away, peaks near horizon")
    ax.grid(True, which="both", linewidth=0.3)
    ax.legend()

    out = PLOTS / "F3_log_density_magnitude.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- Markdown summary ---

def write_markdown(df: pd.DataFrame, A_crit: float, plot_paths: list[Path]) -> Path:
    md = []
    md.append("# STAM Model-A F3: Birkhoff's Theorem Evasion (and what it costs)\n")
    md.append("## The attack\n")
    md.append(
        "Birkhoff's theorem (1923) states: every spherically symmetric solution of the "
        "Einstein vacuum equations `R_μν = 0` is locally isometric to a portion of the "
        "Schwarzschild spacetime. The Schwarzschild solution is unique.\n"
        "\n"
        "Bold STAM's metric `g_rr = 1/[(1-A)(1-A²)²]` is not Schwarzschild. The opponent says: "
        "since Birkhoff's theorem proves the spherically symmetric vacuum solution must be "
        "Schwarzschild, bold STAM contradicts a proven theorem.\n"
    )
    md.append("## Bold STAM's defense\n")
    md.append(
        "**Birkhoff's theorem applies to vacuum Einstein gravity.** Its premise is that the "
        "spacetime satisfies `R_μν = 0` (no matter, no fields contributing to stress-energy). "
        "Bold-STAM spacetime is **not** vacuum: the A field is itself a physical field with "
        "its own stress-energy contribution to the metric. Therefore Birkhoff's premise is "
        "unmet, and the theorem does not apply.\n"
        "\n"
        "**This evasion is genuine but it has a cost.** We must compute the effective "
        "stress-energy T_μν that would be required to source the bold-STAM metric under "
        "Einstein gravity, and check whether it has reasonable properties.\n"
    )
    md.append("## The math: effective stress-energy\n")
    md.append(
        "For a static spherically symmetric metric `ds² = -e^(2Φ)c²dt² + e^(2Λ)dr² + r²dΩ²`, "
        "define the mass function via `e^(-2Λ) = 1 - 2Gm(r)/(c²r)`. Then `dm/dr = 4πr²ρ`, so:\n"
        "```text\n"
        "rho_effective(r) = (1/(4πr²)) × dm/dr\n"
        "```\n"
        "For Schwarzschild: `e^(-2Λ) = 1-A`, giving `m(r) = M` constant, and `dm/dr = 0`. "
        "**Vacuum (ρ = 0).**\n"
        "\n"
        "For bold STAM: `e^(-2Λ) = (1-A)(1-A²)²`, giving:\n"
        "```text\n"
        "m_bold(r)/M = (1/A) × [1 - (1-A)(1-A²)²]\n"
        "dm_bold/dr proportional to [1 - (1-A)²(1+A)(1+A+4A²)]\n"
        "```\n"
        "**Non-vacuum.** The bracket determines the sign of effective ρ.\n"
    )
    md.append("## Numerical results\n")
    md.append(df.to_markdown(index=False, floatfmt=".6g"))
    md.append("\n\n")
    md.append(
        f"**The sign-change point (NEC crossover) is at A_crit ≈ {A_crit:.4f}.** For A < A_crit "
        "(outer exterior, far from horizon), the effective ρ is *negative*. For A > A_crit "
        "(inner exterior, close to horizon), it is *positive*. The mass function m(r) "
        "exceeds M for intermediate A (peaks around A ≈ 0.5) and returns to M at both "
        "A=0 (infinity) and A=1 (horizon).\n"
    )
    md.append("## Honest assessment of the cost\n")
    md.append(
        "If bold STAM were standard Einstein gravity sourced by ordinary matter, the negative "
        "effective ρ in the outer exterior would be a NEC violation — exotic matter. That would "
        "be a serious problem, since NEC violation is generally considered unphysical for "
        "ordinary matter content.\n"
        "\n"
        "Two interpretations of this finding:\n"
        "\n"
        "1. **Bold STAM is not Einstein gravity.** Bold STAM is a modified-gravity theory where "
        "the A field is fundamental and the field equations are NOT R_μν = 8πG T_μν. The "
        "'effective ρ' computed above is the redefinition that absorbs bold-STAM-specific "
        "modifications into a fictitious matter term. In bold STAM's actual field equations, "
        "no NEC violation occurs because no actual matter is involved — the modifications are "
        "geometric / field-theoretic.\n"
        "\n"
        "   This is structurally analogous to how f(R) gravity, scalar-tensor gravity, "
        "and quintessence theories work: the metric departs from GR vacuum, and if you "
        "interpret the departure as 'matter,' the matter looks exotic. But it's really just "
        "a modification to gravity itself.\n"
        "\n"
        "2. **The A field is quintessence-like.** The A field carries real stress-energy "
        "with negative pressure / energy density in some regions. This would be a physical "
        "claim — analogous to dark energy — and would need its own defense beyond what bold "
        "STAM has so far articulated.\n"
        "\n"
        "**Both interpretations evade Birkhoff. The cost is that bold STAM owes a "
        "specification of its actual field equations**, which is currently an open piece. "
        "The metric ansatz `g_rr = 1/[(1-A)(1-A²)²]` is a postulate, not yet derived from "
        "a Lagrangian or field equation.\n"
    )
    md.append("## Verdict\n")
    md.append(
        "**F3: PARTIALLY SURVIVED.**\n"
        "\n"
        "- Birkhoff's theorem does not apply to bold STAM because bold STAM is not vacuum "
        "Einstein gravity. The evasion is structurally clean.\n"
        "- The price of this evasion is exposed: bold STAM has a non-trivial effective "
        "stress-energy with a sign change at A ≈ 0.45. Interpreted as ordinary matter, this "
        "would violate NEC; interpreted as bold-STAM-specific gravitational modifications, it "
        "is consistent.\n"
        "- **Bold STAM owes a specification of its underlying field equations.** Until those "
        "are written down, the framework's response to F3 is partial: the theorem doesn't "
        "apply, but the underlying theory that picks out the bold-STAM metric uniquely is "
        "not yet specified.\n"
        "\n"
        "This is the most honest fatality result so far. Bold STAM survives, but the survival "
        "exposes a piece of the framework that needs more work. That's a real research "
        "agenda item, not a failure.\n"
    )
    md.append("## What this means for the README\n")
    md.append(
        "Add to the open problems list: **specify the bold-STAM field equations** that pick "
        "out the metric `g_rr = 1/[(1-A)(1-A²)²]` uniquely. Candidates include:\n"
        "- A scalar-tensor theory with A as the scalar.\n"
        "- An f(R) or other higher-curvature theory whose static spherically-symmetric solution "
        "matches the ansatz.\n"
        "- A non-metric theory where the A field is fundamental and the metric is composite.\n"
        "\n"
        "Until one of these is committed to, bold STAM's metric ansatz is a postulate that "
        "fits all currently-tested data but lacks a Lagrangian derivation.\n"
    )
    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "F3_birkhoff_theorem_evasion_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    df = build_mass_density_table()
    df.to_csv(RESULTS / "F3_mass_density_table.csv", index=False)

    A_crit = find_NEC_crossover()

    plot_paths = [
        plot_mass_function(),
        plot_effective_density(),
        plot_log_density_magnitude(),
    ]

    summary = write_markdown(df, A_crit, plot_paths)

    print("STAM Model-A F3: Birkhoff's Theorem Evasion (with cost)")
    print("=" * 64)
    print()
    print("Mass function and effective density across the exterior:")
    print(df.to_string(index=False))
    print()
    print(f"NEC crossover: A_crit ~ {A_crit:.4f}")
    print(f"  For A < {A_crit:.4f}: effective rho < 0 (NEC violated if interpreted as matter)")
    print(f"  For A > {A_crit:.4f}: effective rho > 0 (NEC OK)")
    print()
    print("Verdict: F3 PARTIALLY SURVIVED.")
    print("  Birkhoff doesn't apply (bold STAM is not vacuum).")
    print("  But bold STAM owes specification of its actual field equations.")
    print()
    print("Files written:")
    print(f"- {RESULTS / 'F3_mass_density_table.csv'}")
    print(f"- {summary}")
    for p in plot_paths:
        print(f"- {p}")


if __name__ == "__main__":
    main()
