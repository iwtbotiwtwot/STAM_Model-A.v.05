#!/usr/bin/env python3
"""
F1_penrose_singularity_evasion.py

FATALITY 1: Bold-STAM evasion of the Penrose 1965 singularity theorem.
Author: Sean Brady / STAM Model-A continuation

The opponent's attack:
    "Penrose's theorem (1965) is a mathematical theorem of general relativity.
     It proves that a closed trapped surface plus the null energy condition
     plus a non-compact Cauchy surface implies an incomplete null geodesic
     (a singularity). Bold STAM denies that singularities form. Therefore
     bold STAM contradicts a proven theorem."

Bold STAM's defense:
    Bold STAM does NOT deny the theorem. The theorem's PREMISE is unmet in
    bold STAM: no closed *strictly* trapped surface ever forms during
    gravitational collapse. The theorem applies to spacetimes that contain
    such a surface; bold-STAM spacetimes do not.

The mechanism:
    A closed trapped surface requires both null expansions θ_+ and θ_- to
    be STRICTLY negative on a closed 2-surface. For a sphere of areal
    radius r in a spherically symmetric metric ds^2 = -h(A) c^2 dt^2 +
    dr^2/k(A) + r^2 dΩ^2, the outgoing null expansion is

        θ_+(r) = (2/r) × sqrt(h(A) k(A))      [in Schwarzschild-time
                                                 parameterization]

    For Schwarzschild (h = k = 1-A):
        θ_+(r) = (2/r) × (1-A) = 0  at A=1 (marginally trapped)
        θ_+(r) < 0                  for A>1 (strictly trapped, interior)

    For bold STAM (h = 1-A, k = (1-A)(1-A^2)^2):
        θ_+(r) = (2/r) × (1-A)(1-A^2) = 0  at A=1
        θ_+(r) > 0                          for all A < 1
        and there is no A > 1 region (universe ends at A=1)

    Therefore: in bold STAM, no closed surface has θ_+ < 0 strictly.
    The "trapped region" of Schwarzschild simply does not exist in bold
    STAM, because the manifold itself does not extend to A > 1.

Consequence:
    Bold STAM has null geodesic completeness at the boundary. Affine
    parameter to reach A=1 along any radial null geodesic is INFINITE
    (logarithmic divergence), versus FINITE in Schwarzschild. So both
    the premise AND the conclusion of the theorem are absent — fully
    consistent with the theorem, just outside its domain of application.

What this script does:
    1. Computes θ_+(A) for Schwarzschild and bold STAM, plots both.
    2. Computes the affine parameter to reach the horizon along a null
       geodesic in both metrics, showing the divergence in bold STAM.
    3. Writes a clean markdown summary stating the evasion as positive
       content (not a denial) — bold STAM is a different geometric class
       of spacetime than the one Penrose's theorem applies to.
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


# --- Metric coefficients ---

def h_schwarzschild(A: np.ndarray) -> np.ndarray:
    return 1.0 - A


def k_schwarzschild(A: np.ndarray) -> np.ndarray:
    return 1.0 - A


def h_bold(A: np.ndarray) -> np.ndarray:
    return 1.0 - A


def k_bold(A: np.ndarray) -> np.ndarray:
    return (1.0 - A) * (1.0 - A ** 2) ** 2


# --- Outgoing null expansion (Schwarzschild-time parameterization) ---

def theta_plus(A: np.ndarray, h_func, k_func, r_over_Rs: np.ndarray = None) -> np.ndarray:
    """θ_+(A) for outgoing null normal on sphere of areal radius r = Rs/A.

    θ_+ = (2/r) × sqrt(h(A) k(A))

    Returns θ_+ in units of 1/Rs (so factor 1/Rs is implicit in the curve).
    """
    if r_over_Rs is None:
        r_over_Rs = 1.0 / A
    return (2.0 / r_over_Rs) * np.sqrt(h_func(A) * k_func(A))


# --- Affine parameter for radial null geodesic ---

def affine_parameter_to_horizon(A_start: float, h_func, k_func, n_points: int = 4000,
                                A_max: float = 0.99999) -> float:
    """Compute affine parameter λ along ingoing radial null geodesic from A_start to A=A_max.

    For metric ds^2 = -h c^2 dt^2 + dr^2/k:
        dr/dλ = -E sqrt(k/h)  for ingoing (E = energy along geodesic, set =1)
    With dr = -Rs/A^2 dA:
        dλ = Rs/(A^2 sqrt(k/h)) dA / E
    Returns λ in units of Rs/c (so output × Rs/c gives physical seconds).
    """
    A_grid = np.linspace(A_start, A_max, n_points)
    sqrt_kh = np.sqrt(k_func(A_grid) / h_func(A_grid))
    integrand = 1.0 / (A_grid ** 2 * sqrt_kh)
    return float(np.trapezoid(integrand, A_grid))


# --- Tables ---

def build_theta_table() -> pd.DataFrame:
    A_values = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999, 0.9999])
    rows = []
    for A in A_values:
        theta_S = theta_plus(np.array([A]), h_schwarzschild, k_schwarzschild)[0]
        theta_B = theta_plus(np.array([A]), h_bold, k_bold)[0]
        rows.append({
            "A": A,
            "r_over_Rs": 1.0 / A,
            "theta_plus_Schwarzschild_per_Rs": theta_S,
            "theta_plus_bold_STAM_per_Rs": theta_B,
            "ratio_bold_over_Schw": theta_B / theta_S if theta_S > 0 else float("nan"),
        })
    return pd.DataFrame(rows)


def build_affine_table() -> pd.DataFrame:
    A_starts = np.array([0.5, 0.9, 0.99, 0.999, 0.9999])
    rows = []
    for A_s in A_starts:
        lam_S = affine_parameter_to_horizon(A_s, h_schwarzschild, k_schwarzschild,
                                            A_max=min(0.99999, A_s + 0.1))
        # Schwarzschild has constant dr/dλ for ingoing null (with E=1)
        # so affine parameter from r_start to Rs is just (r_start - Rs)/E = (1/A_s - 1) Rs
        lam_S_analytical = (1.0 / A_s) - 1.0  # in units of Rs/c

        lam_B = affine_parameter_to_horizon(A_s, h_bold, k_bold, A_max=0.999999)
        rows.append({
            "A_start": A_s,
            "r_start_over_Rs": 1.0 / A_s,
            "lambda_to_horizon_Schwarzschild_Rs_over_c": lam_S_analytical,
            "lambda_to_horizon_bold_STAM_Rs_over_c": lam_B,
            "lambda_ratio_bold_over_Schw": lam_B / lam_S_analytical if lam_S_analytical > 0 else float("nan"),
        })
    return pd.DataFrame(rows)


# --- Plots ---

def plot_theta_comparison() -> Path:
    A_grid = np.linspace(0.001, 0.999, 1000)
    theta_S = theta_plus(A_grid, h_schwarzschild, k_schwarzschild)
    theta_B = theta_plus(A_grid, h_bold, k_bold)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(A_grid, theta_S, color="tab:red", linewidth=2.5, label="Schwarzschild θ_+ (= (2/r)(1-A))")
    ax.plot(A_grid, theta_B, color="tab:blue", linewidth=2.5, label="Bold STAM θ_+ (= (2/r)(1-A)(1-A²))")
    ax.axhline(0, color="black", linewidth=0.6, linestyle="--", alpha=0.5)
    ax.set_xlabel("A")
    ax.set_ylabel("θ_+ (in units of 1/Rs)")
    ax.set_title("Outgoing null expansion θ_+: bold STAM vanishes more sharply at A=1, "
                 "no interior trapped region")
    ax.grid(True, linewidth=0.3)
    ax.legend()
    ax.set_xlim(0, 1)

    out = PLOTS / "F1_theta_comparison.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_affine_parameter_divergence() -> Path:
    """Show affine parameter to reach horizon as A_start → 1.

    Schwarzschild: λ → 0 (already at horizon).
    Bold STAM: λ stays finite from any fixed start; but the integrand diverges so for
    fixed starting position farther out, λ grows logarithmically as A_max → 1.
    """
    A_max_values = 1.0 - np.logspace(-1, -7, 25)  # how close we get to A=1
    A_start = 0.5
    lam_S = []
    lam_B = []
    for A_max in A_max_values:
        # for Schwarzschild, λ from A_start to A_max is (1/A_start - 1/A_max) in Rs/c
        lam_S_value = 1.0 / A_start - 1.0 / A_max
        lam_S.append(lam_S_value)

        lam_B_value = affine_parameter_to_horizon(A_start, h_bold, k_bold, A_max=A_max)
        lam_B.append(lam_B_value)
    lam_S = np.array(lam_S)
    lam_B = np.array(lam_B)
    one_minus_A_max = 1.0 - A_max_values

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.loglog(one_minus_A_max, np.abs(lam_S), "o-", color="tab:red", label="Schwarzschild — finite limit")
    ax.loglog(one_minus_A_max, lam_B, "o-", color="tab:blue", label="Bold STAM — diverges logarithmically")
    ax.invert_xaxis()
    ax.set_xlabel("1 − A_max  (how close the geodesic is sampled to the boundary)")
    ax.set_ylabel("Affine parameter λ to A_max  (in units of Rs/c)")
    ax.set_title("Bold STAM: null geodesic affine parameter to A=1 diverges (geodesic complete)")
    ax.grid(True, which="both", linewidth=0.3)
    ax.legend()

    out = PLOTS / "F1_affine_parameter_divergence.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_trapped_region_topology() -> Path:
    """Schematic showing where trapped surfaces exist in each theory."""
    A_grid = np.linspace(0.001, 1.499, 1000)

    fig, ax = plt.subplots(figsize=(10, 5))
    # Schwarzschild interior is trapped (A > 1)
    ax.fill_between(A_grid, -5, 5, where=(A_grid > 1.0), color="tab:red", alpha=0.25,
                    label="Schwarzschild trapped region (A > 1)")
    # Schwarzschild horizon
    ax.axvline(1.0, color="tab:red", linewidth=2, linestyle="--", alpha=0.7,
               label="Schwarzschild horizon (marginally trapped)")
    # Bold STAM doesn't extend past A=1 — no manifold there
    ax.fill_between(A_grid, -5, 5, where=(A_grid > 1.0), color="white", alpha=1, zorder=2)
    ax.fill_between(A_grid, -5, 5, where=(A_grid > 1.0), hatch="X", color="none", edgecolor="black",
                    label="Bold STAM: no manifold here (universe ends at A=1)")
    ax.axvline(1.0, color="tab:blue", linewidth=3, alpha=0.7,
               label="Bold STAM A=1 boundary (universe edge)")

    ax.text(0.5, 0, "Bold STAM domain (A < 1)\nNo strict trapped surface exists here\n(θ_+ > 0 always)",
            ha="center", va="center", fontsize=10, bbox=dict(facecolor="lightblue", alpha=0.7))
    ax.text(1.25, 0, "No bold-STAM\nmanifold here", ha="center", va="center", fontsize=10)

    ax.set_xlim(0, 1.5)
    ax.set_ylim(-1, 1)
    ax.set_yticks([])
    ax.set_xlabel("A")
    ax.set_title("Where strict trapped surfaces exist: Schwarzschild interior vs bold STAM (nowhere)")
    ax.legend(loc="upper right", fontsize=9)

    out = PLOTS / "F1_trapped_region_topology.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- Markdown summary ---

def write_markdown(theta_df: pd.DataFrame, affine_df: pd.DataFrame, plot_paths: list[Path]) -> Path:
    md = []
    md.append("# STAM Model-A F1: Penrose-Hawking Singularity Theorem Evasion\n")
    md.append("## The attack\n")
    md.append(
        "Penrose's 1965 singularity theorem states: in any spacetime satisfying\n"
        "1. The null energy condition: `T_μν k^μ k^ν ≥ 0` for all null `k^μ`,\n"
        "2. A non-compact Cauchy surface,\n"
        "3. **The existence of a closed (compact, no-boundary) *strictly* trapped surface**, "
        "meaning a closed 2-surface where both null expansions θ_+ and θ_- are strictly "
        "negative,\n"
        "the spacetime is null-geodesically incomplete — i.e. it contains a singularity.\n"
        "\n"
        "An opponent will say: bold STAM denies that black holes have singularities. But "
        "Penrose proved a singularity must form when matter collapses inside a horizon. "
        "Therefore bold STAM contradicts a proven theorem of mathematical physics.\n"
    )
    md.append("## Bold STAM's defense\n")
    md.append(
        "Bold STAM does NOT deny the theorem. **The theorem's premise (3) is unmet** in "
        "bold STAM. The theorem applies to spacetimes that contain a strict closed trapped "
        "surface. Bold-STAM spacetimes do not.\n"
        "\n"
        "More precisely: bold STAM has a fundamentally different geometric class than a "
        "GR Schwarzschild black hole. In Schwarzschild, the interior region `r < Rs` (i.e. "
        "`A > 1`) contains spheres on which both null expansions are strictly negative — "
        "these are the trapped surfaces that the theorem requires. In bold STAM, **the "
        "manifold does not extend to A > 1 at all** (the universe ends at the boundary). "
        "There is therefore no interior region in which strict trapped surfaces could form.\n"
    )
    md.append("## The math: outgoing null expansion θ_+\n")
    md.append(
        "For a sphere of areal radius `r = Rs/A` in a spherically symmetric metric "
        "`ds² = -h(A)c²dt² + dr²/k(A) + r²dΩ²`, the outgoing null expansion (in "
        "Schwarzschild-time parameterization) is:\n"
        "```text\n"
        "θ_+(r) = (2/r) × sqrt(h(A) k(A))\n"
        "```\n"
        "\n"
        "**Schwarzschild** (h = k = 1-A):\n"
        "```text\n"
        "θ_+ = (2/r)(1-A)\n"
        "    > 0  for A < 1   (untrapped, exterior)\n"
        "    = 0  at A = 1    (marginally trapped, horizon)\n"
        "    < 0  for A > 1   (strictly trapped, interior)\n"
        "```\n"
        "\n"
        "**Bold STAM** (h = 1-A, k = (1-A)(1-A²)²):\n"
        "```text\n"
        "θ_+ = (2/r)(1-A)(1-A²)\n"
        "    > 0  for A < 1   (untrapped)\n"
        "    = 0  at A = 1    (marginally trapped, with double zero — even more degenerate)\n"
        "    [no manifold for A > 1]\n"
        "```\n"
    )
    md.append(theta_df.to_markdown(index=False, floatfmt=".6g"))
    md.append("\n")
    md.append(
        "**Verdict on premise (3):** in bold STAM, no closed surface has both null expansions "
        "strictly negative. The strict-trapped condition is never satisfied. The theorem's "
        "premise is unmet, so its conclusion (incomplete null geodesic) is not forced.\n"
    )
    md.append("## Geodesic completeness — automatic consequence\n")
    md.append(
        "A direct check: compute the affine parameter `λ` along an ingoing radial null "
        "geodesic from a starting position to the horizon. In Schwarzschild this is finite "
        "(geodesic incomplete; reaches r=Rs in bounded affine parameter, then continues "
        "into singularity at r=0). In bold STAM the integral diverges logarithmically at "
        "A=1 (geodesic asymptotes to boundary, never reaches it in finite affine parameter).\n"
    )
    md.append(affine_df.to_markdown(index=False, floatfmt=".6g"))
    md.append("\n")
    md.append(
        "Bold-STAM null geodesics are complete in both directions: they approach but never "
        "reach the boundary. This is consistent with — and a consequence of — the absence of "
        "trapped surfaces in the strict sense.\n"
    )
    md.append("## Where bold STAM and Schwarzschild differ as geometric objects\n")
    md.append(
        "- **Schwarzschild manifold**: extends to A → ∞ (r → 0). Trapped region exists for "
        "A > 1. Singularity at A = ∞ (r = 0). Penrose's theorem applies.\n"
        "- **Bold-STAM manifold**: A < 1 only. No interior. No trapped region. No singularity. "
        "Penrose's theorem premise unmet — theorem does not apply.\n"
        "\n"
        "Bold STAM is not denying Penrose. Bold STAM is constructing a *different geometric "
        "object* that the theorem does not address. The theorem says 'IF a strict trapped "
        "surface exists, THEN incomplete geodesics exist.' Bold STAM says 'no such surface "
        "exists in our manifold.' Both statements are simultaneously true.\n"
    )
    md.append("## What this gives bold STAM\n")
    md.append(
        "- A clean answer to 'how do you avoid singularities Penrose proved must form?': by "
        "having a manifold class outside the theorem's domain. The theorem applies to "
        "spacetimes with strictly trapped surfaces; bold-STAM spacetimes don't have them.\n"
        "- A *positive* statement of the bubble picture: gravitational collapse in bold STAM "
        "asymptotically approaches a marginally-trapped boundary (θ_+ = 0) without ever "
        "forming a strictly trapped region (θ_+ < 0). The bubble surface IS the marginally "
        "trapped boundary; matter accumulates on it, never enters a (nonexistent) interior.\n"
        "- Geodesic completeness as a derived consequence, not a postulate.\n"
    )
    md.append("## What this fatality establishes for the framework\n")
    md.append(
        "**Fatality 1: SURVIVED.** Bold STAM is consistent with Penrose's 1965 theorem. The "
        "theorem applies to spacetimes containing strictly trapped surfaces; bold-STAM "
        "spacetimes don't contain them. There is no contradiction with the theorem; bold "
        "STAM is in a different geometric class.\n"
        "\n"
        "Caveat (honest): we have not yet shown that *gravitational collapse dynamics* in "
        "bold STAM smoothly produce a manifold of this class without ever transiently "
        "forming a strictly trapped region. That requires dynamical-A propagation work "
        "(open). The static-state argument here shows the equilibrium configuration evades "
        "the theorem; the collapse-process argument is a separate piece of work.\n"
    )
    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "F1_penrose_singularity_evasion_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    theta_df = build_theta_table()
    affine_df = build_affine_table()

    theta_df.to_csv(RESULTS / "F1_theta_table.csv", index=False)
    affine_df.to_csv(RESULTS / "F1_affine_table.csv", index=False)

    plot_paths = [
        plot_theta_comparison(),
        plot_affine_parameter_divergence(),
        plot_trapped_region_topology(),
    ]

    summary = write_markdown(theta_df, affine_df, plot_paths)

    print("STAM Model-A F1: Penrose-Hawking Singularity Evasion")
    print("=" * 60)
    print()
    print("Outgoing null expansion theta_+ comparison:")
    print(theta_df.to_string(index=False))
    print()
    print("Affine parameter to horizon along null geodesic:")
    print(affine_df.to_string(index=False))
    print()
    print("Verdict: bold STAM has no strict closed trapped surface.")
    print("Penrose's theorem premise is unmet. Conclusion is not forced.")
    print("Null geodesic affine parameter to A=1 diverges logarithmically.")
    print("Geodesic completeness automatic.")
    print()
    print("Files written:")
    print(f"- {RESULTS / 'F1_theta_table.csv'}")
    print(f"- {RESULTS / 'F1_affine_table.csv'}")
    print(f"- {summary}")
    for p in plot_paths:
        print(f"- {p}")


if __name__ == "__main__":
    main()
