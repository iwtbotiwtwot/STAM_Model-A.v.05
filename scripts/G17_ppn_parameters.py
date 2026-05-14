#!/usr/bin/env python3
"""
G17_ppn_parameters.py

Compute Model-A's PPN parameters (γ, β, and higher-order g_ij coefficients)
from the strong-field metric:

    g_tt = -(1 - A) c²
    g_rr = 1 / [(1 - A)(1 - A²)²]
    g_θθ = r²,  g_φφ = r² sin²θ
    A = 2GM/(c² r) = 2U_Sch (where U_Sch = GM/(c²r) in Schwarzschild coords)

PPN formalism requires the metric in *isotropic* coordinates. We do the
coordinate transformation r → r̃ explicitly (symbolically), expand to
high enough order in u = GM/(c²r̃) to read off both PPN parameters,
and compare to GR Schwarzschild done the same way.

Standard PPN form (isotropic coords, static spherical):
    g_tt = -(1 - 2u + 2β u² + O(u³)) c²
    g_ij = (1 + 2γ u + (3/2 + δ₂) u² + O(u³)) δ_ij

GR Schwarzschild gives γ = β = 1 and δ₂ = 0.

Observational bounds (current):
    Cassini:               |γ − 1| < 2.3 × 10⁻⁵
    Mercury perihelion:    |β − 1| < few × 10⁻⁴
    Lunar laser ranging:   |β − 1| < few × 10⁻⁴
    Higher-order (δ₂):     essentially unconstrained at solar system

Why this is the right test for today: it doesn't touch cosmological
distances. It tests whether Model-A's modified k(A) shows up in
weak-field solar-system observables, and where Model-A first diverges
from GR if at all.
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)


def main() -> None:
    print("G17: PPN parameter calculation for Model-A")
    print("=" * 78)
    print()

    u = sp.Symbol("u", positive=True)  # u = GM/(c² r̃) — isotropic Newtonian potential
    r = sp.Symbol("r", positive=True)  # Schwarzschild radial coordinate
    rt = sp.Symbol("rtilde", positive=True)  # isotropic radial coordinate

    # --- Step 1: write metric coefficients in Schwarzschild coords ---
    # Both theories use A = 2GM/(c²r). Define A symbolically.
    A_sch = sp.Symbol("A")  # symbolic A as expansion variable
    h_GR = 1 - A_sch                          # both theories share g_tt = -(1-A)c²
    k_GR = 1 - A_sch                          # GR Schwarzschild: g_rr = 1/(1-A)
    k_MA = (1 - A_sch) * (1 - A_sch**2)**2    # Model-A: g_rr = 1/[(1-A)(1-A²)²]

    print("Schwarzschild-coordinate metric coefficients (g_rr = 1/k(A)):")
    print(f"  GR:       k_GR(A)      = {sp.simplify(k_GR)}")
    print(f"  Model-A:  k_MA(A)      = {sp.expand(k_MA)}")
    print()

    # Expand 1/k for both, in powers of A
    inv_k_GR = sp.series(1/k_GR, A_sch, 0, 5).removeO()
    inv_k_MA = sp.series(1/k_MA, A_sch, 0, 5).removeO()
    print("Reciprocal expansions in A (= 2u in Schwarzschild coords):")
    print(f"  1/k_GR(A)  = {sp.expand(inv_k_GR)}")
    print(f"  1/k_MA(A)  = {sp.expand(inv_k_MA)}")
    print()

    # In Schwarzschild coords, A = 2 U_Sch. Substitute and rewrite.
    U_sch = sp.Symbol("U_Sch", positive=True)
    inv_k_GR_U = sp.expand(inv_k_GR.subs(A_sch, 2*U_sch))
    inv_k_MA_U = sp.expand(inv_k_MA.subs(A_sch, 2*U_sch))
    print("g_rr in Schwarzschild coordinates expanded in U_Sch = GM/(c²r):")
    print(f"  GR:       g_rr_GR  = {inv_k_GR_U}")
    print(f"  Model-A:  g_rr_MA  = {inv_k_MA_U}")
    print()

    # Compare second-order coefficients in Schwarzschild coords
    # GR coefficient of U_Sch² is 4; Model-A is computed below
    coeff_GR_2 = sp.Poly(inv_k_GR_U, U_sch).coeff_monomial(U_sch**2)
    coeff_MA_2 = sp.Poly(inv_k_MA_U, U_sch).coeff_monomial(U_sch**2)
    print(f"Schwarzschild-coordinate g_rr second-order coefficient (U_Sch²):")
    print(f"  GR:       {coeff_GR_2}")
    print(f"  Model-A:  {coeff_MA_2}")
    print(f"  Ratio Model-A / GR: {sp.simplify(coeff_MA_2 / coeff_GR_2)}")
    print()

    # --- Step 2: solve for the isotropic transformation ---
    # In isotropic coords: G(r̃) r̃² = r², so r = r̃ √G.
    # And: G dr̃² = dr²/k(A)
    # Combined: G' r̃ + 2G(1 - √k) = 0  →  G'/G = 2(√k - 1)/r̃

    print("=" * 78)
    print("Isotropic coordinate transformation r(r̃)")
    print("=" * 78)
    print()
    print("Setting G(r̃) = 1 + α₁ u + α₂ u² + α₃ u³ + ..., expand the radial")
    print("ODE order by order to extract the αᵢ. This determines γ and β.")
    print()

    # Define G as power series with unknown coefficients α₁ α₂ α₃
    a1, a2, a3 = sp.symbols("a1 a2 a3", real=True)
    G = 1 + a1*u + a2*u**2 + a3*u**3
    sqrtG = sp.series(sp.sqrt(G), u, 0, 5).removeO()
    inv_sqrtG = sp.series(1/sp.sqrt(G), u, 0, 5).removeO()
    Gprime_rt = sp.series(-(a1*u + 2*a2*u**2 + 3*a3*u**3), u, 0, 5).removeO()
    # r̃ × dG/dr̃ in u-form: u = GM/(c² r̃), so du/dr̃ = -u/r̃, and dG/dr̃ × r̃ = -(dG/du)·u·r̃/r̃ = -u·dG/du
    rt_Gprime = -u * sp.diff(G, u)
    rt_Gprime = sp.expand(rt_Gprime)
    print(f"r̃ × G' (in terms of u):  {rt_Gprime}")
    print()

    # Now do the same for each theory. A = 2u/√G (since A_sch = 2GM/(c²r) and r = r̃√G).
    A_iso = 2*u * inv_sqrtG
    A_iso = sp.series(A_iso, u, 0, 5).removeO()
    A_iso = sp.expand(A_iso)
    print(f"A in isotropic coords (parametric in α's):  A = {A_iso}")
    print()

    # k(A) for each theory, substituted with A_iso
    k_GR_iso = sp.series(k_GR.subs(A_sch, A_iso), u, 0, 5).removeO()
    k_GR_iso = sp.expand(k_GR_iso)
    k_MA_iso = sp.series(k_MA.subs(A_sch, A_iso), u, 0, 5).removeO()
    k_MA_iso = sp.expand(k_MA_iso)

    sqrt_k_GR = sp.series(sp.sqrt(k_GR_iso), u, 0, 4).removeO()
    sqrt_k_MA = sp.series(sp.sqrt(k_MA_iso), u, 0, 4).removeO()

    # ODE: r̃ G' = 2G(√k - 1)  →  -u dG/du = 2G(√k - 1)
    rhs_GR = sp.expand(sp.series(2*G*(sqrt_k_GR - 1), u, 0, 4).removeO())
    rhs_MA = sp.expand(sp.series(2*G*(sqrt_k_MA - 1), u, 0, 4).removeO())
    lhs = sp.expand(rt_Gprime)

    # Match order by order. Solve the simultaneous system of coefficient
    # equations at u^1, u^2, u^3 = 0.
    def solve_system(eq_full):
        eqs = []
        for power in range(1, 4):
            c = sp.expand(eq_full).coeff(u, power)
            eqs.append(c)
        sol = sp.solve(eqs, [a1, a2, a3], dict=True)
        return sol[0] if sol else {}

    eq_GR = sp.expand(lhs - rhs_GR)
    eq_MA = sp.expand(lhs - rhs_MA)
    a_solved_GR = solve_system(eq_GR)
    a_solved_MA = solve_system(eq_MA)

    print("Isotropic G(r̃) = 1 + α₁ u + α₂ u² + α₃ u³ + ...")
    print()
    print(f"  GR solution:       α₁ = {a_solved_GR.get(a1)},  α₂ = {a_solved_GR.get(a2)},  α₃ = {a_solved_GR.get(a3)}")
    print(f"  Model-A solution:  α₁ = {a_solved_MA.get(a1)},  α₂ = {a_solved_MA.get(a2)},  α₃ = {a_solved_MA.get(a3)}")
    print()

    # --- Step 3: read off γ from α₁ ---
    # PPN: g_ij = (1 + 2γ u + ...) δ_ij  →  γ = α₁ / 2
    gamma_GR = a_solved_GR.get(a1) / 2
    gamma_MA = a_solved_MA.get(a1) / 2

    # --- Step 4: read off β from g_tt ---
    # In isotropic coords: g_tt = -(1 - A) c² = -(1 - 2u/√G) c²
    # Expand 1/√G with α's substituted:
    G_GR = G.subs(a_solved_GR)
    G_MA = G.subs(a_solved_MA)
    inv_sqrt_G_GR = sp.series(1/sp.sqrt(G_GR), u, 0, 4).removeO()
    inv_sqrt_G_MA = sp.series(1/sp.sqrt(G_MA), u, 0, 4).removeO()

    A_iso_GR = sp.expand(2*u * inv_sqrt_G_GR)
    A_iso_MA = sp.expand(2*u * inv_sqrt_G_MA)

    g_tt_over_c2_GR = -(1 - A_iso_GR)
    g_tt_over_c2_MA = -(1 - A_iso_MA)

    # PPN form: -(1 - 2u + 2β u² - ...).  Coefficient of u²: 2β.
    g_tt_GR_expanded = sp.series(g_tt_over_c2_GR, u, 0, 4).removeO()
    g_tt_MA_expanded = sp.series(g_tt_over_c2_MA, u, 0, 4).removeO()
    g_tt_GR_expanded = sp.expand(g_tt_GR_expanded)
    g_tt_MA_expanded = sp.expand(g_tt_MA_expanded)

    coeff_u2_GR = sp.Poly(g_tt_GR_expanded, u).coeff_monomial(u**2)
    coeff_u2_MA = sp.Poly(g_tt_MA_expanded, u).coeff_monomial(u**2)
    # PPN convention: g_tt = -(1 - 2u + 2β u² + ...), so the coefficient of u² in g_tt is +2β.
    # We have g_tt = -(1 - 2u + ...), so we need to read coefficient of u² in -(1-A):
    # -(1 - A) = -1 + A. Coefficient of u² in (-1 + A) is the +A coefficient of u².
    # In PPN: -(1 - 2u + 2β u²) = -1 + 2u - 2β u². So coefficient of u² in g_tt is -2β.
    # We computed coeff_u2 in g_tt_over_c2 = -(1-A) = -1 + A. So if A = 2u + c₂ u² + ...,
    # then g_tt has coefficient c₂ at u². And PPN says coefficient of u² in g_tt is -2β.
    # Therefore: -2β = c₂, i.e., β = -c₂/2.
    beta_GR = -coeff_u2_GR / 2
    beta_MA = -coeff_u2_MA / 2

    print("=" * 78)
    print("PPN PARAMETERS")
    print("=" * 78)
    print()
    print(f"{'Theory':>10s}    {'gamma (PPN)':>14s}    {'beta (PPN)':>14s}    "
          f"{'alpha2 (g_ij U^2 coeff)':>26s}")
    print("-" * 78)
    print(f"{'GR':>10s}    {str(sp.nsimplify(gamma_GR)):>14s}    "
          f"{str(sp.nsimplify(beta_GR)):>14s}    {str(a_solved_GR.get(a2)):>26s}")
    print(f"{'Model-A':>10s}    {str(sp.nsimplify(gamma_MA)):>14s}    "
          f"{str(sp.nsimplify(beta_MA)):>14s}    {str(a_solved_MA.get(a2)):>26s}")
    print()

    # --- Step 5: observational comparison ---
    print("=" * 78)
    print("OBSERVATIONAL COMPARISON")
    print("=" * 78)
    print()
    print(f"  Cassini bound (γ − 1):       |γ_obs − 1| < 2.3 × 10⁻⁵")
    print(f"    GR prediction:             γ − 1 = {sp.simplify(gamma_GR - 1)}")
    print(f"    Model-A prediction:        γ − 1 = {sp.simplify(gamma_MA - 1)}")
    print()
    print(f"  Mercury perihelion bound:    |β − 1| < few × 10⁻⁴")
    print(f"    GR prediction:             β − 1 = {sp.simplify(beta_GR - 1)}")
    print(f"    Model-A prediction:        β − 1 = {sp.simplify(beta_MA - 1)}")
    print()

    # --- Step 6: where Model-A first differs from GR ---
    print("=" * 78)
    print("WHERE MODEL-A FIRST DIFFERS FROM GR")
    print("=" * 78)
    print()
    a2_diff = sp.simplify(a_solved_MA.get(a2) - a_solved_GR.get(a2))
    a3_diff = sp.simplify(a_solved_MA.get(a3) - a_solved_GR.get(a3))
    print(f"  α₁ (g_ij first-order coefficient = 2γ):")
    print(f"    GR:       {a_solved_GR.get(a1)}")
    print(f"    Model-A:  {a_solved_MA.get(a1)}")
    print(f"    Difference: {sp.simplify(a_solved_MA.get(a1) - a_solved_GR.get(a1))}  "
          f"(matches at first order — γ identical)")
    print()
    print(f"  α₂ (g_ij second-order coefficient — third-order PPN territory):")
    print(f"    GR:       {a_solved_GR.get(a2)}")
    print(f"    Model-A:  {a_solved_MA.get(a2)}")
    print(f"    Difference: {a2_diff}")
    if a_solved_GR.get(a2) != 0:
        ratio = sp.simplify(a_solved_MA.get(a2) / a_solved_GR.get(a2))
        print(f"    Ratio Model-A / GR: {ratio}")
    print()

    print(f"  α₃ (g_ij third-order coefficient):")
    print(f"    GR:       {a_solved_GR.get(a3)}")
    print(f"    Model-A:  {a_solved_MA.get(a3)}")
    print(f"    Difference: {a3_diff}")
    print()

    # Schwarzschild-coord ratio (the "thirds" hint)
    print(f"  Schwarzschild-coordinate g_rr second-order coefficient (U_Sch²):")
    print(f"    GR:       {coeff_GR_2}")
    print(f"    Model-A:  {coeff_MA_2}")
    print(f"    Ratio Model-A / GR: {sp.simplify(coeff_MA_2 / coeff_GR_2)}")
    print(f"    Note: this is coordinate-dependent, but the integer ratio is suggestive.")
    print()

    # --- Step 7: write summary ---
    write_summary(gamma_GR, gamma_MA, beta_GR, beta_MA,
                   a_solved_GR, a_solved_MA, coeff_GR_2, coeff_MA_2,
                   a1, a2, a3)
    print(f"Summary: results/G17_ppn_parameters_summary.md")


def write_summary(gamma_GR, gamma_MA, beta_GR, beta_MA,
                   a_solved_GR, a_solved_MA, coeff_GR_2, coeff_MA_2,
                   a1, a2, a3) -> None:
    md = []
    md.append("# G17: PPN Parameter Calculation for Model-A\n")
    md.append("## Setup\n")
    md.append(
        "Compute the parameterized post-Newtonian (PPN) parameters γ and β "
        "for Model-A's strong-field metric and compare to GR Schwarzschild "
        "and to current solar-system observational bounds. Non-distance test "
        "— probes weak-field deviations of Model-A from GR using local "
        "spacetime structure.\n"
        "\n"
        "Model-A metric:\n"
        "```text\n"
        "g_tt = -(1 − A) c²\n"
        "g_rr = 1 / [(1 − A)(1 − A²)²]\n"
        "g_θθ = r²,  g_φφ = r² sin²θ\n"
        "A = 2GM/(c²r)\n"
        "```\n"
        "\n"
        "PPN formalism requires isotropic coordinates. The coordinate "
        "transformation r → r̃ is solved symbolically (sympy) to high enough "
        "order to extract γ and β.\n"
    )

    md.append("## Result\n")
    md.append(f"| Theory | γ (PPN) | β (PPN) | α₂ (g_ij U² coeff) |\n")
    md.append("|---|:---:|:---:|:---:|\n")
    md.append(f"| GR | {sp.nsimplify(gamma_GR)} | {sp.nsimplify(beta_GR)} | {a_solved_GR.get(a2)} |\n")
    md.append(f"| Model-A | {sp.nsimplify(gamma_MA)} | {sp.nsimplify(beta_MA)} | {a_solved_MA.get(a2)} |\n")
    md.append("\n")

    md.append("## Reading\n")
    md.append(
        f"**γ (Cassini-bounded at 10⁻⁵):**\n"
        f"- GR: γ − 1 = {sp.simplify(gamma_GR - 1)}\n"
        f"- Model-A: γ − 1 = {sp.simplify(gamma_MA - 1)}\n"
        "\n"
        "Model-A predicts γ identical to GR. Cassini bound is automatically "
        "satisfied. Light deflection and Shapiro delay match GR exactly at "
        "the level current measurements can distinguish.\n"
        "\n"
        f"**β (Mercury-bounded at 10⁻⁴):**\n"
        f"- GR: β − 1 = {sp.simplify(beta_GR - 1)}\n"
        f"- Model-A: β − 1 = {sp.simplify(beta_MA - 1)}\n"
        "\n"
        "Model-A predicts β identical to GR. Mercury perihelion shift, "
        "lunar laser ranging, and other β-dependent solar-system tests "
        "match GR exactly.\n"
        "\n"
        "**Conclusion on first/second-order PPN:** Model-A passes every "
        "current solar-system PPN constraint automatically. The framework's "
        "modified k(A) does not show up in any PPN-bounded observable.\n"
    )

    md.append("## Where Model-A first differs from GR\n")
    md.append(
        "Model-A's deviation from GR appears in the third-order coefficient "
        "of the spatial isotropic metric (α₂ in `g_ij = (1 + 2γu + α₂ u² + α₃ u³ + ...) δ_ij`):\n"
        "\n"
        f"- GR α₂ = {a_solved_GR.get(a2)}\n"
        f"- Model-A α₂ = {a_solved_MA.get(a2)}\n"
        f"- Difference: {sp.simplify(a_solved_MA.get(a2) - a_solved_GR.get(a2))}\n"
        "\n"
        "This is **third-order PPN territory**, where solar-system "
        "constraints are essentially absent (no current measurement reaches "
        "this precision). Model-A is observationally indistinguishable from "
        "GR at all current solar-system tests, but the deviation is real "
        "and would show up in any future measurement that reaches third-order "
        "PPN sensitivity.\n"
    )

    md.append("## Coordinate-dependent structure (worth noting)\n")
    md.append(
        f"In **Schwarzschild coordinates**, the second-order coefficient of "
        f"`g_rr` (in U_Sch = GM/(c²r)) is:\n"
        f"- GR: {coeff_GR_2}\n"
        f"- Model-A: {coeff_MA_2}\n"
        f"- Ratio: {sp.simplify(coeff_MA_2 / coeff_GR_2)}\n"
        "\n"
        "**Model-A's Schwarzschild-coordinate g_rr has exactly 3× the "
        "second-order coefficient of GR.** This is coordinate-dependent — in "
        "isotropic coordinates the ratio shifts away from exactly 3 — but "
        "the appearance of an integer factor of 3 in the natural "
        "Schwarzschild-coordinate expansion is suggestive in light of the "
        "framework's other 'thirds' structure (ISCO/photon-sphere/horizon "
        "at A = 1/3, 2/3, 1).\n"
        "\n"
        "Whether this is the same '3' appearing structurally in different "
        "places is the open question we've been circling. This script "
        "doesn't answer it — it surfaces another instance and notes it.\n"
    )

    md.append("## Implications\n")
    md.append(
        "- **Model-A is solar-system-safe.** No current PPN test constrains "
        "the framework. Every test that has bounded γ or β to 10⁻⁵–10⁻⁴ "
        "precision is automatically passed.\n"
        "- **Model-A is solar-system-indistinguishable from GR.** The "
        "framework's modified k(A) does not produce any second-order PPN "
        "deviation. The first signature appears at third-order, currently "
        "unconstrained.\n"
        "- **The falsifiable corners are not in the solar system.** As the "
        "strong-field commitment memory already states, the genuine tests "
        "are LIGO QNM ringdown (G1's τ_MA/τ_GR = 1.80), EHT shadow "
        "structure beyond size, and any near-horizon Shapiro-class "
        "measurement at A → 1.\n"
        "- **The 3× ratio in Schwarzschild g_rr** is coordinate-dependent "
        "but worth tracking. If a coordinate-invariant statement of this "
        "structure exists, it would be another '3' to add to the framework's "
        "list — alongside the strong-field thirds and the (4π × 3) "
        "decomposition of A_0. Whether they are the same 3 is open.\n"
    )

    out = RESULTS / "G17_ppn_parameters_summary.md"
    out.write_text("".join(md), encoding="utf-8")


if __name__ == "__main__":
    main()
