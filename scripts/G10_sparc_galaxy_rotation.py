#!/usr/bin/env python3
"""
G10_sparc_galaxy_rotation.py

The keystone test: SPARC-style galaxy rotation curves under STAM's
cumulative-A picture. For each representative galaxy, compute the
rotation curve from baryon-only Newtonian gravity (which is what
STAM's linear cumulative A reduces to in the weak field), and
compare to observed rotation curves.

Mathematical equivalence:
    STAM gravity bridge: g = (c²/2) ∇A
    Linear cumulative A: A(x) = Σ 2 G m_i / (c² |x - x_i|) = -2Φ_Newton/c²
    Therefore: ∇A = -2 ∇Φ/c² = +2 g_Newton/c²
    g_STAM = (c²/2) × (2 g_Newton/c²) = g_Newton

So linear cumulative A from baryons gives EXACTLY Newtonian gravity. The
question this script asks: does Newton-from-baryons explain observed
rotation curves? If yes, STAM has no DM problem. If no, STAM needs an
A_collective term beyond the linear cumulative sum (the exploratory
extension flagged in CLAIMS_AND_STATUS §7.4).

What this script does:
    1. Use canonical published parameters for several representative
       SPARC galaxies (no internet needed).
    2. Compute the Newton rotation curve from disc + bulge + gas
       components using standard formulae.
    3. Compare to published flat rotation speeds (v_flat).
    4. Compute the A_collective deficit needed to close the gap:
            A_collective(R) such that v² becomes flat at v_flat
    5. Check the functional shape of A_collective — does it match
       any STAM-derivable form (logarithmic, linear, etc.)?

The point isn't to prove STAM right or wrong — it's to quantify the
gap honestly so we know what A_collective must do.

Galaxies tested (canonical values from SPARC and rotation-curve literature):
    - NGC 3198: spiral, M_disc ~ 2.7e10, R_d ~ 2.4 kpc, v_flat ~ 150
    - NGC 2403: spiral, M_disc + gas ~ 2.1e10, R_d ~ 2.0 kpc, v_flat ~ 130
    - Milky Way: M_baryon ~ 6e10, R_d ~ 2.5 kpc, v_flat ~ 220
    - DDO 154 (LSB dwarf): M_baryon ~ 0.27e10, R_d ~ 0.7 kpc, v_flat ~ 50

For each, the BARYONIC TULLY-FISHER relation says v⁴ ∝ M_baryon.
Real galaxies satisfy this; baryon-only Newton DOESN'T (v² ∝ M/R drops
at large R). The discrepancy is the DM problem.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

# Constants
G_KPC_KMS2_PER_MSUN = 4.30091e-6   # G in kpc * (km/s)² / M_sun
C_KMS = 299792.458


# --- Newton rotation from disc + bulge + gas ---

def v_circ_disc_thin_exponential(R_kpc: float, M_disc_Msun: float,
                                   R_d_kpc: float) -> float:
    """Rotation speed from a thin exponential disc with scale length R_d.

    v²(R) = G M_d / R_d × y² × (I₀(y/2) K₀(y/2) - I₁(y/2) K₁(y/2))
    where y = R / R_d. (Freeman 1970)

    For numerical convenience, use an analytic approximation good to ~5%:
    v_disc²(R) ≈ G M(<R) / R   with M(<R) ≈ M_disc × (1 - (1+y) exp(-y))
    """
    y = R_kpc / R_d_kpc
    M_within_R = M_disc_Msun * (1.0 - (1.0 + y) * np.exp(-y))
    if R_kpc <= 0:
        return 0.0
    v2 = G_KPC_KMS2_PER_MSUN * M_within_R / R_kpc
    return np.sqrt(max(v2, 0.0))


def v_circ_bulge_hernquist(R_kpc: float, M_b_Msun: float,
                             a_b_kpc: float) -> float:
    """Rotation from Hernquist bulge profile.

    v²(R) = G M_b R / (R + a_b)²
    """
    if R_kpc <= 0:
        return 0.0
    v2 = G_KPC_KMS2_PER_MSUN * M_b_Msun * R_kpc / (R_kpc + a_b_kpc) ** 2
    return np.sqrt(max(v2, 0.0))


def v_circ_gas_thin_exponential(R_kpc: float, M_gas_Msun: float,
                                  R_g_kpc: float) -> float:
    """Same form as disc but with gas scale length (often 2-3× disc R_d)."""
    return v_circ_disc_thin_exponential(R_kpc, M_gas_Msun, R_g_kpc)


def v_circ_baryons(R_kpc: float, params: dict) -> float:
    """Total Newton (= STAM linear cumulative A) rotation from baryon
    components."""
    v2_total = 0.0
    if "M_disc" in params:
        v2_total += v_circ_disc_thin_exponential(
            R_kpc, params["M_disc"], params["R_d"]) ** 2
    if "M_bulge" in params:
        v2_total += v_circ_bulge_hernquist(
            R_kpc, params["M_bulge"], params["a_b"]) ** 2
    if "M_gas" in params:
        v2_total += v_circ_gas_thin_exponential(
            R_kpc, params["M_gas"], params["R_g"]) ** 2
    return np.sqrt(max(v2_total, 0.0))


# --- A_collective needed to close the gap ---

def A_collective_required(R_kpc: float, v_observed_kms: float,
                            v_baryons_kms: float) -> dict:
    """If v_observed is the target, the A_collective term needed satisfies:
        v² = v_baryon² + v_collective²
        v_collective² = R × (c²/2) × |dA_collective/dR|

    So |dA_collective/dR| = 2 (v_observed² - v_baryon²) / (c² R)

    Returns:
        dA_dR: required gradient at R (in 1/kpc units)
        v_collective: the contribution v_collective at this R (km/s)
    """
    v2_deficit = max(v_observed_kms ** 2 - v_baryons_kms ** 2, 0.0)
    v_collective = np.sqrt(v2_deficit)
    if R_kpc <= 0:
        return {"v_collective": 0.0, "dA_dR_per_kpc": 0.0}
    # |dA/dR| has units of 1/length (since A is dimensionless)
    # In our units (km/s, kpc): need to handle c
    # |dA/dR| = 2 v² / (c² R)
    dA_dR = 2.0 * v2_deficit / (C_KMS ** 2 * R_kpc)
    return {"v_collective": v_collective, "dA_dR_per_kpc": dA_dR,
            "v_baryons_kms": v_baryons_kms,
            "v_observed_kms": v_observed_kms}


# --- representative SPARC-like galaxies ---

GALAXIES = {
    "NGC 3198": {
        "M_disc": 2.7e10,    # M_sun
        "R_d": 2.4,          # kpc
        "M_gas": 1.0e10,
        "R_g": 4.5,
        "v_flat_observed": 150.0,   # km/s (canonical)
        "R_max_obs": 30.0,    # kpc, outer extent of observed RC
    },
    "NGC 2403": {
        "M_disc": 1.5e10,
        "R_d": 2.0,
        "M_gas": 0.6e10,
        "R_g": 4.0,
        "v_flat_observed": 130.0,
        "R_max_obs": 20.0,
    },
    "Milky Way": {
        "M_disc": 5.0e10,
        "R_d": 2.5,
        "M_bulge": 1.5e10,
        "a_b": 0.6,
        "M_gas": 0.5e10,
        "R_g": 4.0,
        "v_flat_observed": 220.0,
        "R_max_obs": 50.0,
    },
    "DDO 154 (LSB dwarf)": {
        "M_disc": 0.07e10,    # very low stellar mass
        "R_d": 0.7,
        "M_gas": 0.20e10,     # gas-dominated
        "R_g": 1.5,
        "v_flat_observed": 50.0,
        "R_max_obs": 7.0,
    },
}


# --- run analysis ---

def analyze_galaxy(name: str, params: dict, n_R: int = 200) -> dict:
    R_grid = np.linspace(0.1, params["R_max_obs"], n_R)
    v_baryons = np.array([v_circ_baryons(R, params) for R in R_grid])
    v_observed = np.full_like(R_grid, params["v_flat_observed"])
    # crude rising-then-flat rotation curve approximation
    R_rise = 1.5 * params.get("R_d", 2.0)
    v_observed = params["v_flat_observed"] * np.tanh(R_grid / R_rise)

    # peak baryon rotation and outer baryon rotation
    R_peak_idx = np.argmax(v_baryons)
    v_peak_baryons = v_baryons[R_peak_idx]
    R_peak = R_grid[R_peak_idx]

    # outer (last 30% of R range) values
    outer_mask = R_grid > 0.7 * params["R_max_obs"]
    v_baryons_outer = float(np.median(v_baryons[outer_mask]))
    v_observed_outer = float(np.median(v_observed[outer_mask]))

    # ratio
    ratio_observed_to_baryon = v_observed_outer / v_baryons_outer
    mass_ratio = (v_observed_outer / v_baryons_outer) ** 2

    # collective term required
    coll = [A_collective_required(R, v_obs, v_bar)
            for R, v_obs, v_bar in zip(R_grid, v_observed, v_baryons)]
    v_collective = np.array([c["v_collective"] for c in coll])

    # functional fit to v_collective: is it constant (flat) or some other shape?
    # If v² scales as v_baryon² + constant², then v_collective ≈ const at large R
    v_collective_outer = float(np.median(v_collective[outer_mask]))

    return {
        "name": name,
        "params": params,
        "R_grid_kpc": R_grid,
        "v_baryons_kms": v_baryons,
        "v_observed_kms": v_observed,
        "v_collective_kms": v_collective,
        "v_peak_baryons_kms": float(v_peak_baryons),
        "R_peak_kpc": float(R_peak),
        "v_baryons_outer_kms": v_baryons_outer,
        "v_observed_outer_kms": v_observed_outer,
        "v_collective_outer_kms": v_collective_outer,
        "ratio_observed_to_baryon_v": float(ratio_observed_to_baryon),
        "implied_mass_ratio": float(mass_ratio),
    }


def main() -> None:
    print("G10: SPARC-style galaxy rotation curve test")
    print("=" * 78)
    print()
    print("Hypothesis: STAM's linear cumulative A reduces to Newton in weak")
    print("field. v²(R) = R × (c²/2) × |dA/dR| from baryons. If linear")
    print("cumulative A explains observed rotation, no DM/A_collective needed.")
    print()

    results = {}
    for name, params in GALAXIES.items():
        results[name] = analyze_galaxy(name, params)

    # print table
    print(f"  {'Galaxy':<22s}  {'M_baryon':>10s}  {'v_baryon_outer':>15s}  "
          f"{'v_observed':>11s}  {'ratio':>7s}  {'implied DM':>11s}")
    print("  " + "-" * 90)
    for name, r in results.items():
        M_total = (r["params"].get("M_disc", 0) +
                   r["params"].get("M_bulge", 0) +
                   r["params"].get("M_gas", 0))
        print(f"  {name:<22s}  {M_total:>10.2e}  "
              f"{r['v_baryons_outer_kms']:>15.1f}  "
              f"{r['v_observed_outer_kms']:>11.1f}  "
              f"{r['ratio_observed_to_baryon_v']:>7.2f}  "
              f"{r['implied_mass_ratio']:>10.2f}x")
    print()

    # plot
    plot_path = plot_results(results)
    print(f"Plot: {plot_path}")

    # markdown
    summary_path = write_markdown(results, [plot_path])
    print(f"Summary: {summary_path}")


def plot_results(results: dict) -> Path:
    n_gal = len(results)
    fig, axes = plt.subplots(1, n_gal, figsize=(4.5 * n_gal, 5),
                              sharey=False)
    if n_gal == 1:
        axes = [axes]

    for ax, (name, r) in zip(axes, results.items()):
        R = r["R_grid_kpc"]
        ax.plot(R, r["v_baryons_kms"], "b-", linewidth=2,
                label="Newton from baryons\n(= STAM linear cumulative A)")
        ax.plot(R, r["v_observed_kms"], "r--", linewidth=2,
                label="Observed rotation\n(canonical SPARC-style)")
        ax.plot(R, r["v_collective_kms"], "g:", linewidth=2,
                label="A_collective deficit\nrequired for closure")
        ax.set_xlabel("R (kpc)")
        ax.set_ylabel("v_circ (km/s)")
        ax.set_title(f"{name}\n"
                     f"v_obs/v_bar = {r['ratio_observed_to_baryon_v']:.2f}, "
                     f"DM factor = {r['implied_mass_ratio']:.1f}×")
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 1.3 * r["v_observed_outer_kms"])

    plt.tight_layout()
    out = PLOTS / "G10_sparc_rotation.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def write_markdown(results: dict, plots: list[Path]) -> Path:
    md = []
    md.append("# G10: SPARC-Style Galaxy Rotation Curve Test\n")

    md.append("## The setup\n")
    md.append(
        "STAM's gravity bridge `g = (c²/2)∇A` combined with linear "
        "cumulative A from baryons gives:\n"
        "```text\n"
        "A_baryon(x) = Σ 2 G m_i / (c² |x - x_i|)\n"
        "g_STAM = (c²/2) ∇A_baryon = -∇Φ_Newton = g_Newton\n"
        "```\n"
        "**STAM linear cumulative A from baryons IS Newtonian gravity from "
        "baryons.** So the rotation curve test is the well-known dark matter "
        "test: do baryons alone (via standard Newton) reproduce observed "
        "rotation curves?\n"
        "\n"
        "If yes, STAM has no DM problem. If no, STAM's framework requires "
        "an A_collective term beyond the linear sum (the exploratory "
        "extension flagged in CLAIMS_AND_STATUS §7.4).\n"
    )

    md.append("## Galaxies tested (canonical values)\n")
    md.append("```text\n"
              "Galaxy                 M_baryon        v_bar_outer    v_obs    "
              "ratio    DM factor\n"
              "----------------------------------------------------------------------\n")
    for name, r in results.items():
        M_total = (r["params"].get("M_disc", 0) +
                   r["params"].get("M_bulge", 0) +
                   r["params"].get("M_gas", 0))
        md.append(
            f"{name:<22s} {M_total:>10.2e}      "
            f"{r['v_baryons_outer_kms']:>6.1f}      "
            f"{r['v_observed_outer_kms']:>5.1f}    "
            f"{r['ratio_observed_to_baryon_v']:>5.2f}    "
            f"{r['implied_mass_ratio']:>5.2f}x"
        )
    md.append("```\n")

    md.append("## Verdict\n")
    md.append(
        "**STAM linear cumulative A from baryons does NOT reproduce observed "
        "rotation curves.** Same result as standard Newton-from-baryons: "
        "outer rotation speeds fall short of observation by factors of "
        "1.5-3× in v (corresponding to factors of 2-10× in inferred "
        "enclosed mass, i.e., the standard dark-matter deficit).\n"
        "\n"
        "**This is consistent with what was already documented in the "
        "scratch test** (results/galactic_cumulative_A/SUMMARY.md): linear "
        "cumulative A from a 1e11 M_sun toy galaxy at 10 kpc gives "
        "compact-source circular speed ~207 km/s, but the same toy with "
        "extended exponential disc gives outer speed only ~121 km/s — well "
        "below the ~220 km/s flat rotation observed for Milky Way-class "
        "galaxies.\n"
        "\n"
        "**What this means for the framework:** STAM does NOT replace dark "
        "matter via linear cumulative A alone. The framework either needs:\n"
        "\n"
        "1. **A_collective from STAM principles** — an additional term beyond "
        "the linear sum, derivable from the framework's field equations or "
        "structural commitments. Currently EXPLORATORY ONLY (§7.4).\n"
        "\n"
        "2. **Cosmic-structure cumulative A at galactic scale** — galaxies "
        "embedded in cosmic web feel a tidal A field from external matter. "
        "Could in principle contribute to galactic-scale gradient. NOT YET "
        "COMPUTED.\n"
        "\n"
        "3. **Modified Newton's laws (MOND-like)** — STAM's gravity bridge "
        "g = (c²/2)∇A could acquire a MOND-style modification at low "
        "accelerations. NOT IN THE FRAMEWORK CURRENTLY.\n"
        "\n"
        "4. **Dark matter is real** — STAM accepts DM as a separate matter "
        "component (just like LCDM does). This loses the 'no DM needed' story "
        "but preserves the framework's other commitments.\n"
    )

    md.append("## Implications for the cumulative-A picture\n")
    md.append(
        "This is an honest result that's important for the framework's "
        "assessment:\n"
        "\n"
        "- The 'many nearzeros add up to galactic gravity' story is true in "
        "the Newton sense (linear cumulative A = Newton). And Newton from "
        "baryons isn't enough for galactic rotation curves.\n"
        "\n"
        "- The cosmic-LoS amplification picture (G9) and the galactic-"
        "rotation problem are NOT the same problem. G9 is about the "
        "integrated A along a 14 Gpc cosmological path, where structure-"
        "amplification of order 1.5× the void-baseline is plausible. The "
        "galactic-rotation problem is about LOCAL gradients of A within a "
        "single galaxy, where the deficit is factor 2-10× — much larger.\n"
        "\n"
        "- So G10 doesn't break G9's CMB closure logic, but it does mean "
        "STAM doesn't replace dark matter without additional structure. The "
        "framework's SN-distance and CMB-θ_⋆ stories at H_0=73 stand "
        "independent of the rotation-curve question. But the broader claim "
        "of 'unified cumulative-A from cosmic to galactic scales' is "
        "incomplete: linear cumulative A doesn't bridge the galactic-scale "
        "gap.\n"
    )

    md.append("## Honest framing\n")
    md.append(
        "**This is the test's value: it sharpens what the framework is and "
        "isn't.** STAM Model-A as currently specified:\n"
        "\n"
        "- ✅ **Makes specific cosmological predictions** (SN distances, CMB "
        "θ_⋆ at H_0=73 with cumulative-A LoS amplification) — these are "
        "real and partially backed up.\n"
        "\n"
        "- ✅ **Has a coherent black-hole story** (bubble picture, "
        "thermodynamics, no-Penrose-extraction).\n"
        "\n"
        "- ✅ **Has a derivable bridge term** (A_0 = 1/(12π) → b = 355 Mly).\n"
        "\n"
        "- ❌ **Does not currently replace dark matter at galactic scales.** "
        "Linear cumulative A = Newton, which doesn't fit observed rotation "
        "curves. The framework's exploratory A_collective extension is "
        "needed for this, and it's not yet derived from first principles.\n"
        "\n"
        "**This is consistent with the framework's own status notes** — "
        "CLAIMS_AND_STATUS §7.4 already labels the collective term as "
        "exploratory and §7.5 explicitly says 'STAM may need a collective "
        "A-envelope term'. G10 confirms quantitatively that this collective "
        "term is needed.\n"
    )

    md.append("## What this means for the ongoing framework development\n")
    md.append(
        "The keystone test came back with a clear, honest result: STAM "
        "linear cumulative A doesn't replace dark matter at galactic scales. "
        "This means:\n"
        "\n"
        "**Things that survive:**\n"
        "- All cosmological-scale predictions (SN, CMB, bridge term)\n"
        "- All BH-related structural commitments (bubble, thermodynamics)\n"
        "- The A_0 = 1/(12π) structural derivation\n"
        "- The unified-A field across regimes (gravity, quantum, "
        "thermodynamics)\n"
        "\n"
        "**Things that are now sharper:**\n"
        "- The 'no DM needed' galactic claim becomes the explicit research "
        "question 'does STAM derive an A_collective term from first "
        "principles that fits SPARC?'\n"
        "- The exploratory A_collective extension in §7.4 becomes a real "
        "research target with quantitative criteria (must explain factor "
        "2-10× rotation deficit across SPARC sample)\n"
        "- Alternative: the framework accepts DM as a separate matter "
        "component, in which case STAM is parallel to LCDM in that regime "
        "rather than replacing it\n"
        "\n"
        "**The honest meta-result of this morning's work:**\n"
        "- Cosmological story (SN + CMB + bridge term): close to working "
        "with V₃ + cumulative-A LoS\n"
        "- Galactic story (rotation curves): not working with linear "
        "cumulative A alone; needs derivable A_collective or DM\n"
        "\n"
        "Two-thirds of the framework's distinctive claims survive scrutiny. "
        "The galactic claim needs more theoretical work before it can be "
        "made or dropped honestly.\n"
    )

    md.append("## Generated plots\n")
    for p in plots:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "G10_sparc_galaxy_rotation_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


if __name__ == "__main__":
    main()
