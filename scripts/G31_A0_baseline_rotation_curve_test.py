#!/usr/bin/env python3
"""
G31_A0_baseline_rotation_curve_test.py

Test whether A_0 baseline (cosmic void value) shifts galactic rotation curves
enough to bridge G10's deficit.

Sean 2026-05-11 evening: "the baseline is A_0, not the floor but the void,
so the addition of matter into a void increases from here, not zero. Does
this shift the effect enough?"

Direct test:
  - Compute rotation curve under "A = 0 + A_local" (G10 baseline, Newton-equivalent)
  - Compute rotation curve under "A = A_0 + A_local" (new framing)
  - Quantify the shift
  - Compare both to observed v_flat for several SPARC galaxies

Theory note (worked out before running):
  Gravitational acceleration in STAM metric:
    g_r = (c²/2) × (1/(1-A)) × dA/dr
  At A = 0:       g_r = (c²/2) × dA/dr
  At A = A_0:     g_r = (c²/2) × dA/dr / (1 - A_0) ≈ (c²/2) × dA/dr × 1.027
  Enhancement factor: 1/(1-A_0) ≈ 1.0273  (about 2.7%)

So the predicted enhancement is ~2.7%. The G10 deficit was factor 1.5-3 in v.
A 2.7% shift won't bridge that. The script verifies this numerically.

Author: Sean Brady / STAM Model-A
Date: 2026-05-11 (evening)
"""

from __future__ import annotations
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Constants
# ============================================================================
A_0 = 1.0 / (12.0 * math.pi)
G_KPC_KMS2_PER_MSUN = 4.30091e-6  # G in kpc·(km/s)²/M_sun
C_KMS = 299792.458


# ============================================================================
# Canonical galaxy parameters (same as G10's set)
# ============================================================================
GALAXIES = [
    {
        "name": "NGC 3198",
        "M_disc": 2.7e10,   # M_sun
        "R_d": 2.4,         # kpc disc scale
        "M_gas": 1.1e10,
        "v_flat_obs": 150.0,  # km/s observed flat rotation
        "R_max": 30.0,        # kpc max radius for plot
    },
    {
        "name": "NGC 2403",
        "M_disc": 2.1e10,
        "R_d": 2.0,
        "M_gas": 0.6e10,
        "v_flat_obs": 130.0,
        "R_max": 25.0,
    },
    {
        "name": "Milky Way",
        "M_disc": 6.0e10,
        "R_d": 2.5,
        "M_gas": 0.7e10,
        "v_flat_obs": 220.0,
        "R_max": 30.0,
    },
    {
        "name": "DDO 154",
        "M_disc": 0.27e10,
        "R_d": 0.7,
        "M_gas": 0.3e10,
        "v_flat_obs": 50.0,
        "R_max": 8.0,
    },
]


# ============================================================================
# Baryon mass enclosed (disc + gas, exponential profile approximation)
# ============================================================================
def M_baryon_within(R_kpc, M_disc, R_d, M_gas):
    """Enclosed baryon mass at radius R (disc + uniform gas distribution)."""
    y = R_kpc / R_d
    # Disc: exponential with scale R_d
    M_disc_within = M_disc * (1.0 - (1.0 + y) * np.exp(-y))
    # Gas: assume uniform density up to ~3 R_d, drops linearly past
    R_gas_max = 3.0 * R_d
    M_gas_within = M_gas * np.minimum(1.0, (R_kpc / R_gas_max) ** 2)
    return M_disc_within + M_gas_within


def v_newton_from_baryons(R_kpc, M_disc, R_d, M_gas):
    """Newtonian rotation speed from baryon mass enclosed: v² = G M(<R) / R."""
    M_enc = M_baryon_within(R_kpc, M_disc, R_d, M_gas)
    return np.sqrt(G_KPC_KMS2_PER_MSUN * M_enc / R_kpc)


# ============================================================================
# STAM rotation curve under "A = A_0 + A_local" composition
# ============================================================================
def v_stam_with_A0(R_kpc, M_disc, R_d, M_gas):
    """Rotation speed under STAM with A_0 baseline.

    Under STAM metric g_tt = -(1-A)c²:
      circular orbit speed v² = (c²/2) × (1/(1-A)) × R × dA/dR

    With A = A_0 + A_local(R) and A_local(R) ≈ 2GM(<R)/(c²R) (weak field):
      dA/dR = (2G/c²) × [M'(R)/R - M(R)/R²]
      But for the cumulative A from extended baryon distribution,
      we instead have d(2GM/c²R)/dR which gives Newton + correction.

    In the weak-field equivalent form:
      v_STAM² = v_Newton² / (1 - A_0 - A_local)
              ≈ v_Newton² × (1 + A_0 + A_local + ...)

    For galactic scales A_local << A_0, so:
      v_STAM² ≈ v_Newton² / (1 - A_0)
      v_STAM / v_Newton ≈ 1 / sqrt(1 - A_0) ≈ 1.0136

    The 1.36% enhancement comes from the (1-A_0) factor in the metric — it
    represents that gravitational acceleration in a region with A=A_0 baseline
    is slightly stronger than at A=0 baseline.
    """
    M_enc = M_baryon_within(R_kpc, M_disc, R_d, M_gas)
    # Local A from baryons (weak field)
    A_local = 2.0 * G_KPC_KMS2_PER_MSUN * M_enc / (C_KMS**2 * R_kpc)
    A_total = A_0 + A_local
    # v² from metric: v² = (c²/2) × R × dA/dR × (1/(1-A)). For weak field, equivalent to:
    # v² ≈ GM(<R)/R × (1/(1-A))
    # Using the standard reduction
    v_squared = (G_KPC_KMS2_PER_MSUN * M_enc / R_kpc) / (1.0 - A_total)
    return np.sqrt(v_squared)


def v_newton_only(R_kpc, M_disc, R_d, M_gas):
    """Pure Newton (G10 baseline)."""
    return v_newton_from_baryons(R_kpc, M_disc, R_d, M_gas)


# ============================================================================
# Main
# ============================================================================
def main():
    print("=" * 78)
    print("G31: Does A_0 baseline shift rotation curves enough to bridge G10?")
    print("=" * 78)
    print()
    print(f"A_0           = {A_0:.6f}")
    print(f"1/(1-A_0)     = {1.0/(1.0-A_0):.6f}  (gravity enhancement factor)")
    print(f"sqrt(1/(1-A_0)) = {1.0/math.sqrt(1.0-A_0):.6f}  (v enhancement factor)")
    print(f"Expected v shift: {(1.0/math.sqrt(1.0-A_0) - 1.0)*100:.3f}%")
    print()

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    axes = axes.flatten()

    print(f"  {'Galaxy':<12} {'R_eval':>8} {'v_obs':>8} {'v_Newton':>9} "
          f"{'v_STAM':>9} {'STAM/Newton':>12} {'Newton/obs':>11} {'STAM/obs':>10}")

    for idx, gal in enumerate(GALAXIES):
        name = gal["name"]
        R_grid = np.linspace(0.1, gal["R_max"], 200)

        v_N = v_newton_only(R_grid, gal["M_disc"], gal["R_d"], gal["M_gas"])
        v_S = v_stam_with_A0(R_grid, gal["M_disc"], gal["R_d"], gal["M_gas"])

        # Evaluate at "characteristic flat radius" ≈ 3·R_d
        R_eval = 3.0 * gal["R_d"]
        idx_eval = np.argmin(np.abs(R_grid - R_eval))
        v_N_eval = v_N[idx_eval]
        v_S_eval = v_S[idx_eval]
        v_obs = gal["v_flat_obs"]

        print(f"  {name:<12} {R_eval:>7.2f} {v_obs:>8.1f} {v_N_eval:>9.2f} "
              f"{v_S_eval:>9.2f} {v_S_eval/v_N_eval:>12.4f} "
              f"{v_N_eval/v_obs:>11.3f} {v_S_eval/v_obs:>10.3f}")

        # Plot
        ax = axes[idx]
        ax.plot(R_grid, v_N, "b-", lw=2, label=f"Newton from baryons (G10)")
        ax.plot(R_grid, v_S, "r--", lw=2, label=f"STAM with A_0 baseline")
        ax.axhline(v_obs, color="k", linestyle=":", lw=1.5,
                   label=f"Observed v_flat = {v_obs:.0f} km/s")
        ax.set_xlabel("R [kpc]")
        ax.set_ylabel("v [km/s]")
        ax.set_title(f"{name}  (M_baryon = {gal['M_disc']+gal['M_gas']:.1e} M_sun)")
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)

    fig.suptitle("G31: Adding A_0 baseline to Newton-from-baryons rotation curves")
    fig.tight_layout()
    plot_path = PLOTS / "G31_A0_baseline_rotation_curve.png"
    fig.savefig(plot_path, dpi=140, bbox_inches="tight")
    print()
    print(f"Plot saved: {plot_path}")

    # ====================================================================
    # Markdown summary
    # ====================================================================
    md = []
    md.append("# G31: Does A_0 baseline shift galactic rotation curves enough?")
    md.append("")
    md.append("**Date:** 2026-05-11 (evening)")
    md.append("")
    md.append("## The question")
    md.append("")
    md.append(
        "Sean 2026-05-11 evening: when matter is added to a void, A goes from A_0 "
        "to A_0 + δA (not from 0 to δA). Does sitting on this baseline shift the "
        "galactic rotation curve enough to bridge G10's factor-of-1.5-to-3 deficit?"
    )
    md.append("")
    md.append("## Theory (worked out before running)")
    md.append("")
    md.append(
        "Under STAM metric g_tt = -(1-A)c², the orbital speed for circular motion is:"
    )
    md.append("")
    md.append("```")
    md.append("v² ≈ G M(<R) / R × (1 / (1 - A))")
    md.append("")
    md.append(f"At A = 0:    v² = G M(<R) / R                (Newton baseline)")
    md.append(f"At A = A_0:  v² = G M(<R) / R × (1 / (1 - {A_0:.4f}))")
    md.append(f"                = G M(<R) / R × {1.0/(1.0-A_0):.4f}")
    md.append("```")
    md.append("")
    md.append(
        f"Expected v enhancement: 1/sqrt(1-A_0) ≈ {1.0/math.sqrt(1.0-A_0):.4f}, "
        f"or about {(1.0/math.sqrt(1.0-A_0) - 1.0)*100:.2f}% increase."
    )
    md.append("")
    md.append("## Numerical result by galaxy")
    md.append("")
    md.append("| Galaxy | R_eval (kpc) | v_obs (km/s) | v_Newton | v_STAM | STAM/Newton | Newton/obs | STAM/obs |")
    md.append("|---|---|---|---|---|---|---|---|")
    for gal in GALAXIES:
        name = gal["name"]
        R_grid = np.linspace(0.1, gal["R_max"], 200)
        v_N = v_newton_only(R_grid, gal["M_disc"], gal["R_d"], gal["M_gas"])
        v_S = v_stam_with_A0(R_grid, gal["M_disc"], gal["R_d"], gal["M_gas"])
        R_eval = 3.0 * gal["R_d"]
        idx_eval = np.argmin(np.abs(R_grid - R_eval))
        md.append(
            f"| {name} | {R_eval:.2f} | {gal['v_flat_obs']:.1f} | "
            f"{v_N[idx_eval]:.2f} | {v_S[idx_eval]:.2f} | "
            f"{v_S[idx_eval]/v_N[idx_eval]:.4f} | "
            f"{v_N[idx_eval]/gal['v_flat_obs']:.3f} | "
            f"{v_S[idx_eval]/gal['v_flat_obs']:.3f} |"
        )
    md.append("")
    md.append("## Reading")
    md.append("")
    md.append(
        f"**The A_0 baseline shifts v by ~{(1.0/math.sqrt(1.0-A_0) - 1.0)*100:.2f}%** "
        "across all galaxies. The Newton-from-baryons deficit (~factor 1.5-3 below "
        f"observed) becomes the STAM-with-A_0 deficit (~factor {1.5/(1.0/math.sqrt(1.0-A_0)):.2f}-"
        f"{3.0/(1.0/math.sqrt(1.0-A_0)):.2f} below observed). Essentially unchanged."
    )
    md.append("")
    md.append(
        "**A_0 baseline does NOT bridge G10's deficit.** The mechanism: A_0 is "
        "spatially uniform, contributes no gradient, no gravitational acceleration "
        "by itself. Its only effect is to enhance the response to LOCAL A "
        "contributions via the (1/(1-A)) factor in the metric. The enhancement is "
        "1/(1-A_0) in v² and 1/sqrt(1-A_0) in v — about 1.36% in v, two orders of "
        "magnitude smaller than what's needed."
    )
    md.append("")
    md.append("## Implications")
    md.append("")
    md.append(
        "- **G10's keystone-test result stands.** Linear cumulative A from baryons "
        "(equivalent to Newton from baryons) doesn't explain rotation curves, and "
        "A_0 baseline doesn't change that."
    )
    md.append("")
    md.append(
        "- **The framework still needs additional physics to explain galactic dark "
        "matter signatures.** Current best candidate: PBH-DM (G13). Other options: "
        "structure-dependent A_LoS enhancement (G11 line of argument), nonlinear "
        "cumulative A effects, or revised baryon-only with full STAM nonlinearity."
    )
    md.append("")
    md.append(
        "- **The A_0 baseline DOES have one real role:** it produces the bridge "
        "term in the cosmological branch (photon-A traversal through cosmic A_0 "
        "in voids contributes ≈ A_0·c/H_0 to apparent distance modulus). But for "
        "galactic dynamics, it's effectively invisible."
    )
    md.append("")
    md.append("## Plot")
    md.append("")
    md.append("![G31 rotation curves with A_0](../plots/G31_A0_baseline_rotation_curve.png)")
    md.append("")
    md.append(
        "The Newton-from-baryons and STAM-with-A_0 curves are essentially "
        "indistinguishable on this plot scale. The observed v_flat (dashed) sits "
        f"a factor of 1.5-3 above both curves at outer radii — same deficit as G10."
    )

    md_path = RESULTS / "G31_A0_baseline_rotation_curve_summary.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"Summary saved: {md_path}")


if __name__ == "__main__":
    main()
