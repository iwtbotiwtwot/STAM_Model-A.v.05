#!/usr/bin/env python3
"""
F5c_gw170817_arrival_test.py

Numerical test of Model-B's strong claim — "GW propagates at c always;
light is slowed by (1/c) * integral of A along the path" — against the
GW170817 multi-messenger observation: gamma-ray counterpart arrived
~1.74 s after the GW signal from a source ~40 Mpc away in NGC 4993.

Author: Sean Brady / STAM Model-B continuation 2026-05-08

Setup:
    Model-B (per F5b summary):
        Delta_t_GW    = D / c                        (GW intrinsically at c)
        Delta_t_light = D / c + (1/c) * int A(s) ds  (light Shapiro-slowed by A)
    Differential:
        Delta_t = (1/c) * int A(s) ds
    For matter-sourced A (Framework C: nabla^2 A = (8 pi G / c^2) rho_matter)
    a point mass M sources A(r) = 2 G M / (c^2 r), so each mass along the
    line of sight contributes a Shapiro-style integral.

What this script does:
    1. Estimate int A(s) ds along the GW170817 photon path using:
        - Milky Way: M_MW = 1e12 M_sun at galactic center,
          impact parameter b = R_sun_GC * sin(70 deg)
          -> integral = (2GM/c^2) * asinh(D/b)
        - NGC 4993: M = 4e10 M_sun, kilonova offset 2 kpc,
          radial-outward path -> integral = (2GM/c^2) * ln(D/r_off)
    2. Compute predicted delta_t = int(A) / c.
    3. Compare to 1.74 s observed.
    4. Report what average A would be required to MATCH observation,
       and contrast with realistic A values along the LOS.

Honest expectation:
    The Milky Way alone contributes ~10^16 m of integrated A, predicting a
    GW-vs-light delay of ~3 years — about 8 orders of magnitude larger
    than the observed 1.74 s. The NGC 4993 contribution adds another
    ~45 days. Strong-form Model-B (GW skips local A entirely; light feels
    full matter-sourced A) is INCONSISTENT with GW170817 unless A is
    redefined.

Why the test matters:
    In standard GR both GW and light experience the same Shapiro delay,
    so the differential is zero (the observed 1.74 s is attributed to jet
    launch / breakout astrophysics). Any theory that gives the messengers
    different propagation speeds through galactic potentials must answer
    to GW170817's tight multi-messenger arrival window.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# --- constants (SI) ---
C = 2.99792458e8
G = 6.67430e-11
M_SUN = 1.98892e30
PC = 3.0857e16
KPC = 1.0e3 * PC
MPC = 1.0e6 * PC
SECONDS_PER_DAY = 86400.0
SECONDS_PER_YEAR = 365.25 * SECONDS_PER_DAY

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

# --- GW170817 observed ---
DT_OBSERVED = 1.74                # s, gamma-ray vs GW arrival
D_SOURCE = 40.0 * MPC             # ~40 Mpc to NGC 4993

# --- Milky Way ---
M_MW = 1.0e12 * M_SUN             # total mass (incl. dark halo)
R_SUN_GAL = 8.0 * KPC
GC_ANGLE_DEG = 70.0               # angular separation Sun-GC vs Sun-NGC4993
B_GC = R_SUN_GAL * np.sin(np.deg2rad(GC_ANGLE_DEG))   # impact parameter to GC

# --- NGC 4993 ---
M_NGC = 4.0e10 * M_SUN            # stellar + DM mass (rough)
R_KILONOVA = 2.0 * KPC            # kilonova offset from NGC 4993 center


# --- A field integrals along line of sight ---

def integrate_A_milky_way(D: float = D_SOURCE) -> dict:
    """Integral of A_MW along photon path with impact parameter b to GC.

    A_MW(s) = 2GM/(c^2 * sqrt(b^2 + s^2))  (point-mass approximation)
    integral_0^D A ds = (2GM/c^2) * asinh(D/b).
    """
    coeff = 2.0 * G * M_MW / (C ** 2)
    integral = coeff * np.arcsinh(D / B_GC)
    return {
        "label": "Milky Way",
        "integral_m": float(integral),
        "delay_s": float(integral / C),
    }


def integrate_A_ngc(D: float = D_SOURCE) -> dict:
    """Integral of A_NGC for a photon emerging radially from a kilonova at
    offset r_off from NGC 4993 center, integrated outward to source distance.

    integral_{r_off}^D (2GM/(c^2 r)) dr = (2GM/c^2) * ln(D / r_off).
    """
    coeff = 2.0 * G * M_NGC / (C ** 2)
    integral = coeff * np.log(D / R_KILONOVA)
    return {
        "label": "NGC 4993",
        "integral_m": float(integral),
        "delay_s": float(integral / C),
    }


def total_predicted() -> dict:
    mw = integrate_A_milky_way()
    ngc = integrate_A_ngc()
    total_int = mw["integral_m"] + ngc["integral_m"]
    total_delay = total_int / C
    A_required = DT_OBSERVED * C / D_SOURCE
    A_mw_local = 2 * G * M_MW / (C ** 2 * R_SUN_GAL)
    return {
        "mw": mw,
        "ngc": ngc,
        "total_integral_m": total_int,
        "total_delay_s": total_delay,
        "ratio_predicted_to_observed": total_delay / DT_OBSERVED,
        "A_required_for_observed_match": A_required,
        "A_mw_at_sun_for_reference": A_mw_local,
    }


# --- plots ---

def plot_delay_comparison(result: dict) -> Path:
    fig, ax = plt.subplots(figsize=(10, 6))
    labels = [
        "Observed\n(GW170817)",
        "Model-B prediction\nMW only",
        "Model-B prediction\nNGC 4993 only",
        "Model-B prediction\nTOTAL",
    ]
    delays = [
        DT_OBSERVED,
        result["mw"]["delay_s"],
        result["ngc"]["delay_s"],
        result["total_delay_s"],
    ]
    colors = ["tab:green", "tab:orange", "tab:purple", "tab:red"]

    bars = ax.bar(labels, delays, color=colors, edgecolor="black")
    ax.set_yscale("log")
    ax.set_ylabel("GW-vs-light arrival delay (s)")
    ax.set_title("GW170817 arrival delay: observed vs Model-B prediction")
    ax.grid(True, axis="y", alpha=0.3, which="both")
    ax.axhline(DT_OBSERVED, color="tab:green", linestyle="--", alpha=0.6,
               label=f"Observed: {DT_OBSERVED} s")

    for bar, d in zip(bars, delays):
        if d > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, d * 1.6, f"{d:.2e} s",
                    ha="center", fontsize=9)

    ax.legend()
    out = PLOTS / "F5c_gw170817_arrival_delay.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_A_along_path(result: dict) -> Path:
    """A(s) along the line of sight, observer (s=0) to source (s=D)."""
    s_mpc = np.logspace(-3, np.log10(D_SOURCE / MPC), 5000)
    s = s_mpc * MPC

    A_mw = 2.0 * G * M_MW / (C ** 2 * np.sqrt(B_GC ** 2 + s ** 2))
    dist_from_ngc = np.maximum(D_SOURCE - s, R_KILONOVA)
    A_ngc = 2.0 * G * M_NGC / (C ** 2 * dist_from_ngc)
    A_total = A_mw + A_ngc

    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.loglog(s_mpc, A_mw, color="tab:orange", label="A from Milky Way")
    ax.loglog(s_mpc, A_ngc, color="tab:purple", label="A from NGC 4993")
    ax.loglog(s_mpc, A_total, color="tab:red", linewidth=2, label="Total A(s)")
    ax.axhline(result["A_required_for_observed_match"], color="tab:green",
               linestyle="--", linewidth=1.5,
               label=("Average A required for 1.74 s observed: "
                      f"{result['A_required_for_observed_match']:.2e}"))
    ax.set_xlabel("Distance from Earth along line of sight (Mpc)")
    ax.set_ylabel("A(s)  [dimensionless]")
    ax.set_title("Model-B's matter-sourced A field along GW170817 photon path")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)

    out = PLOTS / "F5c_A_along_path.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- markdown summary ---

def write_markdown(result: dict, plot_paths: list[Path]) -> Path:
    md = []
    md.append("# F5c: GW170817 Multi-Messenger Test of Model-B's "
              "'GW = c, light slowed by A' Claim\n")

    md.append("## The question\n")
    md.append(
        "F5b's summary asserts that GW170817's ~1.74 s arrival offset "
        "between the GW signal and the gamma-ray counterpart is consistent "
        "with Model-B **by construction** — GWs propagate at c intrinsically, "
        "and light is slowed only by `(1/c) * integral A ds` along the path. "
        "This script tests that claim numerically: take the strong form of "
        "the claim at face value, compute the predicted delay over the "
        "40 Mpc path to NGC 4993, and compare to 1.74 s.\n"
    )

    md.append("## Setup\n")
    md.append(
        "Model-B propagation:\n"
        "```text\n"
        "Delta_t_GW    = D / c                          (no A coupling)\n"
        "Delta_t_light = D / c + (1/c) * int A(s) ds    (Shapiro-slowed)\n"
        "Differential  = (1/c) * int A(s) ds\n"
        "```\n"
        "\n"
        "A is matter-sourced (Framework C), so each mass along the line of "
        "sight contributes A(r) = 2GM/(c^2 r). Two dominant contributions:\n"
        "\n"
        f"- **Milky Way**: M = {M_MW/M_SUN:.1e} M_sun at the galactic center; "
        f"photon path with impact parameter b = R_sun * sin({GC_ANGLE_DEG:.0f} deg) "
        f"= {B_GC/KPC:.2f} kpc.\n"
        f"  Integral closed form: (2GM/c^2) * asinh(D/b).\n"
        "\n"
        f"- **NGC 4993**: M = {M_NGC/M_SUN:.1e} M_sun; kilonova at "
        f"{R_KILONOVA/KPC:.1f} kpc offset, photon traveling radially outward.\n"
        "  Integral closed form: (2GM/c^2) * ln(D/r_off).\n"
        "\n"
        "Both integrals are logarithmically divergent for an isolated point "
        "mass; we cut off at the source distance D = 40 Mpc.\n"
    )

    md.append("## Results\n")
    md.append(
        "```text\n"
        f"MW   integral A ds = {result['mw']['integral_m']:.3e} m   "
        f"-> Delta_t = {result['mw']['delay_s']:.3e} s  "
        f"({result['mw']['delay_s']/SECONDS_PER_YEAR:.2f} yr)\n"
        f"NGC  integral A ds = {result['ngc']['integral_m']:.3e} m   "
        f"-> Delta_t = {result['ngc']['delay_s']:.3e} s  "
        f"({result['ngc']['delay_s']/SECONDS_PER_DAY:.1f} days)\n"
        f"TOTAL              = {result['total_integral_m']:.3e} m   "
        f"-> Delta_t = {result['total_delay_s']:.3e} s  "
        f"({result['total_delay_s']/SECONDS_PER_YEAR:.2f} yr)\n"
        "```\n"
        "\n"
        f"Observed GW170817 delay: **{DT_OBSERVED} s**.\n"
        "\n"
        f"Ratio predicted / observed: "
        f"**{result['ratio_predicted_to_observed']:.2e}**.\n"
    )

    md.append("## What value of A would match the observation?\n")
    md.append(
        "Inverting `Delta_t = (1/c) * <A> * D`:\n"
        "```text\n"
        f"<A>_required = c * Delta_t / D = "
        f"{result['A_required_for_observed_match']:.3e}\n"
        f"A_MW at the Sun (reference)   = "
        f"{result['A_mw_at_sun_for_reference']:.3e}\n"
        "```\n"
        "\n"
        f"The path-average A required to reproduce 1.74 s is about "
        f"{result['A_mw_at_sun_for_reference'] / result['A_required_for_observed_match']:.1e} "
        "times **smaller** than the local Milky Way potential A at the Sun. "
        "For a photon to traverse the Galaxy and have an effective "
        "path-average A this small, light must NOT couple to the local "
        "matter-sourced A as written.\n"
    )

    md.append("## Verdict\n")
    md.append(
        "**The strong-form Model-B claim — 'GW always c; light slowed by "
        "the full matter-sourced A integral' — predicts a GW-vs-light "
        "arrival offset of ~years for GW170817. Observed: 1.74 s. The "
        "model is INCONSISTENT with the data by roughly eight orders of "
        "magnitude.**\n"
        "\n"
        "The reason standard GR matches GW170817 cleanly is that BOTH "
        "messengers experience identical Shapiro delay through galactic "
        "potentials, so the differential cancels and only the "
        "astrophysical jet-launch ~1.7 s remains. Any framework that "
        "decouples GW from the gravitational potential while keeping "
        "light coupled will, generically, predict large differential "
        "delays from local potentials.\n"
        "\n"
        "F5b's 'passes by construction' claim does not survive a numerical "
        "evaluation under the construction it specifies.\n"
    )

    md.append("## What Model-B has to do to actually pass GW170817\n")
    md.append(
        "Three structural options, each with consequences:\n"
        "\n"
        "**Option A — A is *not* the full Newtonian potential.** Restrict "
        "A to a cosmological background (or a mode that decouples from "
        "local virialized matter), so light traveling through galactic "
        "potentials sees ordinary GR Shapiro delay (felt equally by GW), "
        "and only a tiny extra A-induced delay from cosmological "
        "background remains. **Cost:** breaks the README's local A "
        "results — GPS clock-rate offset, solar-system Shapiro delay, "
        "the gravity bridge `g = (c^2/2) grad A` — all of which assume "
        "A = 2GM/(c^2 r) for local masses.\n"
        "\n"
        "**Option B — GW also feels A, identically to light.** Restore "
        "messenger equivalence; the differential delay collapses to the "
        "astrophysical jet-launch ~1.7 s. **Cost:** Model-B's "
        "narrative claim that 'GW = space, light = matter on space' "
        "becomes purely interpretive — there is no observable predictive "
        "distinction between the two messengers' propagation.\n"
        "\n"
        "**Option C — A as a perturbation, not the full potential.** "
        "Define A as the deviation from a fiducial GR background, so "
        "local matter still sources the GR metric (and both messengers "
        "feel it equally), with A capturing only an additional "
        "small-amplitude scalar mode. **Cost:** real theoretical work — "
        "specifying the coupling structure is essentially the open "
        "problem flagged in the README ('Deriving the full field "
        "equations').\n"
        "\n"
        "Until one of these is committed and the local A-based "
        "calculations re-derived consistently, GW170817's 1.74 s "
        "remains an unresolved tension for Model-B as currently stated.\n"
    )

    md.append("## Caveats and assumptions\n")
    md.append(
        "- Milky Way modeled as a point mass at the galactic center. "
        "A realistic NFW halo would change the prefactor by O(1) but "
        "not the order of magnitude.\n"
        "- NGC 4993 modeled as a point mass at the galaxy center. The "
        "kilonova offset of 2 kpc is consistent with observed offsets "
        "of short-GRB hosts; varying it by an order of magnitude "
        "moves the NGC delay by a factor of ~ln(10) = 2.3, not enough "
        "to change the verdict.\n"
        "- Intergalactic structure (Local Group, Virgo Cluster, "
        "intervening haloes) is neglected. Including them would only "
        "increase the predicted delay.\n"
        "- The line-of-sight geometry uses a single representative angle "
        "(70 deg from Sun-GC). Realistic geometry gives a similar "
        "impact parameter to within a factor of 2.\n"
        "- Cosmological expansion (proper vs comoving distance over "
        "40 Mpc) is a few-percent correction at this redshift and "
        "doesn't affect the conclusion.\n"
    )

    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "F5c_gw170817_arrival_test_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    result = total_predicted()
    plots = [plot_delay_comparison(result), plot_A_along_path(result)]
    summary = write_markdown(result, plots)

    print("F5c: GW170817 multi-messenger test of Model-B")
    print("=" * 64)
    print()
    print("Geometry:")
    print(f"  Source distance:     D = {D_SOURCE/MPC:.1f} Mpc")
    print(f"  Milky Way mass:      M_MW = {M_MW/M_SUN:.1e} M_sun")
    print(f"  GC impact parameter: b   = {B_GC/KPC:.2f} kpc")
    print(f"  NGC 4993 mass:       M_NGC = {M_NGC/M_SUN:.1e} M_sun")
    print(f"  Kilonova offset:     r_off = {R_KILONOVA/KPC:.1f} kpc")
    print()
    print("Line-of-sight integral A ds:")
    print(f"  MW:    {result['mw']['integral_m']:.3e} m")
    print(f"  NGC:   {result['ngc']['integral_m']:.3e} m")
    print(f"  total: {result['total_integral_m']:.3e} m")
    print()
    print("Predicted Model-B GW-vs-light delay = (integral A ds) / c:")
    print(f"  MW:    {result['mw']['delay_s']:.3e} s "
          f"({result['mw']['delay_s']/SECONDS_PER_YEAR:.2f} yr)")
    print(f"  NGC:   {result['ngc']['delay_s']:.3e} s "
          f"({result['ngc']['delay_s']/SECONDS_PER_DAY:.1f} days)")
    print(f"  total: {result['total_delay_s']:.3e} s "
          f"({result['total_delay_s']/SECONDS_PER_YEAR:.2f} yr)")
    print()
    print(f"Observed GW170817 delay: {DT_OBSERVED:.2f} s")
    print(f"Ratio predicted / observed: "
          f"{result['ratio_predicted_to_observed']:.2e}")
    print()
    print("Average A required to reproduce 1.74 s:")
    print(f"  <A>_required = {result['A_required_for_observed_match']:.3e}")
    print(f"  A_MW at Sun  = {result['A_mw_at_sun_for_reference']:.3e}  "
          "(local-potential reference)")
    print()
    print("VERDICT:")
    print("  Strong-form Model-B (GW=c always; light slowed by full")
    print("  matter-sourced A) is INCONSISTENT with GW170817 by ~8")
    print("  orders of magnitude. The 1.74 s offset is NOT explained")
    print("  as written; F5b's 'passes by construction' assertion needs")
    print("  structural amendment (Options A/B/C in summary).")
    print()
    print("Files written:")
    print(f"  {summary}")
    for p in plots:
        print(f"  {p}")


if __name__ == "__main__":
    main()
