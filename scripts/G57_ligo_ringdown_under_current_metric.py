#!/usr/bin/env python3
"""
G57_ligo_ringdown_under_current_metric.py

LIGO ringdown prediction under the CURRENT k(A) commitment (2026-05-12 evening),
contrasted with the obsolete G1 prediction.

Current strong-field metric (per README and
results/strong_field_metric_commitment_2026_05_12.md):

    k(A) = (1 - A)                          for A <= 2/3   (outside PS, exact GR)
    k(A) = 27(1 - A)^3 (9A^2 - 10A + 3)     for 2/3 < A < 1 (final-shell quartic Hermite)

Key consequence for LIGO ringdown:

    At the photon sphere (A = 2/3): k(A) = 1 - 2/3 = 1/3 -- EXACTLY Schwarzschild.
    Eikonal QNM at the photon sphere depends on h(r_c) * k(r_c), both of which
    are now GR-identical. Therefore:

        tau_STAM(spinless) / tau_GR(spinless) = 1.000   (eikonal)

    This replaces G1's prediction of 1.80, which used the obsolete
    k(A) = (1-A)(1-A^2)^2 form with k(r_c) = 25/243 ~ 0.103.

Smoothness of new k at photon sphere:
    k and dk/dA and d^2k/dA^2 all match the outside-PS GR form at A = 2/3
    (C^2 at PS, by construction). So first-order WKB correction beyond
    eikonal is ALSO GR-identical. Departures from GR appear only at
    higher-order WKB and from the region A > 2/3 (inside PS).

Spinning case: per HANDOFF_2026_05_12_strong_field.md, the static-bubble +
rotating-matter reading gives tau_STAM / tau_GR_Kerr ~ 0.87 at typical
Kerr remnant spin (a ~ 0.67). That's a 13% deficit, well within current
LIGO precision (+/- 20-30%). Contingent on STAM-Kerr extension being
derived; currently a directional prediction, not a sharp number.
"""

from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Constants
C = 2.99792458e8                 # m/s
G = 6.67430e-11                  # m^3/(kg s^2)
M_SUN = 1.98892e30               # kg

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

A_PS = 2.0 / 3.0  # photon sphere
A_HORIZON = 1.0


# --- metric components in geometric units (c = G = 1, M = 1) ---

def A_of_r(r, M=1.0):
    return 2.0 * M / r


def h(r, M=1.0):
    """g_tt coefficient. h(A) = 1 - A, identical for Schwarzschild and STAM."""
    return 1.0 - A_of_r(r, M)


def k_schwarzschild(r, M=1.0):
    return 1.0 - A_of_r(r, M)


def k_stam_old(r, M=1.0):
    """Obsolete (pre-2026-05-12 evening): k(A) = (1-A)(1-A^2)^2.
    Kept for direct comparison with G1.
    """
    A = A_of_r(r, M)
    return (1.0 - A) * (1.0 - A**2) ** 2


def k_stam_current_A(A):
    """Current k as a function of A directly."""
    if A <= A_PS:
        return 1.0 - A
    return 27.0 * (1.0 - A) ** 3 * (9.0 * A**2 - 10.0 * A + 3.0)


def k_stam_current(r, M=1.0):
    """Current commitment (2026-05-12 evening): piecewise.
        k = (1-A)                        for A <= 2/3   (outside PS, exact GR)
        k = 27(1-A)^3 (9A^2 - 10A + 3)   for A > 2/3    (final-shell quartic Hermite)
    """
    return k_stam_current_A(A_of_r(r, M))


# --- photon-sphere quantities (h-only, so identical for all three) ---

def photon_sphere_radius(M=1.0):
    return 3.0 * M


def Omega_c(M=1.0):
    r_c = photon_sphere_radius(M)
    return np.sqrt(h(r_c, M)) / r_c


def lambda_lyapunov(k_func, M=1.0):
    """At the photon sphere of h = 1 - 2M/r (r_c = 3M):
        |V''_eff(r_c)| / E^2 = 2 / (3 M^2)
        lambda^2 = (h_c * k_c / 2) * (2 / (3 M^2)) = (h_c * k_c) / (3 M^2)
    """
    r_c = photon_sphere_radius(M)
    h_c = h(r_c, M)
    k_c = k_func(r_c, M)
    V_dd_over_E2 = 2.0 / (3.0 * M**2)
    return np.sqrt((h_c * k_c / 2.0) * V_dd_over_E2)


# --- QNM eikonal ---

def qnm_eikonal(k_func, l, n, M=1.0):
    omega_R = l * Omega_c(M)
    lam = lambda_lyapunov(k_func, M)
    omega_I = -(n + 0.5) * lam
    return {
        "omega_R": omega_R,
        "omega_I": omega_I,
        "lambda": lam,
        "Q": omega_R / (2.0 * abs(omega_I)),
    }


def to_SI(omega_geom, M_solar):
    M_kg = M_solar * M_SUN
    geom_M_in_seconds = G * M_kg / C**3
    return omega_geom / geom_M_in_seconds


def ringdown_obs(k_func, M_solar, l=2, n=0):
    q = qnm_eikonal(k_func, l, n)
    f_Hz = to_SI(q["omega_R"], M_solar) / (2.0 * np.pi)
    tau_s = 1.0 / abs(to_SI(q["omega_I"], M_solar))
    return {
        "f_Hz": f_Hz,
        "tau_s": tau_s,
        "tau_ms": tau_s * 1e3,
        "Q": q["Q"],
        "omega_R": q["omega_R"],
        "omega_I": q["omega_I"],
        "lambda": q["lambda"],
    }


BENCHMARKS = [
    ("GW150914 remnant", 62.0),
    ("GW170729 remnant", 80.0),
    ("GW190521 remnant", 142.0),
    ("Stellar-mass BBH (typical)", 30.0),
    ("Intermediate-mass (~1000 M_sun)", 1000.0),
    ("Sgr A*", 4.3e6),
    ("M87*", 6.5e9),
]


# --- plots ---

def plot_k_comparison():
    fig, ax = plt.subplots(figsize=(11, 6))
    A_grid = np.linspace(0.001, 0.9995, 2000)

    k_S = 1.0 - A_grid
    k_old = (1.0 - A_grid) * (1.0 - A_grid**2) ** 2
    k_new = np.where(
        A_grid <= A_PS,
        1.0 - A_grid,
        27.0 * (1.0 - A_grid) ** 3 * (9.0 * A_grid**2 - 10.0 * A_grid + 3.0),
    )

    ax.plot(A_grid, k_S, "b-", linewidth=2, label="Schwarzschild  k = 1 - A")
    ax.plot(A_grid, k_old, "r--", linewidth=1.8, alpha=0.75,
            label="STAM (obsolete, G1)  k = (1-A)(1-A^2)^2")
    ax.plot(A_grid, k_new, "g-", linewidth=2.6,
            label=("STAM (current)  k = (1-A) for A<=2/3; "
                   "27(1-A)^3(9A^2-10A+3) for A>2/3"))
    ax.axvline(1.0 / 3.0, color="gray", linestyle=":", alpha=0.4,
               label="ISCO (A=1/3)")
    ax.axvline(A_PS, color="gray", linestyle="--", alpha=0.6,
               label="Photon sphere (A=2/3)")
    ax.set_xlabel("A = R_s / r")
    ax.set_ylabel("k(A)   (g_rr = 1 / k)")
    ax.set_title(
        "k(A) under three metric commitments\n"
        "Current STAM = exact GR outside the photon sphere "
        "(C^2 at PS, C^1 at horizon)"
    )
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.02, 1.02)
    plt.tight_layout()
    out = PLOTS / "G57_k_comparison.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_ringdown_waveforms():
    M_solar = 62.0
    r_S = ringdown_obs(k_schwarzschild, M_solar)
    r_old = ringdown_obs(k_stam_old, M_solar)
    r_new = ringdown_obs(k_stam_current, M_solar)

    t_ms = np.linspace(0, 30, 5000)
    t_s = t_ms * 1e-3

    h_S = np.exp(-t_s / r_S["tau_s"]) * np.cos(2 * np.pi * r_S["f_Hz"] * t_s)
    h_old = np.exp(-t_s / r_old["tau_s"]) * np.cos(2 * np.pi * r_old["f_Hz"] * t_s)
    h_new = np.exp(-t_s / r_new["tau_s"]) * np.cos(2 * np.pi * r_new["f_Hz"] * t_s)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    axes[0].plot(t_ms, h_S, "b-", linewidth=1.5,
                 label=f"Schwarzschild  tau = {r_S['tau_ms']:.2f} ms")
    axes[0].plot(t_ms, h_old, "r--", linewidth=1.2, alpha=0.7,
                 label=f"STAM obsolete (G1)  tau = {r_old['tau_ms']:.2f} ms")
    axes[0].plot(t_ms, h_new, "g-", linewidth=1.5, alpha=0.85,
                 label=f"STAM current  tau = {r_new['tau_ms']:.2f} ms")
    axes[0].axhline(0, color="black", linewidth=0.4)
    axes[0].set_ylabel("Ringdown amplitude (normalized)")
    axes[0].set_title(
        f"GW150914-like ringdown ({M_solar} M_sun remnant)\n"
        "Current STAM ringdown = GR exactly (spinless eikonal)"
    )
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    env_S = np.exp(-t_s / r_S["tau_s"])
    env_old = np.exp(-t_s / r_old["tau_s"])
    env_new = np.exp(-t_s / r_new["tau_s"])
    axes[1].semilogy(t_ms, env_S, "b-", linewidth=2, label="Schwarzschild envelope")
    axes[1].semilogy(t_ms, env_old, "r--", linewidth=1.6,
                     label="STAM obsolete envelope (G1)")
    axes[1].semilogy(t_ms, env_new, "g-", linewidth=2,
                     label="STAM current envelope")
    axes[1].set_xlabel("Time after merger (ms)")
    axes[1].set_ylabel("|envelope|")
    axes[1].legend()
    axes[1].grid(True, which="both", alpha=0.3)
    axes[1].set_ylim(1e-6, 2)

    plt.tight_layout()
    out = PLOTS / "G57_ringdown_waveforms.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_inside_ps_departure():
    """Inside PS: show how the current k drops relative to GR, and where
    higher-overtone / late-inspiral observables could pick up the difference.
    """
    fig, ax = plt.subplots(figsize=(11, 6))
    A_inside = np.linspace(A_PS, 0.999, 500)
    k_GR = 1.0 - A_inside
    k_new = 27.0 * (1.0 - A_inside) ** 3 * (9.0 * A_inside**2 - 10.0 * A_inside + 3.0)
    k_old = (1.0 - A_inside) * (1.0 - A_inside**2) ** 2

    ax.plot(A_inside, k_GR, "b-", linewidth=2, label="GR  k = 1 - A")
    ax.plot(A_inside, k_old, "r--", linewidth=1.6, alpha=0.7,
            label="STAM obsolete (G1)")
    ax.plot(A_inside, k_new, "g-", linewidth=2.2,
            label="STAM current (quartic Hermite)")
    ax.axvline(A_PS, color="gray", linestyle="--", alpha=0.7, label="PS (A=2/3)")
    ax.axvline(A_HORIZON, color="black", linestyle=":", alpha=0.7,
               label="Horizon (A=1)")
    ax.set_xlabel("A inside photon sphere")
    ax.set_ylabel("k(A)")
    ax.set_title("Inside-PS metric departure from GR (current vs obsolete)\n"
                 "Probed by higher overtones, late-inspiral chirp, LISA EMRI")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(A_PS - 0.01, 1.0)
    plt.tight_layout()
    out = PLOTS / "G57_inside_ps_departure.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- markdown summary ---

def write_md(ratio_old, ratio_new, rows_print, k_at_ps, lambda_vals):
    md = []
    md.append("# G57 — LIGO Ringdown Under Current k(A) Commitment\n")
    md.append("**Date: 2026-05-12 (continuation).** Test re-run under the strong-field "
              "metric committed in [strong_field_metric_commitment_2026_05_12.md]"
              "(strong_field_metric_commitment_2026_05_12.md) and described in the "
              "current [README](../README.md). Replaces [G1](../scripts/G1_qnm_ringdown.py) "
              "for the spinless QNM ringdown prediction.\n")

    md.append("## What changed\n")
    md.append("G1 was written when the framework's strong-field metric was "
              "`k(A) = (1-A)(1-A^2)^2`. Under the current commitment:\n")
    md.append("```text")
    md.append("k(A) = (1 - A)                          for A <= 2/3       (outside PS, exact GR)")
    md.append("k(A) = 27(1-A)^3 (9A^2 - 10A + 3)       for 2/3 < A < 1    (final-shell quartic Hermite)")
    md.append("```\n")
    md.append("**At the photon sphere (A = 2/3):**")
    md.append("")
    md.append("| Metric                | k(A_PS)              |")
    md.append("|---|---:|")
    md.append(f"| Schwarzschild         | 1/3 = {k_at_ps['S']:.6f} |")
    md.append(f"| STAM obsolete (G1)    | 25/243 = {k_at_ps['old']:.6f} |")
    md.append(f"| **STAM current**      | **1/3 = {k_at_ps['new']:.6f}**  *(= GR exactly)* |")
    md.append("")
    md.append("Because both `h` and `k` are now identical to Schwarzschild at the photon "
              "sphere, the eikonal Lyapunov exponent `lambda^2 = h*k / (3 M^2)` is also "
              "identical to Schwarzschild's:\n")
    md.append("```text")
    md.append(f"lambda_Schwarzschild  = {lambda_vals['S']:.6f} / M")
    md.append(f"lambda_STAM_obsolete  = {lambda_vals['old']:.6f} / M   (was 5/9 smaller -> hence the historic 1.80)")
    md.append(f"lambda_STAM_current   = {lambda_vals['new']:.6f} / M   (= Schwarzschild)")
    md.append("```\n")

    md.append("## Spinless eikonal prediction\n")
    md.append("| Metric             | tau_STAM / tau_GR (spinless, eikonal) |")
    md.append("|---|---:|")
    md.append(f"| G1 (obsolete)       | **{ratio_old:.4f}**   (the historic 9/5 = 1.80) |")
    md.append(f"| Current             | **{ratio_new:.4f}**   (exact GR) |")
    md.append("")
    md.append("The framework's earlier signature prediction (1.80x longer ringdown) is "
              "**retired** under the current k(A). Spinless STAM produces the Schwarzschild "
              "eikonal QNM exactly because the metric is exactly GR at and outside the "
              "photon sphere.\n")
    md.append("More strongly: the quartic Hermite F(Sigma) was constructed with `F(2) = 1, "
              "F'(2) = 0, F''(2) = 0` (C^2 at PS), so `k`, `dk/dA`, and `d^2k/dA^2` all "
              "agree with the outside-PS GR form at A = 2/3. First-order WKB beyond "
              "eikonal is also GR-identical at the photon sphere. Sub-leading departures "
              "require higher-order WKB corrections or the inside-PS region.\n")

    md.append("## Predictions for representative LIGO events (spinless eikonal)\n")
    md.append("```text")
    md.append(f"{'Event':<35s}  {'M (M_sun)':>10s}  "
              f"{'f (Hz)':>10s}  {'tau_GR (ms)':>12s}  "
              f"{'tau_old (ms)':>13s}  {'tau_new (ms)':>13s}")
    md.append("-" * 100)
    for row in rows_print:
        md.append(
            f"{row['label']:<35s}  {row['M']:>10.2g}  "
            f"{row['f_S']:>10.1f}  {row['tau_S']:>12.4g}  "
            f"{row['tau_old']:>13.4g}  {row['tau_new']:>13.4g}"
        )
    md.append("```")
    md.append("Frequencies are identical across all three metrics because they depend only "
              "on `h`, which is unchanged. Damping time tau is what the metric modification "
              "of g_rr would have shifted; under the current k, it doesn't.\n")

    md.append("## Where the current k still differs from GR (inside the photon sphere)\n")
    md.append("Inside the photon sphere (A > 2/3), the current k drops sharply and reaches "
              "zero with order D = 3 at the horizon. Sampling:\n")
    md.append("")
    md.append("| A | k_GR | k_STAM(current) | ratio |")
    md.append("|---:|---:|---:|---:|")
    for A in [0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99]:
        kG = 1.0 - A
        kN = 27.0 * (1.0 - A) ** 3 * (9.0 * A**2 - 10.0 * A + 3.0)
        md.append(f"| {A:.2f} | {kG:.4f} | {kN:.4e} | {kN/kG:.4f} |")
    md.append("")
    md.append("Eikonal QNM is photon-sphere-localized and doesn't see this. What does see it:\n")
    md.append("- **Higher overtones (n >= 1)**: increasingly sensitive to the geometry just "
              "inside r_c. High-SNR overtone extraction in LIGO O5+ could probe this directly.\n")
    md.append("- **Higher-order WKB / full Regge-Wheeler**: integrates over the radial "
              "profile, picking up the inside-PS region. Needed for sharp sub-percent "
              "predictions on n=0 itself.\n")
    md.append("- **Late-inspiral chirp shape**: probes A close to 2/3 from the merger side.\n")
    md.append("- **LISA EMRI ringdowns**: probe intermediate A with high precision.\n")

    md.append("## Spinning case (LIGO observed remnants)\n")
    md.append("LIGO BBH remnants are spinning Kerr-like, typically a ~ 0.5-0.7. Per "
              "[HANDOFF_2026_05_12_strong_field.md](../memory/HANDOFF_2026_05_12_strong_field.md), "
              "the static-bubble + rotating-matter ontology gives:\n")
    md.append("```text")
    md.append("tau_STAM / tau_GR_Kerr  ~  0.87   at a ~ 0.67   (13% deficit)")
    md.append("```\n")
    md.append("This is within current LIGO ringdown precision (+/- 20-30%). It is a "
              "**directional prediction**, not a sharp number -- the explicit STAM-Kerr "
              "metric has not been derived. Two readings of the spinning extension (naive "
              "A_Kerr_BL vs static-bubble) currently give different predictions. Closing "
              "that gap is the next strong-field deliverable.\n")

    md.append("## Status summary\n")
    md.append(f"- **Spinless eikonal**: tau_STAM / tau_GR = {ratio_new:.4f} (was {ratio_old:.4f} under G1)\n")
    md.append("- **First-order WKB beyond eikonal**: also exact GR (C^2 smoothness at PS)\n")
    md.append("- **Spinning Kerr ringdown**: ~ 13% deficit estimate, contingent on STAM-Kerr extension\n")
    md.append("- **G1 result (1.80) is preserved** in [G1_qnm_ringdown.py](../scripts/G1_qnm_ringdown.py) "
              "as the historical record of the obsolete k(A) form\n")
    md.append("- **G53 candidate-n exploration** is also superseded; current k is not in the "
              "(1-A)(1-A^2)^n family but a piecewise form with explicit outside-PS = GR\n")

    md.append("## Implications\n")
    md.append("1. **Spinless LIGO consistency is now automatic.** The framework no longer "
              "has a sharp tau-ratio tension to defend against current LIGO ringdown bounds.\n")
    md.append("2. **The strong-field STAM-vs-GR wedge is pushed inside the photon sphere.** "
              "All photon-sphere-localized observables (n=0 ringdown frequency and damping "
              "in eikonal, EHT shadow, light bending up to the PS) are now GR-exact in "
              "the spinless case. STAM departures live inside r_c.\n")
    md.append("3. **The framework's main strong-field empirical handle is now Kerr.** "
              "Spinless ringdown can no longer distinguish STAM from GR; spinning ringdown "
              "(via the static-bubble + rotating-matter reading) is the testable channel. "
              "STAM-Kerr extension is gating the framework's near-term LIGO comparison.\n")

    md.append("## Files\n")
    md.append("- [scripts/G57_ligo_ringdown_under_current_metric.py](../scripts/G57_ligo_ringdown_under_current_metric.py)\n")
    md.append("- [plots/G57_k_comparison.png](../plots/G57_k_comparison.png)\n")
    md.append("- [plots/G57_ringdown_waveforms.png](../plots/G57_ringdown_waveforms.png)\n")
    md.append("- [plots/G57_inside_ps_departure.png](../plots/G57_inside_ps_departure.png)\n")

    out = RESULTS / "G57_ligo_ringdown_under_current_metric_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main():
    print("=" * 80)
    print("G57: LIGO ringdown under current k(A) commitment")
    print("=" * 80)
    print()

    r_c = photon_sphere_radius(1.0)
    k_S = k_schwarzschild(r_c)
    k_old = k_stam_old(r_c)
    k_new = k_stam_current(r_c)

    print("Photon-sphere k(A) at A = 2/3:")
    print(f"  Schwarzschild       k = 1/3      = {k_S:.6f}")
    print(f"  STAM obsolete (G1)  k = 25/243   = {k_old:.6f}")
    print(f"  STAM current        k = 1/3      = {k_new:.6f}  <-- exact GR")
    print()

    lam_S = lambda_lyapunov(k_schwarzschild)
    lam_old = lambda_lyapunov(k_stam_old)
    lam_new = lambda_lyapunov(k_stam_current)
    print("Lyapunov exponent at photon sphere (per M):")
    print(f"  lambda_Schwarzschild = {lam_S:.6f}")
    print(f"  lambda_STAM_old      = {lam_old:.6f}")
    print(f"  lambda_STAM_current  = {lam_new:.6f}")
    print()

    ratio_old = lam_S / lam_old
    ratio_new = lam_S / lam_new
    print("Eikonal tau_STAM / tau_GR (spinless):")
    print(f"  G1 obsolete:    tau_STAM / tau_GR = {ratio_old:.4f}   (= 9/5 = 1.80)")
    print(f"  Current:        tau_STAM / tau_GR = {ratio_new:.4f}   (exact GR)")
    print()

    print("Spinless benchmark table (eikonal):")
    print(f"  {'Event':<35} {'M(M_sun)':>10} {'f(Hz)':>10} "
          f"{'tau_GR(ms)':>12} {'tau_old(ms)':>13} {'tau_new(ms)':>13}")
    print("  " + "-" * 96)
    rows_print = []
    for label, M in BENCHMARKS:
        rS = ringdown_obs(k_schwarzschild, M)
        rO = ringdown_obs(k_stam_old, M)
        rN = ringdown_obs(k_stam_current, M)
        print(f"  {label:<35} {M:>10.2g} {rS['f_Hz']:>10.1f} "
              f"{rS['tau_ms']:>12.4g} {rO['tau_ms']:>13.4g} "
              f"{rN['tau_ms']:>13.4g}")
        rows_print.append({
            "label": label, "M": M, "f_S": rS["f_Hz"],
            "tau_S": rS["tau_ms"], "tau_old": rO["tau_ms"],
            "tau_new": rN["tau_ms"],
        })
    print()
    print("Note: frequencies are identical across all three metrics (h is unchanged).")
    print()

    print("Inside the photon sphere (A > 2/3), current k drops sharply:")
    print(f"  {'A':>6} {'k_GR':>12} {'k_STAM(new)':>16} {'ratio':>10}")
    for A in [0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99]:
        kG = 1.0 - A
        kN = 27.0 * (1.0 - A) ** 3 * (9.0 * A**2 - 10.0 * A + 3.0)
        print(f"  {A:>6.2f} {kG:>12.6f} {kN:>16.6e} {kN/kG:>10.4f}")
    print()
    print("Eikonal QNM doesn't see this region; higher overtones, late-inspiral chirp,")
    print("and LISA EMRIs would.")
    print()

    # Plots
    print("Generating plots...")
    p1 = plot_k_comparison()
    p2 = plot_ringdown_waveforms()
    p3 = plot_inside_ps_departure()
    print(f"  {p1}")
    print(f"  {p2}")
    print(f"  {p3}")
    print()

    md = write_md(
        ratio_old, ratio_new, rows_print,
        k_at_ps={"S": k_S, "old": k_old, "new": k_new},
        lambda_vals={"S": lam_S, "old": lam_old, "new": lam_new},
    )
    print(f"Summary: {md}")
    print()
    print("Bottom line:")
    print("  Old prediction (G1, obsolete):   tau_STAM / tau_GR = 1.80")
    print("  Current prediction (G57):        tau_STAM / tau_GR = 1.000 (spinless eikonal)")
    print("  Spinning Kerr (per handoff):     tau_STAM / tau_GR_Kerr ~ 0.87 (contingent on STAM-Kerr)")


if __name__ == "__main__":
    main()
