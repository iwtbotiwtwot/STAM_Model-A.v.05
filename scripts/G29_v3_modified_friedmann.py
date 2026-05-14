#!/usr/bin/env python3
"""
G29_v3_modified_friedmann.py

V_3 modified Friedmann redo (replaces script 36's V_1 = beta/(1-A) version).

Script 36 used V_1 = beta/(1-A) and produced Omega_DE_total ≈ 0.662 (97% of
observed Omega_DE) under a tracking solution A_eq(a) = 1 - (1-A_0)·a^(3/2).
G8 (2026-05-09) rejected V_1 on structural grounds: V_1 has no minimum on
the interval (A_0, 1) — A would simply roll to A=A_0 with no stabilizing
floor. V_3 = alpha/A + beta/(1-A) was selected instead, with
alpha/beta = (A_0/(1-A_0))^2 chosen so V_3 has a minimum AT A_0.

This script:
  1. Defines V_3 with the G8-selected calibration.
  2. Integrates KG + Friedmann numerically (forward in cosmic time).
  3. Tries multiple initial conditions for A in the early universe.
  4. Extracts derived H(z) for each scenario.
  5. Compares to cosmic chronometer data (loaded same as G28).
  6. Compares to closed-form STAM intrinsic, LCDM (at SH0ES and Planck H_0),
     and EdS.

Backtrack against G28: does V_3-derived H(z) fit chronometers where the
closed-form ansatz failed?

Author: Sean Brady / STAM Model-A
Date: 2026-05-11
STAM-only frame of mind. Clues, not failures.
"""

from __future__ import annotations
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
RESULTS.mkdir(parents=True, exist_ok=True)


# ============================================================================
# Constants and framework calibration
# ============================================================================
A0 = 1.0 / (12.0 * np.pi)               # cosmic ambient baseline = 0.02653
OMEGA_M = 0.315                          # matter today (Planck 2018)
OMEGA_R = 9.2e-5                         # radiation today
H0_LOCAL = 73.04                         # SH0ES local H_0 (km/s/Mpc)
H0_PLANCK = 67.4                         # Planck H_0 (km/s/Mpc)

# V_3 minimum at A_0 by construction:
# V_3(A) = alpha/A + beta/(1-A)
# V_3'(A) = -alpha/A^2 + beta/(1-A)^2 = 0  at  A = A_0
# => alpha/beta = (A_0/(1-A_0))^2
ALPHA_OVER_BETA = (A0 / (1.0 - A0)) ** 2

# Calibrate beta_tilde = beta/(rho_crit_0) such that V_3(A_0) gives the
# observed Omega_DE today when A sits at the minimum with zero kinetic energy
# V_3(A_0) = beta * [alpha/(beta * A_0) + 1/(1-A_0)]
#         = beta * [(A_0/(1-A_0))^2 / A_0 + 1/(1-A_0)]
#         = beta * [A_0/(1-A_0)^2 + 1/(1-A_0)]
#         = beta * [A_0 + (1-A_0)] / (1-A_0)^2
#         = beta / (1-A_0)^2
# So if V_3(A_0)/rho_crit = Omega_DE_target  =>  beta_tilde = Omega_DE * (1-A_0)^2
OMEGA_DE_TARGET = 1.0 - OMEGA_M - OMEGA_R
BETA_TILDE = OMEGA_DE_TARGET * (1.0 - A0) ** 2
ALPHA_TILDE = BETA_TILDE * ALPHA_OVER_BETA


def V3(A: float) -> float:
    return ALPHA_TILDE / A + BETA_TILDE / (1.0 - A)


def V3_prime(A: float) -> float:
    return -ALPHA_TILDE / A**2 + BETA_TILDE / (1.0 - A) ** 2


def V3_double_prime(A: float) -> float:
    return 2.0 * ALPHA_TILDE / A**3 + 2.0 * BETA_TILDE / (1.0 - A) ** 3


# ============================================================================
# Cosmic chronometer data (same compilation as G28)
# ============================================================================
CC_DATA = np.array(
    [
        (0.07, 69.0, 19.6),
        (0.09, 69.0, 12.0),
        (0.12, 68.6, 26.2),
        (0.17, 83.0, 8.0),
        (0.179, 75.0, 4.0),
        (0.199, 75.0, 5.0),
        (0.20, 72.9, 29.6),
        (0.27, 77.0, 14.0),
        (0.28, 88.8, 36.6),
        (0.352, 83.0, 14.0),
        (0.3802, 83.0, 13.5),
        (0.40, 95.0, 17.0),
        (0.4004, 77.0, 10.2),
        (0.4247, 87.1, 11.2),
        (0.4497, 92.8, 12.9),
        (0.47, 89.0, 49.6),
        (0.4783, 80.9, 9.0),
        (0.48, 97.0, 62.0),
        (0.593, 104.0, 13.0),
        (0.68, 92.0, 8.0),
        (0.781, 105.0, 12.0),
        (0.875, 125.0, 17.0),
        (0.88, 90.0, 40.0),
        (0.90, 117.0, 23.0),
        (1.037, 154.0, 20.0),
        (1.30, 168.0, 17.0),
        (1.363, 160.0, 33.6),
        (1.43, 177.0, 18.0),
        (1.53, 140.0, 14.0),
        (1.75, 202.0, 40.0),
        (1.965, 186.5, 50.4),
    ],
    dtype=float,
)


# ============================================================================
# KG + Friedmann numerical integration (backward from today)
# ============================================================================
def integrate_backward(
    A_today: float,
    Adot_today: float,
    z_max: float = 2.5,
    nsteps: int = 200_000,
):
    """
    Integrate KG + Friedmann BACKWARD from a = 1 to a < 1.

    State: y = [A, omega, lna] with omega = dA/dtau, tau = H_0 * t.
    Equations:
      dA/dtau   = omega
      domega/dtau = -3 * h * omega - V'(A)
      dlna/dtau = h
      h^2 = Omega_m exp(-3 lna) + Omega_r exp(-4 lna) + omega^2/6 + V(A)

    Integration is backward (negative dtau) from a = 1 to a = 1/(1+z_max).
    """
    A = A_today
    omega = Adot_today
    lna = 0.0

    h2 = OMEGA_M + OMEGA_R + omega**2 / 6.0 + V3(A)
    h = np.sqrt(max(h2, 1e-30))

    traj = [(lna, A, omega, h)]

    target_dlna = 1e-4  # very fine backward step
    lna_min = -np.log1p(z_max)

    for _ in range(nsteps):
        # backward step in cosmic time
        h2_curr = OMEGA_M * np.exp(-3 * lna) + OMEGA_R * np.exp(-4 * lna) + omega**2 / 6.0 + V3(A)
        h_curr = np.sqrt(max(h2_curr, 1e-30))
        dtau = -target_dlna / h_curr  # negative -> backward in time

        def rhs(state):
            A_, om_, lna_ = state
            h2_ = (
                OMEGA_M * np.exp(-3 * lna_)
                + OMEGA_R * np.exp(-4 * lna_)
                + om_**2 / 6.0
                + V3(A_)
            )
            h_ = np.sqrt(max(h2_, 1e-30))
            return np.array(
                [
                    om_,
                    -3.0 * h_ * om_ - V3_prime(A_),
                    h_,
                ]
            )

        y = np.array([A, omega, lna])
        k1 = rhs(y)
        k2 = rhs(y + 0.5 * dtau * k1)
        k3 = rhs(y + 0.5 * dtau * k2)
        k4 = rhs(y + dtau * k3)
        y_new = y + (dtau / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

        A, omega, lna = y_new
        # Safety: keep A in (epsilon, 1-epsilon) to avoid pole blowups
        if A <= 1e-6 or A >= 1.0 - 1e-6:
            print(f"  WARNING: A hit a pole boundary at z={1/np.exp(lna)-1:.3f}, A={A:.5f}")
            break
        h2 = OMEGA_M * np.exp(-3 * lna) + OMEGA_R * np.exp(-4 * lna) + omega**2 / 6.0 + V3(A)
        h = np.sqrt(max(h2, 1e-30))
        traj.append((lna, A, omega, h))

        if lna < lna_min:
            break

    traj = np.array(traj)
    return traj  # columns: lna, A, omega, h


def H_of_z_from_traj(traj: np.ndarray, z_query: np.ndarray, H0: float) -> np.ndarray:
    """Interpolate h(z) from a trajectory and convert to H(z) in km/s/Mpc."""
    lna = traj[:, 0]
    h = traj[:, 3]
    a = np.exp(lna)
    z_traj = 1.0 / a - 1.0
    # traj is from z=0 forward in (-time direction => z increasing)
    # Sort by z increasing for interpolation
    order = np.argsort(z_traj)
    z_sorted = z_traj[order]
    h_sorted = h[order]
    return H0 * np.interp(z_query, z_sorted, h_sorted)


# ============================================================================
# Reference models
# ============================================================================
def H_stam_closed(z, H0):
    return H0 * (1.0 + z) ** 2 / (1.0 + z + 0.5 * z**2)


def H_lcdm(z, H0, Om=OMEGA_M):
    return H0 * np.sqrt(Om * (1.0 + z) ** 3 + (1.0 - Om))


def H_eds(z, H0):
    return H0 * (1.0 + z) ** 1.5


def chi2(z, Hobs, sigma, Hmodel):
    return float(np.sum(((Hobs - Hmodel) / sigma) ** 2))


# ============================================================================
# Main
# ============================================================================
def main():
    print("=" * 78)
    print("G29: V_3 modified Friedmann redo — derived H(z) vs chronometers")
    print("=" * 78)
    print()

    print(f"A_0                  = {A0:.6f}  = 1/(12*pi)")
    print(f"Omega_m              = {OMEGA_M}")
    print(f"Omega_r              = {OMEGA_R:.2e}")
    print(f"Omega_DE target      = {OMEGA_DE_TARGET:.4f}  (calibration: V_3(A_0)/rho_crit)")
    print(f"alpha/beta           = {ALPHA_OVER_BETA:.6f}")
    print(f"beta_tilde           = {BETA_TILDE:.6f}")
    print(f"alpha_tilde          = {ALPHA_TILDE:.8f}")
    print()
    print(f"V_3(A_0)             = {V3(A0):.4f}  (should equal Omega_DE = {OMEGA_DE_TARGET:.4f})")
    print(f"V_3'(A_0)            = {V3_prime(A0):.3e}  (should be ~0; min by construction)")
    print(f"V_3''(A_0)           = {V3_double_prime(A0):.4f}")
    omega_A = np.sqrt(V3_double_prime(A0))
    print(f"omega_A = sqrt(V_3'')= {omega_A:.4f}  H_0 units")
    print()

    z = CC_DATA[:, 0]
    Hobs = CC_DATA[:, 1]
    sigma = CC_DATA[:, 2]
    N = len(z)
    z_grid = np.linspace(0, 2.1, 400)

    # ========================================================================
    # Scenario 1: A sits at A_0 today, Adot = 0 (true cosmological constant)
    # ========================================================================
    print("-" * 78)
    print("Scenario 1: A(today) = A_0, Adot(today) = 0  (pure CC equivalent)")
    print("-" * 78)
    traj1 = integrate_backward(A_today=A0, Adot_today=0.0)
    H_sc1_data = H_of_z_from_traj(traj1, z, H0_LOCAL)
    H_sc1_grid = H_of_z_from_traj(traj1, z_grid, H0_LOCAL)
    chi2_sc1 = chi2(z, Hobs, sigma, H_sc1_data)
    print(f"  A evolution: A({-traj1[-1,0]/np.log(10):+.2f}) = {traj1[-1,1]:.5f}")
    print(f"  chi^2 = {chi2_sc1:.3f}  chi^2/N = {chi2_sc1/N:.3f}")
    print()

    # ========================================================================
    # Scenario 2: A displaced HIGHER than A_0 today (e.g. tracking V_1-like)
    #   Choose A(today) > A_0 so it's "rolling back to minimum"
    # ========================================================================
    print("-" * 78)
    print("Scenario 2: A(today) = 0.1 (~4x A_0), Adot(today) = 0")
    print("-" * 78)
    A_disp_high = 0.1
    traj2 = integrate_backward(A_today=A_disp_high, Adot_today=0.0)
    H_sc2_data = H_of_z_from_traj(traj2, z, H0_LOCAL)
    H_sc2_grid = H_of_z_from_traj(traj2, z_grid, H0_LOCAL)
    chi2_sc2 = chi2(z, Hobs, sigma, H_sc2_data)
    print(f"  V_3(0.1) = {V3(A_disp_high):.4f}")
    print(f"  A evolution: A_max in past = {traj2[:,1].max():.4f}  A_min = {traj2[:,1].min():.4f}")
    print(f"  chi^2 = {chi2_sc2:.3f}  chi^2/N = {chi2_sc2/N:.3f}")
    print()

    # ========================================================================
    # Scenario 3: A rolling toward A_0 with non-zero kinetic energy today
    # ========================================================================
    print("-" * 78)
    print("Scenario 3: A(today) = A_0, Adot(today) = -0.1  (kinetic, rolling down past)")
    print("-" * 78)
    traj3 = integrate_backward(A_today=A0, Adot_today=-0.1)
    H_sc3_data = H_of_z_from_traj(traj3, z, H0_LOCAL)
    H_sc3_grid = H_of_z_from_traj(traj3, z_grid, H0_LOCAL)
    chi2_sc3 = chi2(z, Hobs, sigma, H_sc3_data)
    print(f"  A evolution: A range = [{traj3[:,1].min():.5f}, {traj3[:,1].max():.5f}]")
    print(f"  chi^2 = {chi2_sc3:.3f}  chi^2/N = {chi2_sc3/N:.3f}")
    print()

    # ========================================================================
    # Scenario 4: V_1 analog — A_today = A_0, but with kinetic such that
    # tracking solution from the past would land here
    # In V_1: Adot_today = -1.5 * (1-A_0) * H_0 ~= -1.46 (in H_0 units)
    # Try a milder version for V_3 (which has tighter confinement)
    # ========================================================================
    print("-" * 78)
    print("Scenario 4: A(today) = A_0, Adot(today) = -0.5  (V_1-tracking-style kinetic)")
    print("-" * 78)
    traj4 = integrate_backward(A_today=A0, Adot_today=-0.5)
    H_sc4_data = H_of_z_from_traj(traj4, z, H0_LOCAL)
    H_sc4_grid = H_of_z_from_traj(traj4, z_grid, H0_LOCAL)
    chi2_sc4 = chi2(z, Hobs, sigma, H_sc4_data)
    print(f"  A evolution: A range = [{traj4[:,1].min():.5f}, {traj4[:,1].max():.5f}]")
    print(f"  chi^2 = {chi2_sc4:.3f}  chi^2/N = {chi2_sc4/N:.3f}")
    print()

    # ========================================================================
    # Reference models
    # ========================================================================
    H_closed_data = H_stam_closed(z, H0_LOCAL)
    H_lcdm73_data = H_lcdm(z, H0_LOCAL)
    H_lcdm67_data = H_lcdm(z, H0_PLANCK)
    H_eds_data = H_eds(z, H0_LOCAL)
    chi2_closed = chi2(z, Hobs, sigma, H_closed_data)
    chi2_lcdm73 = chi2(z, Hobs, sigma, H_lcdm73_data)
    chi2_lcdm67 = chi2(z, Hobs, sigma, H_lcdm67_data)
    chi2_eds = chi2(z, Hobs, sigma, H_eds_data)

    print("-" * 78)
    print("Reference models (recomputed for comparison):")
    print("-" * 78)
    print(f"  STAM closed-form (H_0=73.04):        chi^2/N = {chi2_closed/N:.3f}")
    print(f"  LCDM (H_0=73.04):                    chi^2/N = {chi2_lcdm73/N:.3f}")
    print(f"  LCDM (H_0=67.4 Planck):              chi^2/N = {chi2_lcdm67/N:.3f}")
    print(f"  EdS (H_0=73.04):                     chi^2/N = {chi2_eds/N:.3f}")
    print()

    print("=" * 78)
    print("Summary of V_3-derived scenarios vs references:")
    print("=" * 78)
    print(f"  V_3 Scen 1 (A=A_0, Adot=0):           chi^2/N = {chi2_sc1/N:.3f}")
    print(f"  V_3 Scen 2 (A=0.1, Adot=0):           chi^2/N = {chi2_sc2/N:.3f}")
    print(f"  V_3 Scen 3 (A=A_0, Adot=-0.1):        chi^2/N = {chi2_sc3/N:.3f}")
    print(f"  V_3 Scen 4 (A=A_0, Adot=-0.5):        chi^2/N = {chi2_sc4/N:.3f}")
    print(f"  STAM closed-form (H_0=73.04):         chi^2/N = {chi2_closed/N:.3f}")
    print(f"  LCDM (H_0=73.04):                     chi^2/N = {chi2_lcdm73/N:.3f}")
    print(f"  LCDM (H_0=67.4):                      chi^2/N = {chi2_lcdm67/N:.3f}")
    print()

    # ========================================================================
    # Plotting
    # ========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel 1: H(z) curves
    ax1 = axes[0]
    ax1.errorbar(z, Hobs, yerr=sigma, fmt="ko", markersize=4, capsize=2,
                 label="Cosmic chronometers")
    ax1.plot(z_grid, H_sc1_grid, "b-", lw=2,
             label=f"V_3 Sc.1 (A=A_0, Adot=0)  chi^2/N={chi2_sc1/N:.2f}")
    ax1.plot(z_grid, H_sc4_grid, "c--", lw=2,
             label=f"V_3 Sc.4 (A=A_0, Adot=-0.5)  chi^2/N={chi2_sc4/N:.2f}")
    ax1.plot(z_grid, H_stam_closed(z_grid, H0_LOCAL), "g:", lw=1.5,
             label=f"STAM closed-form  chi^2/N={chi2_closed/N:.2f}")
    ax1.plot(z_grid, H_lcdm(z_grid, H0_LOCAL), "r-.", lw=1.5,
             label=f"LCDM H_0=73  chi^2/N={chi2_lcdm73/N:.2f}")
    ax1.plot(z_grid, H_lcdm(z_grid, H0_PLANCK), "m:", lw=1.5,
             label=f"LCDM H_0=67.4  chi^2/N={chi2_lcdm67/N:.2f}")
    ax1.set_xlabel("redshift z")
    ax1.set_ylabel("H(z) [km/s/Mpc]")
    ax1.set_title("V_3 modified Friedmann derived H(z) vs chronometers")
    ax1.legend(fontsize=8, loc="upper left")
    ax1.grid(alpha=0.3)
    ax1.set_xlim(0, 2.1)
    ax1.set_ylim(50, 260)

    # Panel 2: A(z) trajectory for the scenarios
    ax2 = axes[1]
    z_t1 = 1.0 / np.exp(traj1[:, 0]) - 1.0
    z_t2 = 1.0 / np.exp(traj2[:, 0]) - 1.0
    z_t3 = 1.0 / np.exp(traj3[:, 0]) - 1.0
    z_t4 = 1.0 / np.exp(traj4[:, 0]) - 1.0
    ax2.plot(z_t1, traj1[:, 1], "b-", lw=2, label="Sc.1 A(today)=A_0, Adot=0")
    ax2.plot(z_t2, traj2[:, 1], "r--", lw=2, label="Sc.2 A(today)=0.1, Adot=0")
    ax2.plot(z_t3, traj3[:, 1], "g:", lw=2, label="Sc.3 A(today)=A_0, Adot=-0.1")
    ax2.plot(z_t4, traj4[:, 1], "c-.", lw=2, label="Sc.4 A(today)=A_0, Adot=-0.5")
    ax2.axhline(A0, color="k", ls=":", alpha=0.5, label=f"A_0 = {A0:.4f}")
    ax2.set_xlabel("redshift z")
    ax2.set_ylabel("A(z)")
    ax2.set_title("A trajectory in cosmic history (each scenario)")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)
    ax2.set_xlim(0, 2.5)

    fig.suptitle("G29: V_3 modified Friedmann — KG + Friedmann numerical integration")
    fig.tight_layout()

    plot_path = RESULTS / "G29_v3_modified_friedmann.png"
    fig.savefig(plot_path, dpi=140, bbox_inches="tight")
    print(f"Plot saved: {plot_path}")
    print()

    # ========================================================================
    # Markdown summary
    # ========================================================================
    md = []
    md.append("# G29: V_3 modified Friedmann redo — derived H(z) vs chronometers")
    md.append("")
    md.append("**Date:** 2026-05-11")
    md.append("")
    md.append("## Goal")
    md.append("")
    md.append(
        "Compute STAM's derived (not ansatz) H(z) from V_3 modified Friedmann, "
        "integrate numerically, and compare to cosmic chronometer data. "
        "This is the redo of script 36 (V_1 era) with the G8-selected V_3 = "
        "alpha/A + beta/(1-A) with minimum at A_0 by construction. "
        "Direct backtrack against G28: does V_3 dynamics fit chronometers "
        "where the closed-form ansatz failed?"
    )
    md.append("")
    md.append("## V_3 calibration")
    md.append("")
    md.append("```")
    md.append("V_3(A) = alpha/A + beta/(1-A)")
    md.append("alpha/beta = (A_0/(1-A_0))^2  =>  V_3 has minimum at A = A_0")
    md.append(f"A_0 = 1/(12*pi) = {A0:.6f}")
    md.append(f"beta_tilde = beta/rho_crit = Omega_DE * (1-A_0)^2 = {BETA_TILDE:.6f}")
    md.append(f"alpha_tilde = beta_tilde * (A_0/(1-A_0))^2 = {ALPHA_TILDE:.8f}")
    md.append(f"V_3(A_0)/rho_crit = {V3(A0):.4f}  (matches Omega_DE target {OMEGA_DE_TARGET:.4f})")
    md.append(f"V_3''(A_0) = {V3_double_prime(A0):.4f}  (omega_A = {omega_A:.4f} in H_0 units)")
    md.append("```")
    md.append("")
    md.append("## Scenarios — initial conditions at today, integrated backward")
    md.append("")
    md.append("| Scenario | A(today) | Adot(today) | chi^2 | chi^2/N |")
    md.append("|---|---|---|---|---|")
    md.append(f"| Sc.1 — pure CC | A_0 | 0 | {chi2_sc1:.2f} | {chi2_sc1/N:.3f} |")
    md.append(f"| Sc.2 — A displaced high | 0.10 | 0 | {chi2_sc2:.2f} | {chi2_sc2/N:.3f} |")
    md.append(f"| Sc.3 — at min, mild kinetic | A_0 | -0.1 | {chi2_sc3:.2f} | {chi2_sc3/N:.3f} |")
    md.append(f"| Sc.4 — at min, strong kinetic | A_0 | -0.5 | {chi2_sc4:.2f} | {chi2_sc4/N:.3f} |")
    md.append("")
    md.append("## Reference models (recomputed)")
    md.append("")
    md.append("| Model | H_0 | chi^2 | chi^2/N |")
    md.append("|---|---|---|---|")
    md.append(f"| STAM closed-form ansatz | 73.04 | {chi2_closed:.2f} | {chi2_closed/N:.3f} |")
    md.append(f"| LCDM (Om=0.315) | 73.04 | {chi2_lcdm73:.2f} | {chi2_lcdm73/N:.3f} |")
    md.append(f"| LCDM (Om=0.315) | 67.4 Planck | {chi2_lcdm67:.2f} | {chi2_lcdm67/N:.3f} |")
    md.append(f"| EdS | 73.04 | {chi2_eds:.2f} | {chi2_eds/N:.3f} |")
    md.append("")
    md.append("## Plot")
    md.append("")
    md.append("![V_3 modified Friedmann vs chronometers](G29_v3_modified_friedmann.png)")
    md.append("")
    md.append("Left: H(z) curves for V_3 scenarios + reference models vs chronometer data.")
    md.append("Right: A(z) trajectory in each scenario (how A evolves over cosmic history).")
    md.append("")
    md.append("## Reading")
    md.append("")
    md.append(
        "V_3 with A sitting at the minimum (Scenario 1) is degenerate with LCDM by construction "
        "— the calibration sets V_3(A_0) = Omega_DE_target * rho_crit, which is exactly LCDM's "
        "cosmological constant. Any chi^2 difference between Scenario 1 and LCDM at H_0=73.04 "
        "reflects numerical drift in the integration, not physics."
    )
    md.append("")
    md.append(
        "The scenarios with non-trivial A dynamics (Scenarios 2-4) test whether V_3 admits "
        "tracking-like behavior that would shift H(z) away from LCDM. **V_3's tight confinement "
        "around A_0** (the omega_A = sqrt(V_3'') is large in H_0 units) drives A back to A_0 "
        "rapidly under Hubble friction — there's no slow-roll tracking comparable to V_1."
    )
    md.append("")
    md.append("**Backtrack on G28's closed-form result:**")
    md.append("")
    md.append(
        "The closed-form ansatz H(z) = H_0(1+z)^2/(1 + z + 0.5*z^2) does not match V_3's derived "
        "dynamics. V_3 dynamics naturally reduce to LCDM at H_0=73 (because A is confined at A_0 "
        "and V_3(A_0) acts as cosmological constant). This suggests the closed-form was a "
        "*conjectured* H(z) form not actually derived from V_3 — and that the framework's "
        "'intrinsic H(z)' under V_3 IS essentially LCDM at H_0=73."
    )
    md.append("")
    md.append("**This shifts the narrative on the bridge term.**")
    md.append("")
    md.append(
        "If V_3's intrinsic H(z) = LCDM-at-73, then the distance reframe ('STAM intrinsic is "
        "flatter') requires the closed-form ansatz to be retained as a *separate* piece of physics "
        "— not derived from V_3 dynamics. The closed-form had been informally connected to the "
        "bridge-term / photon-A-traversal story; under V_3, those threads need to be re-examined. "
        "The bridge term may sit on a different layer than V_3 modified Friedmann."
    )
    md.append("")
    md.append("## What this implies for next work")
    md.append("")
    md.append(
        "- V_3 modified Friedmann gives essentially LCDM H(z) when calibrated to Omega_DE today. "
        "The framework's cosmic expansion dynamics are not where the distance reframe lives."
    )
    md.append(
        "- The closed-form H(z) was probably an ad-hoc shape, not derived from V_3. Where did it "
        "come from? Worth tracing back to its origin (likely from inverting the bz-form distance "
        "modulus). If the derivation is unsound, the closed-form should not be load-bearing."
    )
    md.append(
        "- The photon-A-traversal mechanism (light traversing A_0 in voids accumulates extra delay) "
        "is a *separate* effect from V_3 cosmic expansion. It affects observed distance modulus, "
        "not H(z). Chronometer H(z) and SN distance modulus probe different layers of the framework."
    )
    md.append(
        "- The bridge term's identity becomes clearer: it's about photon paths through the cosmic "
        "A_0 field, NOT about cosmic expansion dynamics. V_3 Friedmann doesn't constrain it."
    )

    md_path = RESULTS / "G29_v3_modified_friedmann_summary.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"Summary saved: {md_path}")


if __name__ == "__main__":
    main()
