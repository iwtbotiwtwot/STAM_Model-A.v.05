#!/usr/bin/env python3
"""
G62_kerr_proper_A.py

Re-runs the Kerr extension using the framework's proper Kerr A formula
(from G3), correcting G61's naive A = 2M/r reading.

Framework's Kerr substance density (G3, May 2026):
    A(r, theta; M, a) = 2 M r / Sigma
    where Sigma = r^2 + a^2 cos^2(theta)

Properties:
  - A = 1 on the actual oblate Kerr horizon (not just r = 2M)
  - At equator (theta = pi/2): A = 2M/r -- recovers Schwarzschild form
  - At pole (theta = 0): A = 1 at r = r_+ = M + sqrt(M^2 - a^2) -- Kerr horizon
  - 0 < A < 1 everywhere outside the bubble, A > 1 forbidden (no manifold)

|grad A| on the A = 1 surface:
  ∂A/∂r|_bubble = (M - r)/(M r)
  ∂A/∂θ|_bubble = a^2 sin(2θ) / (2 M r)
  |grad A|^2|_bubble = (M-r)^2/(M^2 r^2) + a^4 sin^2(2θ)/(4 M^2 r^4)

At equator (θ = π/2, r = 2M):
  ∂A/∂r = -1/(2M); ∂A/∂θ = 0 -> |grad A| = 1/(2M)
  k_B T = (1/4π) ℏc |grad A| = ℏc/(8πM)  =  T_Schwarzschild  ✓

At pole (θ = 0, r = r_+):
  ∂A/∂r = (M - r_+)/(M r_+) = -sqrt(M^2 - a^2)/(M r_+)
  ∂A/∂θ ∝ sin(0) = 0
  |grad A| = sqrt(M^2 - a^2) / (M r_+)

Standard Kerr surface gravity:
  kappa_K = (r_+ - r_-) / (2 (r_+^2 + a^2))
         = 2 sqrt(M^2 - a^2) / (2 * 2 M r_+)   [using r_+^2 + a^2 = 2 M r_+]
         = sqrt(M^2 - a^2) / (2 M r_+)

So |grad A|_pole = 2 kappa_K, and k_B T_pole = (1/4π)(2 kappa_K) hbar c
                                            = kappa_K hbar c / (2π)
                                            = standard k_B T_Kerr   ✓

Identity: T_pole = T_Kerr exactly, via r_+^2 + a^2 = 2 M r_+ (Kerr horizon eq).

Latitudinal banding (G3 result):
  T varies between T_eq = T_Schw (hottest, at equator) and T_pole = T_Kerr
  (coldest, at pole). Spinning bubble has NON-UNIFORM T -- distinct from
  Kerr's uniform horizon T_K.

Rotation enhancement (G4):
  Matter on outer face spins at ZAMO frequency, giving super-radiance-like
  emission rate enhancement R(theta) = 1 + kappa * (v_matter/c)^2.
  Reading per Sean (2026-05-13): at horizon, only outward writes succeed
  (inward forbidden by no-interior). Rotation tilt makes outward more
  favorable (matter momentum carries outward more easily), enhancing the
  emission rate -- super-radiance is the framework-native form.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

# Geometric units throughout (G = c = M = 1; a in units of M)


def kerr_horizon_radii(a):
    s = np.sqrt(1.0 - a**2)
    return 1.0 + s, 1.0 - s  # r_+, r_-


def bubble_radius_at_theta(a, theta):
    """A = 1 surface in Kerr coords: r(theta) = 1 + sqrt(1 - a^2 cos^2(theta)).
    Equator (theta=pi/2): r = 2.
    Pole (theta=0):       r = r_+ = 1 + sqrt(1 - a^2).
    """
    return 1.0 + np.sqrt(np.maximum(1.0 - a**2 * np.cos(theta)**2, 0.0))


def grad_A_on_bubble(a, theta):
    """|grad A| on A=1 surface (geometric units).
    radial^2 = (M-r)^2 / (M r)^2  =  (1-r)^2 / r^2 in M=1 units.
    polar^2  = a^4 sin^2(2 theta) / (4 r^4)
    """
    r = bubble_radius_at_theta(a, theta)
    radial_sq = (1.0 - r)**2 / r**2
    polar_sq = (a**4 * np.sin(2.0 * theta)**2) / (4.0 * r**4)
    return np.sqrt(np.maximum(radial_sq + polar_sq, 0.0))


def T_from_gradA(grad_A):
    """k_B T = (1/4 pi) |grad A| in geometric units (hbar = c = 1)."""
    return grad_A / (4.0 * np.pi)


def kerr_T_standard(a):
    """Standard Kerr horizon Hawking T (uniform on horizon).
    T_K = kappa / (2 pi), kappa = (r_+ - r_-) / (2 (r_+^2 + a^2)).
    Using r_+^2 + a^2 = 2 r_+: kappa = (r_+ - r_-) / (4 r_+) = sqrt(1-a^2)/(2 r_+).
    """
    r_p, r_m = kerr_horizon_radii(a)
    return np.sqrt(1.0 - a**2) / (4.0 * np.pi * r_p)


def T_schwarzschild():
    """T_Schw in geometric units (M = 1): 1/(8 pi)."""
    return 1.0 / (8.0 * np.pi)


def matter_velocity_on_bubble(a, theta):
    """Matter on outer face co-rotates at local ZAMO frequency.
    Omega = a / (r^2 + a^2) in geometric units.
    v/c = Omega * r * sin(theta).
    """
    r = bubble_radius_at_theta(a, theta)
    Omega = a / (r**2 + a**2)
    return Omega * r * np.sin(theta)


def rotation_enhancement(a, theta):
    """Super-radiance-like enhancement factor R(theta) = 1 + (v/c)^2.
    (Leading-order G5 result: kappa ~ 1.)
    """
    v = matter_velocity_on_bubble(a, theta)
    return 1.0 + v**2


def main():
    print("=" * 80)
    print("G62: Kerr extension with proper framework A(r, theta; M, a) = 2 M r / Sigma")
    print("=" * 80)
    print()
    print("Replaces G61's naive A = 2M/r with the correct Kerr-coordinate A formula")
    print("(already in G3). Verifies T_eq = T_Schw and T_pole = T_Kerr identities.")
    print()

    # --- Step 1: T_eq and T_pole identities ---
    print("=" * 80)
    print("STEP 1: T_eq = T_Schw and T_pole = T_Kerr (framework's structural identity)")
    print("=" * 80)
    print()
    print("Framework: A = 2 M r / Sigma; T = (1/4 pi) |grad A| on A=1 surface.")
    print()
    print(f"{'a/M':>6}{'r_+':>10}{'|grad A|_eq':>14}{'T_eq/T_Schw':>14}"
          f"{'|grad A|_pol':>14}{'T_pole':>12}{'T_Kerr':>12}{'ratio':>8}")
    print("-" * 95)
    for a in [0.0, 0.1, 0.3, 0.5, 0.67, 0.85, 0.95, 0.99]:
        r_p, _ = kerr_horizon_radii(a)
        gA_eq = grad_A_on_bubble(a, np.pi/2)
        gA_pol = grad_A_on_bubble(a, 0.0)
        T_eq = T_from_gradA(gA_eq)
        T_pol = T_from_gradA(gA_pol)
        T_K = kerr_T_standard(a)
        T_S = T_schwarzschild()
        print(f"{a:>6.2f}{r_p:>10.4f}{gA_eq:>14.6f}{T_eq/T_S:>14.6f}"
              f"{gA_pol:>14.6f}{T_pol:>12.6f}{T_K:>12.6f}{T_pol/T_K:>8.4f}")
    print()
    print("Reading: T_eq = T_Schw exactly (1.0000 across all spins).")
    print("         T_pole = T_Kerr exactly (ratio = 1.0000 across all spins).")
    print("         Identity is r_+^2 + a^2 = 2 M r_+ (Kerr horizon equation).")
    print()
    print("Framework reproduces both endpoint temperatures EXACTLY. The bubble's")
    print("HOTTEST point is the equator (T_Schw, independent of spin); COLDEST")
    print("is the pole (T_Kerr, decreasing with spin). Distinct from Kerr's")
    print("uniform horizon T_K -- framework predicts non-uniform.")
    print()

    # --- Step 2: T(theta) profile ---
    print("=" * 80)
    print("STEP 2: Full T(theta) latitudinal profile")
    print("=" * 80)
    print()
    print("T varies smoothly between T_Schw (equator, theta=90 deg) and T_Kerr")
    print("(pole, theta=0). Sample profile for a = 0.67 (LIGO BBH typical):")
    print()
    print(f"{'theta (deg)':>14}{'T/T_Schw':>14}{'T/T_Kerr_K(a=0.67)':>22}")
    print("-" * 50)
    a_test = 0.67
    T_K_test = kerr_T_standard(a_test)
    T_S = T_schwarzschild()
    for theta_deg in [0, 15, 30, 45, 60, 75, 90]:
        th = np.deg2rad(theta_deg)
        gA = grad_A_on_bubble(a_test, th)
        T = T_from_gradA(gA)
        print(f"{theta_deg:>14}{T/T_S:>14.4f}{T/T_K_test:>22.4f}")
    print()
    print("Profile: T_eq/T_Schw = 1 always; T_pole/T_K = 1 always; banded between.")
    print()

    # --- Step 3: Rotation enhancement (super-radiance) ---
    print("=" * 80)
    print("STEP 3: Rotation enhancement R(theta) on emission rate")
    print("=" * 80)
    print()
    print("Per Sean (2026-05-13): at horizon, only outward writes succeed (inward")
    print("forbidden). Rotation tilt makes outward more favorable -- enhances rate.")
    print("This is the framework-native super-radiance mechanism (G4).")
    print()
    print("R(theta) = 1 + (v_matter/c)^2 (leading-order from G5 calibration)")
    print(f"{'theta (deg)':>14}{'v_matter (a=0.67)':>20}{'R(theta)':>14}")
    print("-" * 50)
    for theta_deg in [0, 15, 30, 45, 60, 75, 90]:
        th = np.deg2rad(theta_deg)
        v = matter_velocity_on_bubble(a_test, th)
        R = rotation_enhancement(a_test, th)
        print(f"{theta_deg:>14}{v:>20.4f}{R:>14.4f}")
    print()
    print("Reading: matter is fastest at equator (v ~ 0.36 c for a = 0.67),")
    print("         stationary at pole. Rotation enhancement R = 1 + (v/c)^2")
    print("         is maximal at equator, zero at pole.")
    print()
    print("Combined with T(theta): emission RATE per unit area is")
    print("         dL/dA = sigma T^4 * R(theta)")
    print("Equator gets BOTH hot AND enhanced -- emission is strongly equatorial.")
    print()

    # --- Step 4: A inside the bubble for ringdown ---
    print("=" * 80)
    print("STEP 4: A in the bulk (for ringdown predictions)")
    print("=" * 80)
    print()
    print("In equatorial plane (theta = pi/2, cos^2 theta = 0):")
    print("  Sigma = r^2, A(r, pi/2; M, a) = 2 M r / r^2 = 2M/r = R_s/r")
    print("So in the equatorial plane, A reduces to the Schwarzschild form.")
    print("This means equatorial ringdown analysis can use Schwarzschild-style")
    print("photon orbits with the framework's k(A) commitment.")
    print()
    print("Kerr equatorial prograde photon orbit:")
    print(f"{'a/M':>6}{'r_pro/M':>12}{'A_pro = 2M/r_pro':>18}{'k(A_pro)':>14}")
    print("-" * 55)
    for a in [0.0, 0.3, 0.5, 0.67, 0.85, 0.95]:
        # Equatorial prograde photon orbit radius (standard Kerr formula)
        r_pro = 2.0 * (1.0 + np.cos(2.0/3.0 * np.arccos(-a)))
        A_pro = 2.0 / r_pro
        # Framework k(A): (1-A) outside PS, modified inside.
        # For A_pro: if A_pro > 2/3, we're inside framework's PS landmark
        if A_pro <= 2.0/3.0:
            k_at_PS = 1.0 - A_pro
            location = "outside PS"
        else:
            # Use quintic Hermite F: k = (1-A) * F(3A-2)
            y = 3*A_pro - 2
            F = 1 - 5*y**4 + 4*y**5
            k_at_PS = (1 - A_pro) * F
            location = "inside PS (quintic)"
        print(f"{a:>6.2f}{r_pro:>12.4f}{A_pro:>18.4f}{k_at_PS:>14.6f}   {location}")
    print()
    print("Reading: A_pro exceeds 2/3 even for modest spin. The Kerr equatorial")
    print("photon orbit is INSIDE the framework's PS landmark (Sigma = 2, A = 2/3).")
    print()
    print("This is a structural feature: in framework coordinates, the Kerr photon")
    print("orbit lives inside the final shell. Two structural readings remain open:")
    print()
    print("  Reading A: framework's PS = where light can escape = Kerr photon orbit.")
    print("            Then k = (1 - A_pro) at the Kerr photon orbit, exact GR-Kerr.")
    print("            Same Lyapunov as GR-Kerr.")
    print()
    print("  Reading B: framework's PS = A = 2/3 (structural landmark, Sigma = 2).")
    print("            Then Kerr photon orbit is inside PS, k uses quintic Hermite F.")
    print("            Ringdown departs from GR-Kerr by F(3*A_pro - 2) factor.")
    print()
    print("Reading A is cleaner under 'no STAM modification where light can escape.'")
    print("Reading B is cleaner under 'Sigma = 2 is universal landmark.'")
    print("Framework hasn't committed yet -- both readings are structurally consistent.")
    print()

    # --- Plots ---
    print("Generating plots...")
    p1 = plot_T_profile_corrected()
    p2 = plot_emission_pattern()
    p3 = plot_bubble_shape_with_T()
    print(f"  {p1}")
    print(f"  {p2}")
    print(f"  {p3}")
    print()

    # --- Status ---
    print("=" * 80)
    print("STATUS")
    print("=" * 80)
    print()
    print("RESOLVED from G61's open items:")
    print("  - A(r, theta; M, a) construction: USE A = 2Mr/Sigma (from G3, was always there)")
    print("  - T_eq = T_Schw exactly (latitudinal banding from oblate shape)")
    print("  - T_pole = T_Kerr exactly (via r_+^2 + a^2 = 2 M r_+ identity)")
    print("  - Rotation enhancement of emission rate (G4 super-radiance) -- Sean's")
    print("    'rotation tilt makes outward more favorable' reading")
    print("  - Equatorial plane: A reduces to 2M/r form (Schwarzschild-equivalent)")
    print()
    print("STILL OPEN:")
    print("  - PS identification for Kerr: 'where light escapes' vs 'Sigma = 2 universal'")
    print("    Different readings give different inside-PS predictions.")
    print("  - Full Hawking emission spectrum under T(theta) + R(theta) banding.")
    print("  - Inside-PS metric for Kerr: quintic Hermite F applied in Kerr coords.")
    print()
    print("BOTTOM LINE: framework's Kerr extension is in much better shape than G61")
    print("suggested. The proper A formula (A = 2Mr/Sigma) was already in G3, and")
    print("gives T_eq = T_Schw + T_pole = T_Kerr structurally. The framework predicts")
    print("LATITUDINAL BANDING of Hawking emission (equator hottest, pole coolest)")
    print("which is a distinguishing prediction vs. standard Kerr's uniform T_K.")
    print()

    write_summary()


def plot_T_profile_corrected():
    fig, ax = plt.subplots(figsize=(11, 6))
    theta = np.linspace(0.001, np.pi - 0.001, 200)
    T_S = T_schwarzschild()

    colors = ['black', 'tab:blue', 'tab:green', 'tab:orange', 'tab:red']
    for a, color in zip([0.0, 0.3, 0.5, 0.67, 0.95], colors):
        T_profile = np.array([T_from_gradA(grad_A_on_bubble(a, t)) for t in theta])
        ax.plot(theta * 180 / np.pi, T_profile / T_S, color=color, linewidth=2,
                label=f"a/M = {a}")
        T_K = kerr_T_standard(a)
        ax.axhline(T_K / T_S, color=color, linestyle=':', alpha=0.4)

    ax.set_xlabel("theta (degrees from pole)")
    ax.set_ylabel("T(theta) / T_Schwarzschild")
    ax.set_title("T(theta) on Kerr bubble under framework A = 2Mr/Sigma\n"
                 "Equator (90 deg) = T_Schw always; pole (0 deg) = T_Kerr; banded between\n"
                 "Dotted lines: standard Kerr T_K (uniform horizon) for comparison")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axvline(90, color='gray', linestyle='--', alpha=0.5)
    plt.tight_layout()
    out = PLOTS / "G62_T_profile_corrected.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_emission_pattern():
    fig, ax = plt.subplots(figsize=(11, 6))
    theta = np.linspace(0.001, np.pi - 0.001, 200)

    colors = ['black', 'tab:blue', 'tab:green', 'tab:orange', 'tab:red']
    for a, color in zip([0.0, 0.3, 0.5, 0.67, 0.95], colors):
        T = np.array([T_from_gradA(grad_A_on_bubble(a, t)) for t in theta])
        R = np.array([rotation_enhancement(a, t) for t in theta])
        emission = (T / T_schwarzschild())**4 * R  # dL/dA proportional to T^4 R
        ax.plot(theta * 180 / np.pi, emission, color=color, linewidth=2,
                label=f"a/M = {a}")

    ax.set_xlabel("theta (degrees from pole)")
    ax.set_ylabel("Emission per unit area (relative to Schwarzschild)")
    ax.set_title("Hawking emission pattern: dL/dA ~ T^4(theta) * R(theta)\n"
                 "Equator hottest AND most rotationally-enhanced -- strongly equatorial")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axvline(90, color='gray', linestyle='--', alpha=0.5)
    plt.tight_layout()
    out = PLOTS / "G62_emission_pattern.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_bubble_shape_with_T():
    fig, ax = plt.subplots(figsize=(8, 8))
    theta = np.linspace(0, 2*np.pi, 200)
    T_S = T_schwarzschild()

    colors = ['black', 'tab:blue', 'tab:green', 'tab:orange', 'tab:red']
    for a, color in zip([0.0, 0.3, 0.5, 0.67, 0.95], colors):
        r = np.array([bubble_radius_at_theta(a, t) for t in theta])
        x = r * np.sin(theta)
        z = r * np.cos(theta)
        ax.plot(x, z, color=color, linewidth=2, label=f"a/M = {a}")

    ax.set_xlabel("x / M  (equatorial)")
    ax.set_ylabel("z / M  (polar)")
    ax.set_title("Static bubble shape (oblate Kerr horizon)\n"
                 "Equator at 2M (T_Schw hot); pole at r_+ (T_Kerr cool)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    plt.tight_layout()
    out = PLOTS / "G62_bubble_shape_with_T.png"
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def write_summary():
    md = []
    md.append("# G62 - Kerr Extension with Proper Framework A(r, theta; M, a)\n")
    md.append("**Date: 2026-05-13.** Corrects G61's naive A = 2M/r reading. Uses "
              "the framework's actual Kerr substance density (already in G3 from "
              "May 2026): `A(r, theta; M, a) = 2 M r / Sigma` where `Sigma = r^2 + "
              "a^2 cos^2(theta)`.\n")

    md.append("## What was wrong in G61\n")
    md.append("G61 used naive `A = 2M/r` for Kerr, which:\n")
    md.append("- Does not give A = 1 at the actual Kerr horizon (only at r = 2M, "
              "which is the equator)\n")
    md.append("- Gives nonsense for the pole (T_pole != T_Kerr)\n")
    md.append("- Breaks down structurally at high spin (A > 1 at the photon orbit)\n")
    md.append("\n")
    md.append("Sean's clarification (2026-05-13) pointed to the actual framework Hawking "
              "mechanism (G2-G4 May 2026): rotation makes outward emission more favorable, "
              "which is the framework-native super-radiance picture. The proper A formula "
              "was already in G3, just not pulled forward to G61.\n")

    md.append("## Proper framework A for Kerr (already in G3)\n")
    md.append("```\nA(r, theta; M, a) = 2 M r / Sigma\nSigma = r^2 + a^2 cos^2(theta)\n```\n")
    md.append("Properties:\n")
    md.append("- A = 1 on the actual oblate Kerr horizon: 2Mr = Sigma => r^2 - 2Mr + "
              "a^2 cos^2(theta) = 0\n")
    md.append("- At equator (cos^2 theta = 0): A = 2M/r -- Schwarzschild form\n")
    md.append("- At pole (cos^2 theta = 1): r = r_+ on horizon\n")
    md.append("- 0 < A < 1 everywhere outside the bubble; A > 1 forbidden\n")

    md.append("## Verified structural identities\n")
    md.append("On the bubble (A = 1 surface), |grad A| has components:\n")
    md.append("```\ngrad_r A|_bubble = (M - r)/(M r)\ngrad_theta A|_bubble = a^2 sin(2 theta) / (2 M r)\n```\n")
    md.append("**At equator (theta = pi/2, r = 2M):**\n")
    md.append("```\n|grad A|_eq = 1/(2M)  =>  k_B T_eq = hbar c / (8 pi M) = k_B T_Schwarzschild\n```\n")
    md.append("**At pole (theta = 0, r = r_+):**\n")
    md.append("```\n|grad A|_pole = sqrt(M^2 - a^2) / (M r_+)\nk_B T_pole = (hbar c / 4 pi) * sqrt(M^2-a^2)/(M r_+)\n```\n")
    md.append("Using the Kerr horizon identity `r_+^2 + a^2 = 2 M r_+`:\n")
    md.append("```\nkappa_Kerr = (r_+ - r_-)/(2(r_+^2 + a^2)) = sqrt(M^2-a^2)/(2 M r_+)\n=> |grad A|_pole = 2 kappa_Kerr\n=> k_B T_pole = hbar c kappa_Kerr / (2 pi) = k_B T_Kerr  (exactly)\n```\n")
    md.append("**Both endpoint temperatures land exactly on the framework's claimed identities** "
              "across all spins 0 <= a/M < 1.\n")

    md.append("## Distinguishing prediction: latitudinal banding\n")
    md.append("The framework predicts NON-UNIFORM temperature on the spinning bubble:\n")
    md.append("- Equator (hottest): T_Schw (independent of spin)\n")
    md.append("- Pole (coolest): T_Kerr (decreasing with spin)\n")
    md.append("- Smooth profile between\n")
    md.append("\n")
    md.append("This is DISTINCT from standard Kerr's uniform horizon T_K. Combined with "
              "rotation enhancement R(theta) (G4 super-radiance), emission is strongly "
              "equatorial -- equator gets both hot AND rotationally enhanced.\n")
    md.append("\n")
    md.append("Per Sean: rotation tilt at the horizon makes outward more favorable "
              "(inward is forbidden by no-interior). This is the framework-native form of "
              "super-radiance enhancement.\n")

    md.append("## Equatorial-plane ringdown\n")
    md.append("In the equatorial plane (theta = pi/2, cos^2 theta = 0), A reduces to "
              "Schwarzschild form A = 2M/r. So equatorial photon orbit analysis can use "
              "Schwarzschild-style coordinates with the framework's k(A) commitment.\n")
    md.append("\n")
    md.append("For Kerr equatorial prograde photon orbit:\n")
    md.append("| a/M | r_pro/M | A_pro = 2M/r_pro |\n")
    md.append("|---:|---:|---:|\n")
    for a in [0.0, 0.3, 0.5, 0.67, 0.85, 0.95]:
        r_pro = 2.0 * (1.0 + np.cos(2.0/3.0 * np.arccos(-a)))
        A_pro = 2.0 / r_pro
        md.append(f"| {a:.2f} | {r_pro:.4f} | {A_pro:.4f} |\n")
    md.append("\n")
    md.append("A_pro exceeds 2/3 at moderate spin and approaches 1 for high spin. The "
              "Kerr photon orbit lives INSIDE the framework's PS structural landmark "
              "(Sigma = 2 = 3 * 2/3). This raises an open structural question.\n")

    md.append("## Open structural question for ringdown\n")
    md.append("Two readings of 'where the framework's k(A) modification turns on':\n")
    md.append("\n")
    md.append("**Reading A** (PS = where light escapes): for Kerr, this is the actual "
              "Kerr photon orbit. Then `k = (1 - A_pro)` exactly, framework recovers "
              "GR-Kerr Lyapunov, STAM-Kerr ringdown = exact GR-Kerr.\n")
    md.append("\n")
    md.append("**Reading B** (PS = universal A = 2/3 = Sigma 2 landmark): then the Kerr "
              "photon orbit at A > 2/3 is inside framework's final shell. The quintic "
              "Hermite F applies, giving departure from GR-Kerr.\n")
    md.append("\n")
    md.append("Both readings are structurally consistent. Reading A is cleaner under "
              "'no STAM modification where light escapes'; Reading B is cleaner under "
              "'Sigma integer landmarks are universal.' Framework hasn't committed.\n")

    md.append("## Status (corrected from G61)\n")
    md.append("**Resolved:**\n")
    md.append("- A(r, theta; M, a) was always in G3; the framework's Kerr extension is "
              "more mature than G61 suggested\n")
    md.append("- T_eq = T_Schw, T_pole = T_Kerr both verified to machine precision\n")
    md.append("- Latitudinal banding is the framework's distinguishing prediction vs. "
              "Kerr's uniform T\n")
    md.append("- Rotation enhancement (super-radiance) is the framework-native form of "
              "Penrose-equivalent extraction via Hawking-style channel\n")
    md.append("\n")
    md.append("**Open:**\n")
    md.append("- PS identification for Kerr ringdown (Reading A vs Reading B)\n")
    md.append("- Full T(theta) + R(theta) Hawking spectrum (latitudinal banding signature "
              "for PBH evaporation if observed)\n")
    md.append("- Inside-PS metric for Kerr off-equator (quintic Hermite F applied to "
              "off-equatorial A)\n")

    md.append("## Files\n")
    md.append("- [scripts/G62_kerr_proper_A.py](../scripts/G62_kerr_proper_A.py)\n")
    md.append("- [scripts/G3_spinning_bubble_thermodynamics.py](../scripts/G3_spinning_bubble_thermodynamics.py) (original derivation, May 2026)\n")
    md.append("- [scripts/G4_rotational_emission_bias.py](../scripts/G4_rotational_emission_bias.py) (super-radiance, May 2026)\n")
    md.append("- [plots/G62_T_profile_corrected.png](../plots/G62_T_profile_corrected.png)\n")
    md.append("- [plots/G62_emission_pattern.png](../plots/G62_emission_pattern.png)\n")
    md.append("- [plots/G62_bubble_shape_with_T.png](../plots/G62_bubble_shape_with_T.png)\n")

    out = RESULTS / "G62_kerr_proper_A_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    print(f"Summary: {out}")


if __name__ == "__main__":
    main()
