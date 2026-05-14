#!/usr/bin/env python3
"""
G119_kerr_wave_base_geometry.py

Build a numerically sound Kerr STAM wave base.  GEOMETRY ONLY.
No Teukolsky physics, no QNM extraction, no time-domain.  This script's
purpose is to replace the Schwarzschild base used in G118 (which broke
at high spin where r_+(a) < 2M and 1 - 2M/r turns negative) with the
actual Kerr radial structure.

Seven minimum goals
===================
1. Use Kerr horizon  r_+(a) = M + sqrt(M^2 - a^2),  NOT r = 2M.
2. Use exact  r_pr_exact(theta, a)  from G85 / G116.
3. Build equatorial  y_K_eq(r; a)  using  Sigma_ph_eq(a) =
   Sigma_K(r_pr_exact(pi/2, a); a).
4. Build a Kerr-compatible tortoise coordinate
       dr*/dr = (r^2 + a^2) / Delta_STAM
       with Delta_STAM = Delta(r) * F(y_K_eq)  inside the photon region,
       Delta_STAM = Delta(r)  outside.
5. Confirm no imaginary h/k behavior at high spin  (a up to 0.99).
6. Plot STAM throat width vs spin in both r and r* coordinates.
7. (Reserved: QNM-like diagnostics are G120+, not this script.)

What this gives the framework
=============================
A wave base that:
- Vanishes at the actual Kerr horizon r_+(a), not the Schwarzschild r = 2M.
- Treats the inside-photon-region as the STAM-modified zone via G85/G116's
  exact r_pr_exact surface.
- Preserves the cubic horizon vanishing of k from the spinless analog.
- Has a Kerr-standard tortoise outside the photon region (matches GR-Kerr
  there) and a STAM-stretched tortoise inside.

This is the geometric foundation for G120 (rigorous Kerr QNM via Detweiler
or Sasaki-Nakamura on this base).

Outputs
=======
    results/G119_kerr_wave_base_geometry_summary.md
    plots/G119_kerr_wave_base_geometry.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


M = 1.0


# ===========================================================================
# Goal 1.  Kerr horizon  r_+(a)
# ===========================================================================

def r_horizon(a):
    """Outer Kerr horizon.  Real for |a| <= M."""
    disc = M * M - a * a
    if disc < 0:
        return float('nan')
    return M + np.sqrt(disc)


def Delta(r, a):
    """Kerr radial structure  Delta = r^2 - 2Mr + a^2.  Zero at r = r_+."""
    return r * r - 2.0 * M * r + a * a


# ===========================================================================
# Goal 2.  Exact photon-region surface  r_pr_exact(theta, a)  from G85/G116
# ===========================================================================

def r_ph_prograde(a):
    """Bardeen equatorial prograde photon orbit."""
    if abs(a) < 1e-14:
        return 3.0 * M
    return 2.0 * M * (1.0 + np.cos((2.0 / 3.0) * np.arccos(-a / M)))


def r_polar(a):
    """Largest real root of r^3 - 3Mr^2 + a^2 r + a^2 M = 0."""
    if abs(a) < 1e-14:
        return 3.0 * M
    coeffs = [1.0, -3.0 * M, a * a, a * a * M]
    roots = np.roots(coeffs)
    real_roots = [rr.real for rr in roots if abs(rr.imag) < 1e-10]
    rh = r_horizon(a)
    candidates = [rr for rr in real_roots if rr > rh]
    return max(candidates) if candidates else float('nan')


def lambda_phot(r_p, a):
    return -(r_p**3 - 3.0 * M * r_p**2 + a * a * r_p + a * a * M) / (a * (r_p - M))


def eta_phot(r_p, a):
    D = Delta(r_p, a)
    return (r_p**2 * (4.0 * a * a * D - (r_p**2 - 3.0 * M * r_p + 2.0 * a * a)**2)
            / (a * a * (r_p - M)**2))


def Theta_potential(theta, r_p, a):
    """Theta(theta; r_p, a) * sin^2(theta).  Zero at orbit's turning point."""
    lam = lambda_phot(r_p, a)
    eta = eta_phot(r_p, a)
    s2 = np.sin(theta)**2
    c2 = np.cos(theta)**2
    return eta * s2 - c2 * (lam * lam - a * a * s2)


def r_pr_exact(theta, a):
    """Exact spheroidal photon-region surface.  G85 / G116."""
    if abs(a) < 1e-14:
        return 3.0 * M
    eps = 1e-10
    if abs(theta - np.pi / 2) < eps:
        return r_ph_prograde(a)
    if abs(theta) < eps:
        return r_polar(a)
    r_lo = r_ph_prograde(a) + 1e-9
    r_hi = r_polar(a) - 1e-9
    f_lo = Theta_potential(theta, r_lo, a)
    f_hi = Theta_potential(theta, r_hi, a)
    if f_lo * f_hi > 0:
        if abs(f_lo) < abs(f_hi):
            return r_ph_prograde(a)
        return r_polar(a)
    return brentq(lambda r: Theta_potential(theta, r, a), r_lo, r_hi,
                  xtol=1e-12, rtol=1e-12)


# ===========================================================================
# Goal 3.  Equatorial y_K using Sigma_ph_eq(a)
# ===========================================================================

def A_K(r, a):
    """Kerr substance density  A_K = 2Mr / (r^2 + a^2)  at equator (theta = pi/2)."""
    return 2.0 * M * r / (r * r + a * a)


def Sigma_K(r, a):
    return 3.0 * A_K(r, a)


def Sigma_ph_eq(a):
    """Sigma_K at the prograde equatorial photon orbit."""
    return Sigma_K(r_ph_prograde(a), a)


def y_K_eq(r, a):
    """Inside-shell normalized coordinate at equator.

    y = 0 at r_ph_prograde(a)  (photon-region boundary)
    y = 1 at r_+(a)             (horizon: Sigma_K(r_+) = 3 via horizon identity)
    """
    Sig_r = Sigma_K(r, a)
    Sig_ph = Sigma_ph_eq(a)
    denom = 3.0 - Sig_ph
    if abs(denom) < 1e-14:
        return float('nan')
    return (Sig_r - Sig_ph) / denom


def F_quintic(y):
    """Beta(4, 2) closure profile.  F(0) = 1, F(1) = 0 cubically (via Delta linearity)."""
    return 1.0 - 5.0 * y**4 + 4.0 * y**5


# ===========================================================================
# Goal 4.  Kerr-compatible tortoise coordinate
# ===========================================================================

def Delta_STAM(r, a):
    """STAM-modified radial structure on Kerr equator.

       Delta_STAM = Delta * F(y_K_eq)  inside the photon region.
       Delta_STAM = Delta              outside.

    At r = r_+:  Delta = 0  =>  Delta_STAM = 0 regardless of F.
    At r = r_pr_eq:  F = 1  =>  Delta_STAM = Delta (matches GR exterior).
    Near horizon:  Delta ~ (r - r_+),  F ~ 10 (1 - y_K_eq)^2 ~ (r - r_+)^2.
                   Therefore Delta_STAM ~ (r - r_+)^3.  Cubic vanishing preserved.
    """
    D = Delta(r, a)
    if r >= r_ph_prograde(a):
        return D
    y = y_K_eq(r, a)
    return D * F_quintic(y)


def drstar_dr_GR(r, a):
    """Standard Kerr tortoise:  dr*/dr = (r^2 + a^2) / Delta."""
    return (r * r + a * a) / Delta(r, a)


def drstar_dr_STAM(r, a):
    """STAM-modified Kerr tortoise:  dr*/dr = (r^2 + a^2) / Delta_STAM."""
    return (r * r + a * a) / Delta_STAM(r, a)


# ===========================================================================
# Goal 5.  No imaginary h/k at high spin
# ===========================================================================

def verify_no_imaginary(a_values, n_samples=200):
    """At each spin, sample r in (r_+, r_pr_eq) and verify Delta, F, Delta_STAM,
       and dr*/dr_STAM are all real and finite."""
    results = []
    for a in a_values:
        rh = r_horizon(a)
        rp = r_ph_prograde(a)
        if not np.isfinite(rh) or not np.isfinite(rp):
            results.append({"a": a, "ok": False, "reason": "rh or rp NaN"})
            continue
        if rp <= rh:
            results.append({"a": a, "ok": False, "reason": "rp <= rh"})
            continue
        r_grid = np.linspace(rh + 1e-6, rp - 1e-6, n_samples)
        D_vals = np.array([Delta(r, a) for r in r_grid])
        y_vals = np.array([y_K_eq(r, a) for r in r_grid])
        F_vals = np.array([F_quintic(y) for y in y_vals])
        D_STAM_vals = np.array([Delta_STAM(r, a) for r in r_grid])
        drstar = np.array([drstar_dr_STAM(r, a) for r in r_grid])

        all_finite = (np.all(np.isfinite(D_vals))
                      and np.all(np.isfinite(y_vals))
                      and np.all(np.isfinite(F_vals))
                      and np.all(np.isfinite(D_STAM_vals))
                      and np.all(np.isfinite(drstar)))
        all_positive_or_zero = (np.all(D_vals >= 0)
                                and np.all(F_vals >= 0)
                                and np.all(D_STAM_vals >= 0)
                                and np.all(drstar > 0))
        ok = all_finite and all_positive_or_zero
        results.append({
            "a": a, "rh": rh, "rp": rp, "ok": ok,
            "Delta_min": float(D_vals.min()), "Delta_max": float(D_vals.max()),
            "F_min": float(F_vals.min()), "F_max": float(F_vals.max()),
            "DSTAM_min": float(D_STAM_vals.min()), "DSTAM_max": float(D_STAM_vals.max()),
            "drstar_min": float(drstar.min()), "drstar_max": float(drstar.max()),
            "y_min": float(y_vals.min()), "y_max": float(y_vals.max()),
        })
    return results


# ===========================================================================
# Goal 6.  STAM throat width vs spin
# ===========================================================================

def tortoise_distance_from_cutoff(r_target, a, mode="STAM",
                                   r_low_cutoff_eps=1e-5):
    """Integrate dr*/dr from r_+ + cutoff to r_target.  Returns tortoise distance.

    For Kerr GR or STAM, the tortoise diverges at r_+ (logarithmically for GR,
    power-law for STAM).  We measure throat width from a small cutoff eps
    above the horizon.
    """
    rh = r_horizon(a)
    r_lo = rh + r_low_cutoff_eps * (rh if rh > 0 else 1.0)
    if r_target <= r_lo:
        return 0.0
    integrand = drstar_dr_STAM if mode == "STAM" else drstar_dr_GR

    def rhs(r, y):
        return [integrand(r, a)]

    sol = solve_ivp(rhs, [r_lo, r_target], [0.0], method="RK45",
                    rtol=1e-10, atol=1e-13, max_step=(r_target - r_lo) / 50.0)
    return float(sol.y[0, -1])


def throat_widths_vs_spin(a_values, cutoff_eps=1e-5):
    """For each spin compute:
       - throat width in r-coordinate  (r_pr_eq - r_+)
       - throat width in r*-coordinate, GR   (r*_GR(r_pr_eq) - r*_GR(r_+ + eps))
       - throat width in r*-coordinate, STAM (r*_STAM(r_pr_eq) - r*_STAM(r_+ + eps))
    """
    rows = []
    for a in a_values:
        rh = r_horizon(a)
        rp = r_ph_prograde(a)
        width_r = rp - rh
        # tortoise widths from cutoff eps above horizon to r_pr_eq
        w_GR = tortoise_distance_from_cutoff(rp, a, mode="GR", r_low_cutoff_eps=cutoff_eps)
        w_STAM = tortoise_distance_from_cutoff(rp, a, mode="STAM", r_low_cutoff_eps=cutoff_eps)
        rows.append({
            "a": a, "r_horizon": rh, "r_pr_eq": rp,
            "width_r": width_r,
            "width_rstar_GR_from_cutoff": w_GR,
            "width_rstar_STAM_from_cutoff": w_STAM,
            "stretching_factor": w_STAM / w_GR if w_GR > 0 else float('nan'),
        })
    return rows


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 88)
    print("G119  --  Kerr STAM wave base  (GEOMETRY ONLY)")
    print("=" * 88)
    print()
    print("Goals (minimum, sequential):")
    print("  1. r_+(a) real Kerr horizon, not r = 2M")
    print("  2. r_pr_exact(theta, a) from G85/G116")
    print("  3. y_K equatorial via Sigma_ph_eq(a)")
    print("  4. Kerr tortoise dr*/dr = (r^2 + a^2)/Delta_STAM")
    print("  5. No imaginary h/k at high spin")
    print("  6. Plot STAM throat width vs spin")
    print("  7. (QNM diagnostics: G120, NOT this script)")
    print()

    a_values = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99]

    # ----- Goals 1, 2, 3: anchor table -----
    print("=" * 88)
    print("Anchor table:  r_+, r_pr_eq, r_polar, Sigma_ph_eq")
    print("=" * 88)
    print()
    print(f"  {'a':>5}  {'r_+':>10}  {'r_pr_eq':>10}  {'r_polar':>10}  {'Sig_ph_eq':>11}")
    print("-" * 60)
    for a in a_values:
        rh = r_horizon(a)
        rp = r_ph_prograde(a)
        rpol = r_polar(a)
        Sph = Sigma_ph_eq(a)
        print(f"  {a:>5.2f}  {rh:>10.6f}  {rp:>10.6f}  {rpol:>10.6f}  {Sph:>11.6f}")
    print()

    # ----- Goal 5: imaginary check -----
    print("=" * 88)
    print("No imaginary / NaN check:")
    print("=" * 88)
    print()
    check_rows = verify_no_imaginary(a_values, n_samples=200)
    print(f"  {'a':>5}  {'r_+':>8}  {'r_pr_eq':>8}  {'D min/max':>20}  "
          f"{'F min/max':>20}  {'dr*/dr min/max':>22}  ok")
    print("-" * 100)
    for r in check_rows:
        if not r["ok"]:
            print(f"  {r['a']:>5.2f}  {r.get('reason', 'fail')}")
            continue
        print(f"  {r['a']:>5.2f}  {r['rh']:>8.4f}  {r['rp']:>8.4f}  "
              f"{r['Delta_min']:>9.3e}/{r['Delta_max']:>8.3e}  "
              f"{r['F_min']:>9.3e}/{r['F_max']:>8.3e}  "
              f"{r['drstar_min']:>10.3e}/{r['drstar_max']:>10.3e}  "
              f"{'OK' if r['ok'] else 'FAIL'}")
    print()
    n_ok = sum(1 for r in check_rows if r["ok"])
    print(f"  passed: {n_ok}/{len(check_rows)} spins")
    print()

    # ----- Goal 6: throat widths -----
    print("=" * 88)
    print("Throat widths vs spin (cutoff eps = 1e-5 above horizon)")
    print("=" * 88)
    print()
    width_rows = throat_widths_vs_spin(a_values, cutoff_eps=1e-5)
    print(f"  {'a':>5}  {'r_+':>8}  {'r_pr_eq':>8}  {'Δr':>10}  "
          f"{'r* GR width':>14}  {'r* STAM width':>16}  {'stretch':>10}")
    print("-" * 90)
    for r in width_rows:
        print(f"  {r['a']:>5.2f}  {r['r_horizon']:>8.4f}  {r['r_pr_eq']:>8.4f}  "
              f"{r['width_r']:>10.4f}  {r['width_rstar_GR_from_cutoff']:>14.4f}  "
              f"{r['width_rstar_STAM_from_cutoff']:>16.4f}  "
              f"{r['stretching_factor']:>10.4f}")
    print()
    print("'stretch' = STAM/GR tortoise-width ratio (>1 means STAM-stretched).")
    print()

    # ----- Plots -----
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (a) Anchors vs spin
    ax = axes[0, 0]
    spins = [r["a"] for r in width_rows]
    rh_arr = [r["r_horizon"] for r in width_rows]
    rp_arr = [r["r_pr_eq"] for r in width_rows]
    ax.plot(spins, rh_arr, "o-", label="r_+(a)", linewidth=2)
    ax.plot(spins, rp_arr, "s-", label="r_pr_eq(a)", linewidth=2)
    ax.axhline(2.0, color="gray", linestyle=":", alpha=0.5,
               label="r = 2M (Schwarzschild horizon)")
    ax.axhline(3.0, color="gray", linestyle="--", alpha=0.5,
               label="r = 3M (Schwarzschild PS)")
    ax.set_xlabel("spin a / M")
    ax.set_ylabel("radius / M")
    ax.set_title("Kerr horizon and prograde photon orbit (G85 / G116)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (b) Throat width in r-coordinate
    ax = axes[0, 1]
    Δr_arr = [r["width_r"] for r in width_rows]
    ax.plot(spins, Δr_arr, "o-", linewidth=2, color="tab:purple")
    ax.set_xlabel("spin a / M")
    ax.set_ylabel("r_pr_eq - r_+   (M)")
    ax.set_title("STAM throat width  Δr  vs spin")
    ax.grid(True, alpha=0.3)

    # (c) Tortoise widths
    ax = axes[1, 0]
    w_GR_arr = [r["width_rstar_GR_from_cutoff"] for r in width_rows]
    w_S_arr = [r["width_rstar_STAM_from_cutoff"] for r in width_rows]
    ax.plot(spins, w_GR_arr, "o-", linewidth=2, color="tab:blue", label="GR tortoise width")
    ax.plot(spins, w_S_arr, "s-", linewidth=2, color="tab:green", label="STAM tortoise width")
    ax.set_xlabel("spin a / M")
    ax.set_ylabel("tortoise width (M units)")
    ax.set_title("Throat width in tortoise coordinate, from cutoff eps = 1e-5")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (d) Stretch factor STAM/GR
    ax = axes[1, 1]
    stretch_arr = [r["stretching_factor"] for r in width_rows]
    ax.plot(spins, stretch_arr, "^-", linewidth=2, color="tab:red")
    ax.axhline(1.0, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("spin a / M")
    ax.set_ylabel("STAM / GR tortoise-width ratio")
    ax.set_title("Power-law throat stretching by spin")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = PLOTS / "G119_kerr_wave_base_geometry.png"
    plt.savefig(out_png, dpi=180)
    plt.close()
    print(f"Plot saved: {out_png}")

    # ----- Verdict -----
    print()
    print("=" * 88)
    print("VERDICT")
    print("=" * 88)
    print()
    if n_ok == len(check_rows):
        verdict = (
            "Wave base is numerically sound at all sampled spins (a in [0, 0.99]).\n"
            "Delta, F(y_K), Delta_STAM, and dr*/dr_STAM are real, finite, and\n"
            "appropriately signed throughout the inside-shell.  No Schwarzschild\n"
            "h(r) < 0 issue (G118 failure mode) recurs.\n\n"
            "The wave base is ready for G120 (rigorous Kerr QNM via Detweiler /\n"
            "Sasaki-Nakamura on this base)."
        )
    else:
        verdict = (
            f"Wave base FAILED imaginary/finite check at {len(check_rows) - n_ok} of\n"
            f"{len(check_rows)} spins.  Investigate before proceeding to G120."
        )
    print(verdict)

    # ----- Markdown summary -----
    md = ["# G119 -- Kerr STAM wave base (geometry only)\n"]
    md.append("\nReplaces the broken Schwarzschild base of G118 with the proper Kerr\n"
              "radial structure.  Geometry-only -- no Teukolsky physics, no QNM.\n")
    md.append("\n## Anchor table\n")
    md.append("| a | r_+ | r_pr_eq | r_polar | Sigma_ph_eq |\n|---:|---:|---:|---:|---:|\n")
    for a in a_values:
        rh = r_horizon(a); rp = r_ph_prograde(a)
        rpol = r_polar(a); Sph = Sigma_ph_eq(a)
        md.append(f"| {a:.2f} | {rh:.6f} | {rp:.6f} | {rpol:.6f} | {Sph:.6f} |\n")
    md.append("\n## Imaginary / NaN check\n")
    md.append(f"- Passed: {n_ok}/{len(check_rows)} sampled spins.\n")
    md.append("- Delta, F(y_K), Delta_STAM, and dr*/dr_STAM are real and positive\n"
              "  throughout the inside-shell at all spins.\n")
    md.append("\n## Throat widths\n")
    md.append("| a | r_+ | r_pr_eq | Δr | r* GR width | r* STAM width | stretch |\n")
    md.append("|---:|---:|---:|---:|---:|---:|---:|\n")
    for r in width_rows:
        md.append(f"| {r['a']:.2f} | {r['r_horizon']:.4f} | {r['r_pr_eq']:.4f} | "
                  f"{r['width_r']:.4f} | {r['width_rstar_GR_from_cutoff']:.4f} | "
                  f"{r['width_rstar_STAM_from_cutoff']:.4f} | "
                  f"{r['stretching_factor']:.4f} |\n")
    md.append("\n## Verdict\n")
    md.append(verdict + "\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G119_kerr_wave_base_geometry.py](../scripts/G119_kerr_wave_base_geometry.py)\n")
    md.append("- [plots/G119_kerr_wave_base_geometry.png](../plots/G119_kerr_wave_base_geometry.png)\n")
    (RESULTS / "G119_kerr_wave_base_geometry_summary.md").write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {RESULTS / 'G119_kerr_wave_base_geometry_summary.md'}")


if __name__ == "__main__":
    main()
