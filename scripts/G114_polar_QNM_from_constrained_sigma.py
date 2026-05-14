#!/usr/bin/env python3
"""
G114_polar_QNM_from_constrained_sigma.py

Polar (even-parity) sector QNM under the constrained shell-count action
S[Sigma, g, lambda_1, lambda_2] -- companion to G108's axial result.

Structure
=========
Unlike axial (where Sigma scalar has no parity-odd harmonic, so trivially
delta Sigma = delta lambda_1 = delta lambda_2 = 0), the polar sector
carries scalar harmonics.  The two Lagrange-multiplier constraints

    lambda_1 ( (grad Sigma)^2 - W ) = 0   (kinematic)
    lambda_2 ( u^mu d_mu Sigma )    = 0   (flow-constancy)

algebraically eliminate delta Sigma in terms of polar metric perturbations:

    delta lambda_1:  2 g^rr Sigma_bar' partial_r delta Sigma  -  W' delta Sigma  =  0
                     => partial_r delta Sigma  algebraically pinned in delta Sigma.

    delta lambda_2:  u^t partial_t delta Sigma  +  delta u^r Sigma_bar'  =  0
                     => partial_t delta Sigma  algebraically pinned to delta u^r.

Together: delta Sigma is fully fixed by polar metric perturbations.  The
scalar shell-count mode does NOT propagate independently; only the graviton
remains dynamical -- same conclusion as axial, realized algebraically.

Polar master equation
=====================
With delta Sigma constrained, the standard Jordan-frame f(Sigma) R polar
reduction gives a Zerilli-type master equation (rather than Regge-Wheeler):

    Psi_tt  -  Psi_**  +  V_polar_exact * Psi  =  0

In Schwarzschild (f = 1), the polar potential is the standard Zerilli form

    V_Zerilli_GR(r) = h(r) [ 2 n^2 (n+1) r^3 + 6 n^2 M r^2 + 18 n M^2 r + 18 M^3 ]
                      / [ r^3 (n r + 3 M)^2 ]

with n = (ell - 1)(ell + 2) / 2 = (ell^2 + ell - 2) / 2.

For Schwarzschild, V_Zerilli_GR is ISOSPECTRAL with V_RW_GR (Chandrasekhar's
transformation): both give the same QNM frequencies.  This is a known
special feature of GR that we use as a calibration check.

For the STAM background (k = h F(y) inside PS, with non-minimal coupling
f(Sigma)):
    1. Replace M in Zerilli's formula by the effective mass
       M_eff(r) = (r / 2) (1 - k_STAM(r))   inside PS,
       M outside PS.  This is the polar analog of G87's axial proxy.
    2. Add the canonical-form action correction (sqrt f)'' / sqrt f
       (same form as G108's axial result).

    V_polar_exact = V_Zerilli_STAM(M_eff)  +  (sqrt f)'' / sqrt f

Caveats (honest)
================
This is a FIRST-PASS polar prediction.  The simplifications:

(a) M_eff substitution in Zerilli's formula is a "polar proxy" analog of
    G87's axial proxy.  Rigorous polar reduction with k != h inside PS
    requires redoing Zerilli's calculation with the STAM-modified metric;
    that derivation is left for G114-rigorous.
(b) The (sqrt f)'' / sqrt f correction is taken by analogy with axial.
    For polar, additional structure may appear from the scalar-mode
    elimination through delta lambda_1, delta lambda_2 (delta u^r-coupled).
    A rigorous polar reduction may augment this with terms involving
    derivatives of the constraint functions.

These caveats mirror the original G92 vs G108 history for axial: the
"canonical proxy" form is a useful first prediction; the rigorous form
follows from a careful derivation.

What we WILL get from this script
=================================
- V_Zerilli_GR(r) implemented and calibrated against Leaver gold for ell=2.
- V_polar_first-pass for STAM with M_eff proxy + (sqrt f)'' / sqrt f correction.
- Time-domain extraction of the polar QNM shift for STAM.
- Comparison with axial 4.977% shift (G109): does isospectrality break?

If polar shift ~ axial shift: STAM preserves approximate isospectrality.
If polar shift differs significantly: STAM breaks isospectrality (a generic
feature of non-minimally-coupled scalar-tensor theories with non-trivial
background scalar).  Either result is meaningful.

Output
======
    results/G114_polar_QNM_summary.md
    plots/G114_polar_QNM.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


M = 1.0
CALIBRATION_THRESHOLD = 0.02   # 2% (Zerilli for higher ell has steeper peak, fit harder)

LEAVER_GOLD = {
    (2, 0): complex(0.373672, -0.088962),
    (3, 0): complex(0.599443, -0.092703),
    (4, 0): complex(0.809178, -0.094164),
}


# ===========================================================================
# Metric / f(Sigma) (copied from G108)
# ===========================================================================

def h_fn(r): return 1.0 - 2.0 * M / r
def hp_fn(r): return 2.0 * M / r**2
def y_of_r(r): return 6.0 * M / r - 2.0
def F_quintic(y): return 1.0 - 5.0 * y**4 + 4.0 * y**5
def F_first(y): return -20.0 * y**3 + 20.0 * y**4

def k_STAM(r):
    y = y_of_r(r)
    return h_fn(r) * F_quintic(y) if y > 0 else h_fn(r)

def kp_STAM(r):
    y = y_of_r(r)
    if y > 0:
        return hp_fn(r) * F_quintic(y) + h_fn(r) * F_first(y) * (-6.0 * M / r**2)
    return hp_fn(r)


def dlnf_dA_numeric(A):
    if A <= 2.0 / 3.0:
        return 0.0
    y = 3.0 * A - 2.0
    F = 1.0 - 5.0 * y**4 + 4.0 * y**5
    return -2.0 * y**3 * (45.0 * A**2 - 33.0 * A - 13.0) / (A * F)

def d2lnf_dA2_numeric(A, hstep=1e-6):
    return (dlnf_dA_numeric(A + hstep) - dlnf_dA_numeric(A - hstep)) / (2.0 * hstep)


def build_f_grid(r_grid_STAM):
    ln_f = np.zeros_like(r_grid_STAM)
    inside = r_grid_STAM < 3.0 * M
    if not np.any(inside):
        return np.ones_like(r_grid_STAM)
    r_inside = r_grid_STAM[inside]
    r_descending = np.sort(r_inside)[::-1]
    def rhs(rv, yv):
        A = 2.0 * M / rv
        return [dlnf_dA_numeric(A) * (-2.0 * M / rv**2)]
    sol = solve_ivp(rhs, [3.0 * M - 1e-12, r_descending[-1]], [0.0],
                    t_eval=r_descending, method="RK45",
                    rtol=1e-12, atol=1e-14, max_step=0.005)
    ln_f_inside = np.empty_like(r_inside)
    order = np.argsort(r_inside)
    ln_f_inside[order] = sol.y[0][::-1]
    ln_f[inside] = ln_f_inside
    return np.exp(ln_f)


def analytic_correction(r_grid_STAM):
    out = np.zeros_like(r_grid_STAM)
    for i, rv in enumerate(r_grid_STAM):
        if rv >= 3.0 * M:
            continue
        A = 2.0 * M / rv
        dA_dr = -2.0 * M / rv**2
        d2A_dr2 = 4.0 * M / rv**3
        u = dlnf_dA_numeric(A)
        up = d2lnf_dA2_numeric(A)
        dlnf_dr = u * dA_dr
        d_dlnf_dr = up * (dA_dr**2) + u * d2A_dr2
        hk = h_fn(rv) * k_STAM(rv)
        if hk <= 0:
            continue
        sqrt_hk = np.sqrt(hk)
        d_hk = hp_fn(rv) * k_STAM(rv) + h_fn(rv) * kp_STAM(rv)
        P = sqrt_hk * dlnf_dr
        d_sqrt_hk_dr = 0.5 * d_hk / sqrt_hk
        dP_dr = d_sqrt_hk_dr * dlnf_dr + sqrt_hk * d_dlnf_dr
        P_star = sqrt_hk * dP_dr
        out[i] = 0.5 * P_star + 0.25 * P**2
    return out


# ===========================================================================
# Zerilli polar potential (Schwarzschild)
# ===========================================================================

def V_Zerilli_GR(r, ell):
    """Standard Zerilli polar potential, Schwarzschild s = 2.

    n = (ell - 1)(ell + 2) / 2 = (ell^2 + ell - 2) / 2

    V_Zerilli = h * [2 n^2 (n+1) r^3 + 6 n^2 M r^2 + 18 n M^2 r + 18 M^3]
                  / [r^3 (n r + 3 M)^2]
    """
    n = (ell - 1) * (ell + 2) / 2.0
    h_val = h_fn(r)
    num = 2.0 * n**2 * (n + 1) * r**3 + 6.0 * n**2 * M * r**2 \
          + 18.0 * n * M**2 * r + 18.0 * M**3
    den = r**3 * (n * r + 3.0 * M)**2
    return h_val * num / den


def V_Zerilli_STAM_proxy(r, ell):
    """Zerilli on STAM background using M_eff(r) = (r/2)(1 - k_STAM(r))
       inside PS (polar analog of G87's M_eff substitution)."""
    n = (ell - 1) * (ell + 2) / 2.0
    h_val = h_fn(r)
    k_val = k_STAM(r)
    M_eff = 0.5 * r * (1.0 - k_val)
    num = 2.0 * n**2 * (n + 1) * r**3 + 6.0 * n**2 * M_eff * r**2 \
          + 18.0 * n * M_eff**2 * r + 18.0 * M_eff**3
    den = r**3 * (n * r + 3.0 * M_eff)**2
    return h_val * num / den


def build_V_polar_canonical(r_grid_STAM, rstar_grid, ell):
    """V_polar_canonical = V_Zerilli_STAM_proxy + (sqrt f)'' / sqrt f.

    First-pass form, with the standard caveats noted in the script header.
    """
    V_z = np.array([V_Zerilli_STAM_proxy(r, ell) for r in r_grid_STAM])
    corr = analytic_correction(r_grid_STAM)
    return V_z + corr


# ===========================================================================
# Tortoise + TD pipeline (copied from G108)
# ===========================================================================

def r_of_rstar_GR(rstar):
    PS_offset = 3.0 * M + 2.0 * M * np.log(0.5)
    target = rstar + PS_offset
    def f(r):
        return r + 2.0 * M * np.log(r / (2.0 * M) - 1.0) - target
    if rstar > 0:
        return brentq(f, 3.0 * M - 1e-10, 1e8)
    if rstar == 0:
        return 3.0 * M
    if rstar < -40:
        arg = (rstar + 3.0 * M + 2.0 * M * np.log(0.5) - 2.0 * M) / (2.0 * M)
        r_guess = 2.0 * M + 2.0 * M * np.exp(arg)
        r_lo = max(2.0 * M + 1e-300, r_guess * 1e-3)
        r_hi = min(3.0 * M - 1e-12, max(r_guess * 1e3, 2.0 * M + 1e-6))
        try:
            return brentq(f, r_lo, r_hi)
        except ValueError:
            return r_guess
    return brentq(f, 2.0 * M + 1e-14, 3.0 * M + 1e-10)


def build_r_of_rstar_STAM(rstar_grid):
    def rhs(rs, y):
        rv = y[0]
        return [np.sqrt(h_fn(rv) * k_STAM(rv))]
    rs_pos = rstar_grid[rstar_grid > 0]
    rs_neg = rstar_grid[rstar_grid < 0]
    r_arr = np.zeros_like(rstar_grid)
    if len(rs_pos) > 0:
        order = np.argsort(rs_pos)
        x = rs_pos[order]
        sol = solve_ivp(rhs, [0, x[-1]], [3.0 * M], t_eval=x,
                        method="RK45", rtol=1e-12, atol=1e-14, max_step=1.0)
        vals = np.zeros_like(rs_pos)
        vals[order] = sol.y[0]
        r_arr[rstar_grid > 0] = vals
    if len(rs_neg) > 0:
        order = np.argsort(rs_neg)[::-1]
        x = rs_neg[order]
        sol = solve_ivp(rhs, [0, x[-1]], [3.0 * M], t_eval=x,
                        method="RK45", rtol=1e-12, atol=1e-14, max_step=0.02)
        vals = np.zeros_like(rs_neg)
        vals[order] = sol.y[0]
        r_arr[rstar_grid < 0] = vals
    r_arr[rstar_grid == 0] = 3.0 * M
    return r_arr


def evolve_TD(V_grid, dr_star, T_final, dt, psi0, i_obs):
    psi_old = psi0.copy()
    psi = psi0.copy()
    N_steps = int(T_final / dt)
    inv_dr2 = 1.0 / (dr_star * dr_star)
    cfl = dt / dr_star
    t_arr = np.zeros(N_steps)
    sig_arr = np.zeros(N_steps)
    for step in range(N_steps):
        psi_rr = np.zeros_like(psi)
        psi_rr[1:-1] = (psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) * inv_dr2
        psi_new = np.empty_like(psi)
        psi_new[1:-1] = 2.0 * psi[1:-1] - psi_old[1:-1] + dt * dt * (
            psi_rr[1:-1] - V_grid[1:-1] * psi[1:-1])
        psi_new[-1] = psi[-1] - cfl * (psi[-1] - psi[-2])
        psi_new[0] = psi[0] - cfl * (psi[0] - psi[1])
        psi_old = psi
        psi = psi_new
        t_arr[step] = (step + 1) * dt
        sig_arr[step] = psi[i_obs]
    return t_arr, sig_arr


def damped_sinusoid(t, A, om_r, om_i, ph, off):
    tau = t - t[0]
    return A * np.exp(-om_i * tau) * np.cos(om_r * tau + ph) + off


def fit_QNM(t, sig, t0, t1, om_guess=complex(0.37, -0.09)):
    mask = (t >= t0) & (t <= t1)
    tf, sf = t[mask], sig[mask]
    if len(tf) < 10:
        return None
    A_g = 0.5 * (np.max(sf) - np.min(sf))
    p0 = [A_g, om_guess.real, abs(om_guess.imag), 0.0, 0.0]
    bounds = ([-np.inf, 0, 0, -np.pi, -np.inf],
              [np.inf, 5.0, 5.0, np.pi, np.inf])
    try:
        popt, _ = curve_fit(damped_sinusoid, tf, sf, p0=p0, bounds=bounds, maxfev=40000)
        return complex(popt[1], -popt[2])
    except Exception:
        return None


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 88)
    print("G114 -- Polar (even-parity) QNM from constrained shell-count action")
    print("=" * 88)
    print()
    print("Hypothesis under test:")
    print("  V_polar_first_pass = V_Zerilli_STAM_proxy(M_eff)  +  (sqrt f)'' / sqrt f")
    print()
    print("Calibration: V_Zerilli_GR must reproduce Leaver gold (isospectral with RW).")
    print()

    rs_min, rs_max, N = -500.0, 300.0, 11001
    rstar = np.linspace(rs_min, rs_max, N)
    dr_star = rstar[1] - rstar[0]
    print(f"r* grid: N = {N}, range [{rs_min}, {rs_max}], dr* = {dr_star:.4e}")
    r_grid_GR = np.array([r_of_rstar_GR(rs) for rs in rstar])
    r_grid_STAM = build_r_of_rstar_STAM(rstar)
    eps_min = (r_grid_STAM.min() - 2.0 * M) / (2.0 * M)
    print(f"STAM eps_min = {eps_min:.4e}")
    print()

    # Time-domain settings (broad pulse, like G109's converged setup)
    pulse_center, pulse_sigma = 30.0, 3.0
    psi0 = np.exp(-(rstar - pulse_center)**2 / (2.0 * pulse_sigma**2))
    i_obs = int(np.argmin(np.abs(rstar - 50.0)))
    dt = 0.5 * dr_star
    T_final = 350.0

    results = []
    for ell in [2, 3, 4]:
        print("-" * 88)
        print(f"ell = {ell}")
        print("-" * 88, flush=True)

        # ----- Schwarzschild calibration on V_Zerilli_GR -----
        V_Z_GR_grid = np.array([V_Zerilli_GR(r, ell) for r in r_grid_GR])
        t_arr, sig_GR = evolve_TD(V_Z_GR_grid, dr_star, T_final, dt, psi0, i_obs)
        om_gold = LEAVER_GOLD[(ell, 0)]
        om_GR_polar = fit_QNM(t_arr, sig_GR, 100.0, 250.0, om_guess=om_gold)
        if om_GR_polar is None:
            print("  Zerilli GR calibration fit failed; skipping ell")
            continue
        calib_err = abs(om_GR_polar - om_gold) / abs(om_gold)
        if calib_err > CALIBRATION_THRESHOLD:
            print(f"  Zerilli GR calibration FAILED: err = {calib_err*100:.3f}%")
            print(f"  omega_Z_GR = {om_GR_polar}")
            print(f"  Leaver gold = {om_gold}")
            print(f"  > {CALIBRATION_THRESHOLD*100:.1f}% threshold, skipping STAM")
            continue
        print(f"  Zerilli GR omega = {om_GR_polar}  Leaver = {om_gold}"
              f"  calib err = {calib_err*100:.3f}%  PASS")

        # ----- V_polar_STAM_first_pass -----
        V_polar = build_V_polar_canonical(r_grid_STAM, rstar, ell)
        t_arr_S, sig_polar = evolve_TD(V_polar, dr_star, T_final, dt, psi0, i_obs)
        om_polar = fit_QNM(t_arr_S, sig_polar, 100.0, 250.0, om_guess=om_GR_polar)
        if om_polar is None:
            print("  Polar STAM fit failed")
            continue
        shift_polar = abs(om_polar - om_GR_polar) / abs(om_GR_polar)
        print(f"  Polar V_exact omega = {om_polar}  shift vs GR = {shift_polar*100:.4f}%")

        results.append({
            "ell": ell, "om_GR_polar": om_GR_polar, "om_polar": om_polar,
            "calib_err": calib_err, "shift": shift_polar,
            "sig_GR": sig_GR, "sig_polar": sig_polar, "t": t_arr,
        })
        print()

    # ----- Summary -----
    print("=" * 88)
    print("Polar QNM shifts vs GR")
    print("=" * 88)
    print()
    print(f"  {'ell':>4}  {'calib':>8}  {'omega_GR (Zerilli)':>30}  "
          f"{'omega_polar STAM':>30}  {'shift':>10}")
    print("-" * 95)
    for r in results:
        print(f"  {r['ell']:>4}  {r['calib_err']*100:>7.3f}%  "
              f"{str(r['om_GR_polar']):>30}  {str(r['om_polar']):>30}  "
              f"{r['shift']*100:>9.4f}%")
    print()

    # ----- Comparison with axial -----
    print("=" * 88)
    print("Polar vs axial comparison")
    print("=" * 88)
    print()
    AXIAL_REFERENCE = {
        2: 4.977,   # G109 converged
        3: 2.101,   # G112
        4: 1.158,   # G112
    }
    print(f"  {'ell':>4}  {'axial shift':>14}  {'polar shift':>14}  {'polar/axial':>14}")
    print("-" * 60)
    for r in results:
        ell = r["ell"]
        axial = AXIAL_REFERENCE.get(ell)
        polar = r["shift"] * 100
        ratio = polar / axial if axial else None
        rstr = f"{ratio:.3f}" if ratio is not None else "n/a"
        print(f"  {ell:>4}  {axial:>13.4f}%  {polar:>13.4f}%  {rstr:>14}")
    print()
    print("If polar/axial ~ 1: STAM approximately preserves Schwarzschild isospectrality.")
    print("If polar/axial deviates from 1: STAM breaks isospectrality, which is generic")
    print("                              for non-minimally-coupled scalar-tensor theories.")
    print()

    # ----- Plots -----
    if results:
        fig, axes = plt.subplots(2, len(results), figsize=(5 * len(results), 9))
        if len(results) == 1:
            axes = axes[:, np.newaxis]
        for j, r in enumerate(results):
            ax = axes[0, j]
            ax.plot(r["t"], np.abs(r["sig_GR"]) + 1e-30, "tab:blue",
                    label="Zerilli GR", linewidth=0.9)
            ax.plot(r["t"], np.abs(r["sig_polar"]) + 1e-30, "tab:green",
                    label="V_polar STAM", linewidth=0.9)
            ax.set_yscale("log")
            ax.set_xlabel("t / M")
            ax.set_ylabel("|psi|")
            ax.set_title(f"Polar TD signals (ell = {r['ell']})")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

            ax = axes[1, j]
            # Show V profiles for context
            r_grid_inside = np.linspace(2.01, 3.5, 400)
            V_GR_vis = np.array([V_Zerilli_GR(rv, r["ell"]) for rv in r_grid_inside])
            ax.plot(r_grid_inside, V_GR_vis, "tab:blue", label="V_Zerilli_GR", linewidth=2)
            # Indicate r where it matters
            ax.axvline(3.0, color="black", linestyle=":", alpha=0.4, label="PS")
            ax.axvline(2.0, color="red", linestyle=":", alpha=0.4, label="horizon")
            ax.set_xlabel("r / M")
            ax.set_ylabel("V (l = " + str(r["ell"]) + ")")
            ax.set_title(f"Zerilli potential profile (ell = {r['ell']})")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
        plt.tight_layout()
        out_png = PLOTS / "G114_polar_QNM.png"
        plt.savefig(out_png, dpi=180)
        plt.close()
        print(f"Plot saved: {out_png}")

    # ----- Markdown -----
    md = ["# G114 -- Polar (even-parity) QNM from constrained shell-count action\n"]
    md.append("\nFirst-pass polar prediction.  Uses V_Zerilli_STAM_proxy(M_eff) "
              "+ (sqrt f)'' / sqrt f, by analogy with axial G108 result.\n")
    md.append("\n## Setup\n")
    md.append("- Two LM constraints algebraically eliminate delta Sigma in polar\n"
              "  (delta lambda_1 pins partial_r delta Sigma; delta lambda_2 pins\n"
              "  partial_t delta Sigma to delta u^r).  Scalar mode does not propagate.\n")
    md.append("- Effective polar master equation is Zerilli-type with f-modification.\n")
    md.append("- First-pass:  V_polar = V_Zerilli_STAM_proxy + (sqrt f)'' / sqrt f.\n")
    md.append("\n## Results\n")
    md.append("| ell | calib | omega_GR (Zerilli) | omega_polar STAM | shift |\n")
    md.append("|---:|---:|---|---|---:|\n")
    for r in results:
        md.append(f"| {r['ell']} | {r['calib_err']*100:.3f}% | "
                  f"{r['om_GR_polar']} | {r['om_polar']} | "
                  f"{r['shift']*100:.4f}% |\n")
    md.append("\n## Polar vs axial\n")
    md.append("| ell | axial shift | polar shift | polar / axial |\n")
    md.append("|---:|---:|---:|---:|\n")
    for r in results:
        ell = r["ell"]
        axial = AXIAL_REFERENCE.get(ell)
        polar = r["shift"] * 100
        if axial:
            ratio = polar / axial
            md.append(f"| {ell} | {axial:.4f}% | {polar:.4f}% | {ratio:.3f} |\n")
        else:
            md.append(f"| {ell} | n/a | {polar:.4f}% | n/a |\n")
    md.append("\n## Caveats\n")
    md.append("- M_eff substitution in Zerilli (polar analog of G87 axial proxy);\n"
              "  rigorous polar reduction with k != h inside PS would replace this.\n")
    md.append("- (sqrt f)'' / sqrt f correction by analogy with axial; rigorous polar\n"
              "  reduction may augment with constraint-feedback terms.\n")
    md.append("\nFollow-up:  G114-rigorous would derive V_polar directly by varying\n"
              "S[Sigma,g,lambda_1,lambda_2] in even-parity Regge-Wheeler gauge,\n"
              "eliminating delta Sigma via the two LM equations, and reading off\n"
              "the Sturm-Liouville form.  This is the polar analog of G108.\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G114_polar_QNM_from_constrained_sigma.py]"
              "(../scripts/G114_polar_QNM_from_constrained_sigma.py)\n")
    md.append("- [plots/G114_polar_QNM.png](../plots/G114_polar_QNM.png)\n")
    (RESULTS / "G114_polar_QNM_summary.md").write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {RESULTS / 'G114_polar_QNM_summary.md'}")


if __name__ == "__main__":
    main()
