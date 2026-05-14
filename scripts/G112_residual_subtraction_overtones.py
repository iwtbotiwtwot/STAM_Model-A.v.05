#!/usr/bin/env python3
"""
G112_residual_subtraction_overtones.py

Overtone extraction via residual subtraction:
  Step 1.  Late-time single-mode fit to obtain calibrated n=0 mode.
  Step 2.  Subtract n=0 template from full signal.
  Step 3.  Early-time single-mode fit to residual with Leaver n=1 prior.

This avoids the global two-mode minimum trap (G110 failure) and the
matrix-pencil drifting onto non-QNM transients (G111 failure).

Setup
=====
Same grid + V_exact as G110/G111.  For each ell = 2, 3, 4 extract
omega_n=0 and omega_n=1 for V_GR (calibration) and V_exact (prediction).

Honesty note
============
Time-domain overtone extraction is known to be hard for n >= 1 because the
n=1 amplitude decays with Im(omega) ~ 3x n=0.  Residual subtraction is the
best time-domain approach short of frequency-domain Leaver continued
fraction.  If calibration still fails for n=1, the right tool is Leaver;
that's a separate G113 (frequency-domain) script.

Outputs
=======
    results/G112_residual_subtraction_overtones_summary.md
    plots/G112_residual_subtraction_overtones.png
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
LEAVER_GOLD = {
    (2, 0): complex(0.373672, -0.088962),
    (2, 1): complex(0.346711, -0.273915),
    (3, 0): complex(0.599443, -0.092703),
    (3, 1): complex(0.582644, -0.281298),
    (4, 0): complex(0.809178, -0.094164),
    (4, 1): complex(0.796632, -0.284334),
}


# ====== Metric / V_exact (same as G110/G111) ======

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


def V_geom_grid(r_grid, ell):
    out = np.empty_like(r_grid)
    for i, rv in enumerate(r_grid):
        h_val = h_fn(rv)
        k_val = k_STAM(rv)
        out[i] = h_val * (
            ell * (ell + 1) / rv**2
            - 2.0 * (1.0 - k_val) / rv**2
            - kp_STAM(rv) / (2.0 * rv)
            - k_val * hp_fn(rv) / (2.0 * h_val * rv)
        )
    return out


def analytic_correction(r_grid_STAM):
    out = np.zeros_like(r_grid_STAM)
    for i, rv in enumerate(r_grid_STAM):
        if rv >= 3.0 * M:
            out[i] = 0.0
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


# ====== Fits ======

def one_mode(t, A, om_r, om_i, ph, off):
    tau = t - t[0]
    return A * np.exp(-om_i * tau) * np.cos(om_r * tau + ph) + off


def fit_one_mode(t, sig, t0, t1, om_guess, tight=True):
    mask = (t >= t0) & (t <= t1)
    tf, sf = t[mask], sig[mask]
    if len(tf) < 10:
        return None, None
    A_g = 0.5 * (np.max(sf) - np.min(sf))
    p0 = [A_g, om_guess.real, abs(om_guess.imag), 0.0, 0.0]
    if tight:
        # 30% bounds around initial guess on omega
        om_r_lo = 0.7 * om_guess.real
        om_r_hi = 1.3 * om_guess.real
        om_i_lo = 0.5 * abs(om_guess.imag)
        om_i_hi = 1.5 * abs(om_guess.imag)
    else:
        om_r_lo, om_r_hi = 0.01, 5.0
        om_i_lo, om_i_hi = 0.001, 5.0
    bounds = ([-np.inf, om_r_lo, om_i_lo, -np.pi, -np.inf],
              [np.inf, om_r_hi, om_i_hi, np.pi, np.inf])
    try:
        popt, _ = curve_fit(one_mode, tf, sf, p0=p0, bounds=bounds, maxfev=40000)
        om = complex(popt[1], -popt[2])
        # Reconstruct the template ON THE FULL GRID using popt's amplitude/phase referenced to tf[0]
        return om, popt
    except Exception:
        return None, None


def evaluate_template(t_full, popt, t0_anchor):
    """Reconstruct A exp(-om_i (t-t0)) cos(om_r (t-t0) + ph) + off on full grid."""
    A, om_r, om_i, ph, off = popt
    tau = t_full - t0_anchor
    return A * np.exp(-om_i * tau) * np.cos(om_r * tau + ph) + off


def main():
    print("=" * 88)
    print("G112 -- residual-subtraction overtone extraction")
    print("=" * 88)
    print()

    rs_min, rs_max, N = -500.0, 300.0, 11001
    rstar = np.linspace(rs_min, rs_max, N)
    dr_star = rstar[1] - rstar[0]
    print(f"r* grid: N = {N}, dr* = {dr_star:.4e}")
    r_grid_GR = np.array([r_of_rstar_GR(rs) for rs in rstar])
    r_grid_STAM = build_r_of_rstar_STAM(rstar)
    eps_min = (r_grid_STAM.min() - 2.0 * M) / (2.0 * M)
    print(f"STAM eps_min = {eps_min:.4e}")
    print()

    pulse_center, pulse_sigma = 15.0, 1.5
    psi0 = np.exp(-(rstar - pulse_center)**2 / (2.0 * pulse_sigma**2))
    i_obs = int(np.argmin(np.abs(rstar - 40.0)))
    dt = 0.5 * dr_star
    T_final = 350.0

    results = []
    for ell in [2, 3, 4]:
        print("-" * 88)
        print(f"ell = {ell}")
        print("-" * 88, flush=True)

        V_GR_grid = h_fn(r_grid_GR) * (ell * (ell + 1) / r_grid_GR**2 - 6.0 * M / r_grid_GR**3)
        V_g = V_geom_grid(r_grid_STAM, ell)
        corr = analytic_correction(r_grid_STAM)
        V_exact = V_g + corr

        t_arr, sig_GR = evolve_TD(V_GR_grid, dr_star, T_final, dt, psi0, i_obs)
        _, sig_exact = evolve_TD(V_exact, dr_star, T_final, dt, psi0, i_obs)

        # ---- n=0 from late-time single-mode fit (clean) ----
        t0_late, t1_late = 100.0, 250.0
        om_GR_n0, popt_GR_n0 = fit_one_mode(t_arr, sig_GR, t0_late, t1_late,
                                             LEAVER_GOLD[(ell, 0)], tight=True)
        om_ex_n0, popt_ex_n0 = fit_one_mode(t_arr, sig_exact, t0_late, t1_late,
                                             LEAVER_GOLD[(ell, 0)], tight=True)
        if om_GR_n0 is None or om_ex_n0 is None:
            print(f"  ell={ell}: n=0 fit failed")
            continue
        calib_n0 = abs(om_GR_n0 - LEAVER_GOLD[(ell, 0)]) / abs(LEAVER_GOLD[(ell, 0)])
        shift_n0 = abs(om_ex_n0 - om_GR_n0) / abs(om_GR_n0)
        print(f"  n=0  GR: {om_GR_n0}  Leaver: {LEAVER_GOLD[(ell,0)]}  calib {calib_n0*100:.3f}%")
        print(f"  n=0  exact: {om_ex_n0}  shift {shift_n0*100:.4f}%")

        # ---- Residual: subtract n=0 template, fit residual at early times ----
        # Anchor template at t = t0_late (where fit was done); evaluate on full grid
        tpl_GR_n0 = evaluate_template(t_arr, popt_GR_n0, t0_late)
        tpl_ex_n0 = evaluate_template(t_arr, popt_ex_n0, t0_late)
        resid_GR = sig_GR - tpl_GR_n0
        resid_ex = sig_exact - tpl_ex_n0

        # n=1 fit on residual, early-time window where overtone is still visible
        # n=1 Im(om) ~ 0.27 -> e-folding time ~ 4 M.  Residual will be < n=0
        # in late time only after the subtraction, so use [pulse_arrival, ~50]
        t0_early, t1_early = 30.0, 70.0
        om_GR_n1, popt_GR_n1 = fit_one_mode(t_arr, resid_GR, t0_early, t1_early,
                                             LEAVER_GOLD[(ell, 1)], tight=True)
        om_ex_n1, popt_ex_n1 = fit_one_mode(t_arr, resid_ex, t0_early, t1_early,
                                             LEAVER_GOLD[(ell, 1)], tight=True)
        if om_GR_n1 is None or om_ex_n1 is None:
            print(f"  ell={ell}: n=1 residual fit failed")
            results.append({"ell": ell, "n0_ok": True,
                            "om_GR_n0": om_GR_n0, "om_ex_n0": om_ex_n0,
                            "shift_n0": shift_n0, "calib_n0": calib_n0,
                            "n1_ok": False})
            continue
        calib_n1 = abs(om_GR_n1 - LEAVER_GOLD[(ell, 1)]) / abs(LEAVER_GOLD[(ell, 1)])
        shift_n1 = abs(om_ex_n1 - om_GR_n1) / abs(om_GR_n1)
        print(f"  n=1  GR resid: {om_GR_n1}  Leaver: {LEAVER_GOLD[(ell,1)]}  calib {calib_n1*100:.3f}%")
        print(f"  n=1  exact resid: {om_ex_n1}  shift {shift_n1*100:.4f}%")

        results.append({
            "ell": ell, "n0_ok": True, "n1_ok": True,
            "om_GR_n0": om_GR_n0, "om_ex_n0": om_ex_n0,
            "om_GR_n1": om_GR_n1, "om_ex_n1": om_ex_n1,
            "shift_n0": shift_n0, "shift_n1": shift_n1,
            "calib_n0": calib_n0, "calib_n1": calib_n1,
        })
        print()

    # ----- Summary -----
    print("=" * 88)
    print("Summary")
    print("=" * 88)
    print()
    print(f"  {'ell':>4}  {'n':>2}  {'calib':>8}  "
          f"{'omega_GR':>26}  {'omega_exact':>26}  {'shift':>8}  status")
    print("-" * 100)
    for r in results:
        ell = r["ell"]
        if r.get("n0_ok"):
            ok0 = "OK" if r["calib_n0"] < 0.01 else "POOR_CALIB"
            print(f"  {ell:>4}  {'0':>2}  {r['calib_n0']*100:>7.3f}%  "
                  f"{str(r['om_GR_n0']):>26}  {str(r['om_ex_n0']):>26}  "
                  f"{r['shift_n0']*100:>7.4f}%  {ok0}")
        if r.get("n1_ok"):
            ok1 = "OK" if r["calib_n1"] < 0.05 else "POOR_CALIB"
            print(f"  {ell:>4}  {'1':>2}  {r['calib_n1']*100:>7.3f}%  "
                  f"{str(r['om_GR_n1']):>26}  {str(r['om_ex_n1']):>26}  "
                  f"{r['shift_n1']*100:>7.4f}%  {ok1}")
        else:
            print(f"  {ell:>4}  {'1':>2}  EXTRACTION FAILED")
    print()

    # ----- Plot -----
    if results:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        ax = axes[0]
        ells = [r["ell"] for r in results if r.get("n0_ok")]
        sh0 = [r["shift_n0"]*100 for r in results if r.get("n0_ok")]
        ax.plot(ells, sh0, "o-", color="tab:blue", markersize=10, linewidth=2, label="n=0 (fundamental)")
        ells1 = [r["ell"] for r in results if r.get("n1_ok") and r.get("calib_n1", 1) < 0.05]
        sh1 = [r["shift_n1"]*100 for r in results if r.get("n1_ok") and r.get("calib_n1", 1) < 0.05]
        if ells1:
            ax.plot(ells1, sh1, "s-", color="tab:orange", markersize=10, linewidth=2, label="n=1 (first overtone)")
        ax.set_xlabel("ell")
        ax.set_ylabel("V_exact shift vs GR  (%)")
        ax.set_title("Axial QNM shift by mode (residual-subtraction method)")
        ax.set_xticks([2, 3, 4])
        ax.legend()
        ax.grid(True, alpha=0.3)

        # calibration plot
        ax = axes[1]
        for r in results:
            if r.get("n0_ok"):
                ax.scatter(r["ell"], r["calib_n0"]*100, c="tab:blue", marker="o", s=80,
                           label="n=0 calib" if r["ell"] == 2 else None)
            if r.get("n1_ok"):
                ax.scatter(r["ell"], r["calib_n1"]*100, c="tab:orange", marker="s", s=80,
                           label="n=1 calib" if r["ell"] == 2 else None)
        ax.axhline(1.0, color="gray", linestyle="--", alpha=0.5, label="1% threshold")
        ax.set_xlabel("ell")
        ax.set_ylabel("Calibration error  (%)")
        ax.set_title("Schwarzschild calibration error")
        ax.set_yscale("log")
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        out_png = PLOTS / "G112_residual_subtraction_overtones.png"
        plt.savefig(out_png, dpi=180)
        plt.close()
        print(f"Plot saved: {out_png}")

    # ----- Markdown -----
    md = ["# G112 - Higher-mode QNM via residual subtraction\n"]
    md.append("\nReplacement for G110 (two-mode fit failed) and G111 (matrix pencil drifted).\n")
    md.append("Standard Prony-type residual subtraction: fit n=0 at late times,\n"
              "subtract template, fit residual at early times for n=1.\n")
    md.append("\n## Results\n")
    md.append("| ell | n | calib | omega_GR | Leaver | omega_exact | shift | status |\n")
    md.append("|---:|---:|---:|---|---|---|---:|---|\n")
    for r in results:
        ell = r["ell"]
        if r.get("n0_ok"):
            ok0 = "OK" if r["calib_n0"] < 0.01 else "POOR_CALIB"
            md.append(f"| {ell} | 0 | {r['calib_n0']*100:.3f}% | {r['om_GR_n0']} | "
                      f"{LEAVER_GOLD[(ell,0)]} | {r['om_ex_n0']} | "
                      f"{r['shift_n0']*100:.4f}% | {ok0} |\n")
        if r.get("n1_ok"):
            ok1 = "OK" if r["calib_n1"] < 0.05 else "POOR_CALIB"
            md.append(f"| {ell} | 1 | {r['calib_n1']*100:.3f}% | {r['om_GR_n1']} | "
                      f"{LEAVER_GOLD[(ell,1)]} | {r['om_ex_n1']} | "
                      f"{r['shift_n1']*100:.4f}% | {ok1} |\n")
        else:
            md.append(f"| {ell} | 1 | --- | EXTRACTION FAILED | {LEAVER_GOLD[(ell,1)]} | --- | --- | FAILED |\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G112_residual_subtraction_overtones.py](../scripts/G112_residual_subtraction_overtones.py)\n")
    md.append("- [plots/G112_residual_subtraction_overtones.png](../plots/G112_residual_subtraction_overtones.png)\n")
    out_md = RESULTS / "G112_residual_subtraction_overtones_summary.md"
    out_md.write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {out_md}")


if __name__ == "__main__":
    main()
