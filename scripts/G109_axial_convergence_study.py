#!/usr/bin/env python3
"""
G109_axial_convergence_study.py

Convergence study for the rigorous axial QNM under V_exact (analytic
(sqrt f)'' / sqrt f, derived in G108).

Question
========
G108 found omega_exact -> shift = 4.977% on a (N=8001, rs_min=-300) grid,
while G94's robustness sweep on the G92 finite-difference correction gave
mean 5.384% +/- 0.118%.  Both are consistent with "around 5%".

This script tests whether omega_exact converges to a single number as we
deepen the throat (rs_min) at roughly fixed dr_star.  The framework's
predicted QNM shift is the converged value.

Setup
=====
We sweep (rs_min, N) pairs holding dr_star ~ 0.075 fixed and rs_max = 300.
For each pair:
    - Build STAM r(r*) on the deeper grid
    - Compute f(Sigma) by ODE on the new grid
    - Build V_geom + analytic action correction (G108 routine)
    - Schwarzschild calibration (must pass < 1%)
    - Time-domain QNM extraction
    - Record omega and shift

Outputs
=======
    results/G109_axial_convergence_summary.md
    plots/G109_axial_convergence.png

Decision rule
=============
If omega_exact converges to within +/- 0.1 percentage points across the
two deepest grids, the converged value IS the framework's locked QNM
prediction.  Otherwise we report the trend and flag that deeper throats
are still needed.
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


# ===========================================================================
# Setup (matches G108)
# ===========================================================================

M = 1.0
ELL = 2
LEAVER_GOLD = complex(0.373672, -0.088962)
CALIBRATION_THRESHOLD = 0.01


def h_fn(r): return 1.0 - 2.0 * M / r
def hp_fn(r): return 2.0 * M / r**2
def A_of_r(r): return 2.0 * M / r
def y_of_r(r): return 6.0 * M / r - 2.0
def F_quintic(y): return 1.0 - 5.0 * y**4 + 4.0 * y**5
def F_first(y): return -20.0 * y**3 + 20.0 * y**4

def k_GR(r): return h_fn(r)
def kp_GR(r): return hp_fn(r)

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

def d2lnf_dA2_numeric(A, h=1e-6):
    return (dlnf_dA_numeric(A + h) - dlnf_dA_numeric(A - h)) / (2.0 * h)


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


def V_geom_grid(r_grid):
    out = np.empty_like(r_grid)
    for i, rv in enumerate(r_grid):
        h_val = h_fn(rv)
        k_val = k_STAM(rv)
        out[i] = h_val * (
            ELL * (ELL + 1) / rv**2
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
            out[i] = 0.0
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


def damped_sinusoid(t, A, omega_r, omega_i, phi, off):
    return A * np.exp(-omega_i * (t - t[0])) * np.cos(omega_r * (t - t[0]) + phi) + off


def fit_QNM(t, sig, t_fit_start, t_fit_end):
    mask = (t >= t_fit_start) & (t <= t_fit_end)
    tf, sf = t[mask], sig[mask]
    if len(tf) < 10:
        return None
    A_g = 0.5 * (np.max(sf) - np.min(sf))
    p0 = [A_g, 0.37, 0.09, 0.0, 0.0]
    bounds = ([-np.inf, 0, 0, -np.pi, -np.inf],
              [np.inf, 2.0, 2.0, np.pi, np.inf])
    try:
        popt, _ = curve_fit(damped_sinusoid, tf, sf, p0=p0,
                            bounds=bounds, maxfev=30000)
        return complex(popt[1], -popt[2])
    except Exception:
        return None


# ===========================================================================
# Sweep
# ===========================================================================

def run_one(rs_min, N, rs_max=300.0, T_final=350.0,
            t_fit_start=100.0, t_fit_end=250.0):
    """Run one grid point of the convergence sweep.  Return dict of results."""
    rstar = np.linspace(rs_min, rs_max, N)
    dr_star = rstar[1] - rstar[0]
    print(f"  grid: N={N}, rs in [{rs_min}, {rs_max}], dr*={dr_star:.4e}",
          flush=True)

    r_grid_GR = np.array([r_of_rstar_GR(rs) for rs in rstar])
    r_grid_STAM = build_r_of_rstar_STAM(rstar)
    eps_min = (r_grid_STAM.min() - 2.0 * M) / (2.0 * M)
    print(f"  STAM eps_min = {eps_min:.4e}", flush=True)

    f_grid = build_f_grid(r_grid_STAM)
    V_g = V_geom_grid(r_grid_STAM)
    corr = analytic_correction(r_grid_STAM)
    V_exact = V_g + corr

    V_GR_grid = h_fn(r_grid_GR) * (ELL * (ELL + 1) / r_grid_GR**2
                                    - 6.0 * M / r_grid_GR**3)

    # Schwarzschild calibration
    psi0 = np.exp(-(rstar - 30.0)**2 / (2.0 * 3.0**2))
    i_obs = int(np.argmin(np.abs(rstar - 50.0)))
    dt = 0.5 * dr_star
    t_arr, sig_GR = evolve_TD(V_GR_grid, dr_star, T_final, dt, psi0, i_obs)
    omega_GR = fit_QNM(t_arr, sig_GR, t_fit_start, t_fit_end)
    if omega_GR is None:
        return {"converged": False, "reason": "GR fit failed"}
    calib_err = abs(omega_GR - LEAVER_GOLD) / abs(LEAVER_GOLD)
    if calib_err > CALIBRATION_THRESHOLD:
        return {"converged": False, "reason": f"calibration {calib_err*100:.2f}%",
                "omega_GR": omega_GR, "calib_err": calib_err}

    # V_exact run
    t_arr_e, sig_e = evolve_TD(V_exact, dr_star, T_final, dt, psi0, i_obs)
    omega_exact = fit_QNM(t_arr_e, sig_e, t_fit_start, t_fit_end)
    if omega_exact is None:
        return {"converged": False, "reason": "exact fit failed",
                "omega_GR": omega_GR, "calib_err": calib_err}
    shift = abs(omega_exact - omega_GR) / abs(omega_GR)

    return {
        "converged": True,
        "N": N, "rs_min": rs_min, "rs_max": rs_max,
        "dr_star": dr_star, "eps_min": eps_min,
        "f_max": float(f_grid.max()),
        "corr_max": float(np.max(np.abs(corr))),
        "omega_GR": omega_GR, "calib_err": calib_err,
        "omega_exact": omega_exact, "shift": shift,
    }


def main():
    print("=" * 88)
    print("G109  --  convergence study for analytic V_exact axial QNM")
    print("=" * 88)
    print()

    # Sweep: deepen rs_min, scale N to keep dr_star ~ 0.075
    sweep = [
        (-300.0,  8001),
        (-500.0, 11001),
        (-1000.0, 18001),
        (-2000.0, 32001),
    ]

    results = []
    for (rs_min, N) in sweep:
        print(f"--- rs_min = {rs_min}, N = {N} ---", flush=True)
        res = run_one(rs_min, N)
        results.append(res)
        if not res.get("converged", False):
            print(f"  ABORT: {res.get('reason')}")
            continue
        print(f"  omega_GR = {res['omega_GR']}  (calib err {res['calib_err']*100:.3f}%)")
        print(f"  omega_exact = {res['omega_exact']}")
        print(f"  shift = {res['shift']*100:.4f}%")
        print(f"  f_max = {res['f_max']:.3e}, max|correction| = {res['corr_max']:.3e}")
        print()

    # Convergence summary
    print("=" * 88)
    print("Convergence summary")
    print("=" * 88)
    print()
    print(f"  {'rs_min':>8}  {'N':>6}  {'eps_min':>10}  {'f_max':>10}  "
          f"{'calib %':>9}  {'omega_exact':>32}  {'shift %':>10}")
    print("-" * 100)
    for r in results:
        if not r.get("converged", False):
            print(f"  {r.get('rs_min', '?'):>8}  ABORT  ({r.get('reason')})")
            continue
        print(f"  {r['rs_min']:>8.0f}  {r['N']:>6d}  {r['eps_min']:>10.2e}  "
              f"{r['f_max']:>10.2e}  {r['calib_err']*100:>9.4f}  "
              f"{str(r['omega_exact']):>32}  {r['shift']*100:>10.4f}")
    print()

    # Decide
    converged_rows = [r for r in results if r.get("converged", False)]
    if len(converged_rows) >= 2:
        last_two = converged_rows[-2:]
        delta = abs(last_two[1]["shift"] - last_two[0]["shift"]) * 100
        print(f"  |shift_deepest - shift_second_deepest| = {delta:.4f} percentage points")
        if delta < 0.1:
            verdict = (
                f"CONVERGED.  Framework's rigorous axial QNM shift = "
                f"{last_two[1]['shift']*100:.3f}%  (at rs_min = {last_two[1]['rs_min']:.0f})."
            )
        elif delta < 0.5:
            verdict = (
                f"MARGINAL.  Shift between two deepest grids = {delta:.3f} pp. "
                f"Best estimate {last_two[1]['shift']*100:.3f}%, but deeper grids "
                f"would tighten the number."
            )
        else:
            verdict = (
                f"NOT YET CONVERGED.  Shift still moving by {delta:.3f} pp "
                f"between rs_min = {last_two[0]['rs_min']:.0f} and "
                f"rs_min = {last_two[1]['rs_min']:.0f}.  Need deeper throats."
            )
    else:
        verdict = "Insufficient converged rows to assess convergence."
    print()
    print(f"  VERDICT: {verdict}")
    print()

    # ----- Plot -----
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    rs_mins = [r["rs_min"] for r in converged_rows]
    shifts = [r["shift"] * 100 for r in converged_rows]
    eps_mins = [r["eps_min"] for r in converged_rows]
    calib_errs = [r["calib_err"] * 100 for r in converged_rows]

    ax = axes[0]
    ax.plot(rs_mins, shifts, "o-", linewidth=2, markersize=8, color="tab:green",
            label="V_exact analytic")
    ax.axhline(5.384, color="tab:orange", linestyle="--", alpha=0.7,
               label="G94 sweep mean (5.384%)")
    ax.fill_between(rs_mins, 5.384 - 0.118, 5.384 + 0.118,
                    color="tab:orange", alpha=0.15, label="G94 +/- 1 sigma")
    ax.set_xlabel("rs_min  (throat depth)")
    ax.set_ylabel("QNM shift vs GR  (%)")
    ax.set_title("Convergence of axial QNM shift with throat depth")
    ax.invert_xaxis()
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.semilogx(eps_mins, shifts, "o-", linewidth=2, markersize=8, color="tab:purple")
    ax.set_xlabel("eps_min = (r_min - 2M)/(2M)")
    ax.set_ylabel("QNM shift vs GR  (%)")
    ax.set_title("Shift vs minimum eps reached")
    ax.invert_xaxis()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = PLOTS / "G109_axial_convergence.png"
    plt.savefig(out_png, dpi=180)
    plt.close()
    print(f"Plot saved: {out_png}")

    # ----- Markdown -----
    md = ["# G109 - Axial QNM convergence study\n"]
    md.append("\nConvergence sweep of the analytic V_exact axial QNM with deepening throat.\n")
    md.append("\n## Sweep\n")
    md.append("| rs_min | N | eps_min | f_max | calib | omega_exact | shift |\n")
    md.append("|---:|---:|---:|---:|---:|---|---:|\n")
    for r in results:
        if not r.get("converged", False):
            md.append(f"| {r.get('rs_min','?')} | ABORT: {r.get('reason')} | | | | | |\n")
            continue
        md.append(f"| {r['rs_min']:.0f} | {r['N']} | {r['eps_min']:.2e} | "
                  f"{r['f_max']:.2e} | {r['calib_err']*100:.3f}% | "
                  f"{r['omega_exact']} | {r['shift']*100:.4f}% |\n")
    md.append("\n## Verdict\n")
    md.append(verdict + "\n")
    md.append("\n## Reference values\n")
    md.append("- G92 / G93 reported (G92 finite-difference): 5.35%\n")
    md.append("- G94 robustness sweep mean (FD): 5.384% +/- 0.118%\n")
    md.append("- G108 analytic (single grid, rs_min = -300): 4.977%\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G109_axial_convergence_study.py](../scripts/G109_axial_convergence_study.py)\n")
    md.append("- [plots/G109_axial_convergence.png](../plots/G109_axial_convergence.png)\n")
    out_md = RESULTS / "G109_axial_convergence_summary.md"
    out_md.write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {out_md}")


if __name__ == "__main__":
    main()
