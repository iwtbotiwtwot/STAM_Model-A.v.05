#!/usr/bin/env python3
"""
G103_greybody_scattering_STAM_throat.py

STAM / Model-A greybody / scattering diagnostic from the near-horizon throat.

Purpose
=======
G90 showed that the STAM final-shell metric changes the near-horizon tortoise
distance:

    (dr*/dr)_STAM / (dr*/dr)_GR = 1/sqrt(F(y))

with a power-law throat near the horizon instead of GR's logarithmic throat.

G91-G94 then showed that the corrected STAM tortoise coordinate gives a
calibrated low-l time-domain QNM diagnostic, and the action-aware f(Σ)R
correction materially changes the axial potential.

G103 asks a cleaner scattering question:

    For real frequencies ω, how do axial-wave greybody transmission
    probabilities differ between:
      1. Schwarzschild GR
      2. STAM metric-only proxy
      3. STAM action-aware potential?

This avoids complex-frequency QNM boundary instabilities and directly probes
the same near-horizon wave-throat structure.

Method
======
Solve the 1D wave equation in tortoise coordinate x = r*:

    ψ''(x) + [ω² - V(x)] ψ(x) = 0

Boundary setup:
    left/horizon side: transmitted ingoing wave normalized to 1
        ψ ~ exp(-iωx)

    right/infinity side:
        ψ ~ A_in exp(-iωx) + A_out exp(+iωx)

Transmission probability:
    T = 1 / |A_in|²

Reflection probability:
    R = |A_out/A_in|²

For a clean real potential, T + R should be close to 1.  This is used as a
sanity check.

Potentials
==========
Schwarzschild:
    V_GR = h [l(l+1)/r² - 6M/r³]

STAM proxy:
    V_proxy = h [l(l+1)/r² - 6M_eff/r³]
    M_eff = (r/2)(1-k)
    k = h F(y) inside PS, k=h outside

STAM action-aware:
    V_action = V_proxy + (sqrt(f))'' / sqrt(f)
    derivatives taken with respect to STAM tortoise coordinate.

Caveat
======
This is still a diagnostic axial potential, not a full greybody theorem from
the final constrained-Σ perturbation derivation.  But it is a cleaner
observable-channel test than the early QNM extraction attempts because it uses
real-frequency scattering and calibrated tortoise coordinates.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.special import lambertw


M = 1.0


def h_fn(r):
    return 1.0 - 2.0*M/r


def F_quintic(y):
    return 1.0 - 5.0*y**4 + 4.0*y**5


def y_of_r(r):
    return 6.0*M/r - 2.0


def k_GR(r):
    return h_fn(r)


def k_STAM(r):
    h = h_fn(r)
    y = y_of_r(r)
    return np.where(y > 0, h*F_quintic(y), h)


def V_RW_proxy(r, ell, k_func):
    h = h_fn(r)
    k = k_func(r)
    M_eff = 0.5*r*(1.0-k)
    return h*(ell*(ell+1)/r**2 - 6.0*M_eff/r**3)


def r_of_x_GR(x):
    """Invert Schwarzschild tortoise coordinate anchored at r*(3M)=0."""
    offset = 3.0*M + 2.0*M*np.log(0.5)
    R = x + offset
    arg = np.exp(R/(2.0*M) - 1.0)
    u = np.real(lambertw(arg))
    return 2.0*M*(1.0 + u)


def build_r_of_x_STAM(xgrid):
    """Build r(x) for STAM by integrating dr/dx = sqrt(h k), anchored r(0)=3M."""
    xgrid = np.asarray(xgrid)
    rarr = np.zeros_like(xgrid, dtype=float)

    def rhs(x, y):
        r = y[0]
        return [float(np.sqrt(max(h_fn(r)*k_STAM(np.array([r]))[0], 0.0)))]

    # positive side: same as GR outside PS, but integrate for consistency
    pos_mask = xgrid >= 0
    if np.any(pos_mask):
        xpos = np.sort(xgrid[pos_mask])
        sol = solve_ivp(rhs, (0.0, float(xpos[-1])), [3.0*M],
                        t_eval=xpos, rtol=1e-10, atol=1e-12, max_step=1.0)
        vals = sol.y[0]
        # map sorted values back
        mapping = dict(zip(np.round(xpos, 12), vals))
        rarr[pos_mask] = [mapping[np.round(x, 12)] for x in xgrid[pos_mask]]

    # negative side
    neg_mask = xgrid < 0
    if np.any(neg_mask):
        xneg = np.sort(xgrid[neg_mask])[::-1]  # from near 0 down
        sol = solve_ivp(rhs, (0.0, float(xneg[-1])), [3.0*M],
                        t_eval=xneg, rtol=1e-10, atol=1e-12, max_step=0.2)
        vals = sol.y[0]
        mapping = dict(zip(np.round(xneg, 12), vals))
        rarr[neg_mask] = [mapping[np.round(x, 12)] for x in xgrid[neg_mask]]

    return rarr


def dlnf_dA(A):
    y = 3*A - 2
    F = F_quintic(y)
    # Outside PS use f=1; inside formula.
    val = -2.0*y**3*(45*A**2 - 33*A - 13)/(A*F)
    return np.where(y > 0, val, 0.0)


def build_f_on_STAM_grid(xgrid, rgrid):
    """Integrate ln f in x, anchored f(0)=1.  Formula active only inside PS."""
    Agrid = 2*M/rgrid
    # dA/dr = -2M/r^2; dr/dx = sqrt(hk)
    drdx = np.sqrt(np.maximum(h_fn(rgrid)*k_STAM(rgrid), 0.0))
    dAdx = (-2*M/rgrid**2)*drdx
    dlnfdx_grid = dlnf_dA(Agrid)*dAdx

    # Since xgrid is uniform sorted, integrate cumulative trapezoid with anchor at x=0.
    ln_f = np.zeros_like(xgrid)
    # find index closest to zero
    i0 = int(np.argmin(np.abs(xgrid)))
    ln_f[i0] = 0.0
    # integrate to right
    for i in range(i0+1, len(xgrid)):
        dx = xgrid[i]-xgrid[i-1]
        ln_f[i] = ln_f[i-1] + 0.5*(dlnfdx_grid[i]+dlnfdx_grid[i-1])*dx
    # integrate to left
    for i in range(i0-1, -1, -1):
        dx = xgrid[i+1]-xgrid[i]
        ln_f[i] = ln_f[i+1] - 0.5*(dlnfdx_grid[i]+dlnfdx_grid[i+1])*dx
    # avoid overflow in deep throat; grids should be moderate.
    return np.exp(np.clip(ln_f, -100, 100))


def action_correction_from_f(xgrid, fgrid):
    sqrtf = np.sqrt(fgrid)
    d1 = np.gradient(sqrtf, xgrid, edge_order=2)
    d2 = np.gradient(d1, xgrid, edge_order=2)
    return d2/np.maximum(sqrtf, 1e-300)


def build_potentials(ell=2, x_min=-160.0, x_max=300.0, N=7000):
    x = np.linspace(x_min, x_max, N)

    # GR on GR tortoise
    r_gr = r_of_x_GR(x)
    V_gr = V_RW_proxy(r_gr, ell, k_GR)

    # STAM on STAM tortoise
    r_st = build_r_of_x_STAM(x)
    V_proxy = V_RW_proxy(r_st, ell, k_STAM)
    fgrid = build_f_on_STAM_grid(x, r_st)
    corr = action_correction_from_f(x, fgrid)
    V_action = V_proxy + corr

    # enforce exterior equality for action correction far outside? f=1 gives corr≈numerical noise.
    # Clean tiny numerical corr outside PS.
    V_action = np.where(r_st >= 3.0*M, V_proxy, V_action)

    return x, r_gr, V_gr, r_st, V_proxy, V_action, fgrid, corr


def scatter_for_omega(omega, xgrid, Vgrid):
    """Integrate from left with unit ingoing horizon amplitude and read A_in/A_out at right."""
    xL = float(xgrid[0])
    xR = float(xgrid[-1])
    V_interp = interp1d(xgrid, Vgrid, kind="linear", bounds_error=False,
                        fill_value=(Vgrid[0], Vgrid[-1]))

    psi0 = np.exp(-1j*omega*xL)
    psip0 = -1j*omega*psi0

    def rhs(x, y):
        psi = y[0] + 1j*y[1]
        psip = y[2] + 1j*y[3]
        V = float(V_interp(x))
        psipp = -(omega**2 - V)*psi
        return [psip.real, psip.imag, psipp.real, psipp.imag]

    sol = solve_ivp(rhs, (xL, xR), [psi0.real, psi0.imag, psip0.real, psip0.imag],
                    rtol=1e-8, atol=1e-10, max_step=0.25)

    if not sol.success:
        return np.nan, np.nan, np.nan, np.nan

    yend = sol.y[:, -1]
    psi = yend[0] + 1j*yend[1]
    psip = yend[2] + 1j*yend[3]

    Eminus = np.exp(-1j*omega*xR)
    Eplus = np.exp(1j*omega*xR)

    # solve [Eminus Eplus; -iw Eminus iw Eplus] [Ain,Aout] = [psi,psip]
    mat = np.array([[Eminus, Eplus],
                    [-1j*omega*Eminus, 1j*omega*Eplus]], dtype=complex)
    Ain, Aout = np.linalg.solve(mat, np.array([psi, psip], dtype=complex))

    T = 1.0/(abs(Ain)**2)
    R = abs(Aout/Ain)**2
    return float(T.real), float(R.real), float((T+R).real), abs(Ain)


def main():
    parser = argparse.ArgumentParser(description="G103 STAM greybody scattering diagnostic")
    parser.add_argument("--outdir", default="/mnt/data/G103_run")
    parser.add_argument("--ell", type=int, default=2)
    parser.add_argument("--x-min", type=float, default=-160.0)
    parser.add_argument("--x-max", type=float, default=300.0)
    parser.add_argument("--N", type=int, default=7000)
    args = parser.parse_args()

    root = Path(args.outdir)
    results = root/"results"; plots = root/"plots"
    results.mkdir(parents=True, exist_ok=True); plots.mkdir(parents=True, exist_ok=True)

    print("="*100)
    print("G103: greybody/scattering diagnostic from STAM throat")
    print("="*100)
    print(f"ell={args.ell}, x=[{args.x_min},{args.x_max}], N={args.N}")
    print("Building potentials...")
    x, r_gr, V_gr, r_st, V_proxy, V_action, fgrid, corr = build_potentials(args.ell, args.x_min, args.x_max, args.N)

    # Frequencies
    omegas = np.concatenate([
        np.linspace(0.05, 0.30, 11),
        np.linspace(0.35, 1.20, 18)
    ])

    rows = []
    for omega in omegas:
        Tg,Rg,Sg,Aing = scatter_for_omega(omega, x, V_gr)
        Tp,Rp,Sp,Ainp = scatter_for_omega(omega, x, V_proxy)
        Ta,Ra,Sa,Aina = scatter_for_omega(omega, x, V_action)
        rows.append({
            "omega": omega,
            "T_GR": Tg, "R_GR": Rg, "unitarity_GR": Sg,
            "T_proxy": Tp, "R_proxy": Rp, "unitarity_proxy": Sp,
            "T_action": Ta, "R_action": Ra, "unitarity_action": Sa,
            "T_proxy_over_GR": Tp/Tg if Tg>0 else np.nan,
            "T_action_over_GR": Ta/Tg if Tg>0 else np.nan,
            "T_action_over_proxy": Ta/Tp if Tp>0 else np.nan,
        })
        print(f"ω={omega:.3f}  T_GR={Tg:.4e}  T_proxy={Tp:.4e}  T_action={Ta:.4e}  "
              f"ratio_action/GR={Ta/Tg if Tg>0 else np.nan:.3f}  unitarity(action)={Sa:.5f}")

    # CSV
    csv_path = results/"G103_greybody_scattering.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(rows[0].keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)

    # Plots
    omega_arr = np.array([r["omega"] for r in rows])
    Tgr = np.array([r["T_GR"] for r in rows])
    Tpr = np.array([r["T_proxy"] for r in rows])
    Tac = np.array([r["T_action"] for r in rows])

    plt.figure(figsize=(10,6))
    plt.plot(omega_arr, Tgr, "o-", label="GR")
    plt.plot(omega_arr, Tpr, "o-", label="STAM proxy")
    plt.plot(omega_arr, Tac, "o-", label="STAM action-aware")
    plt.yscale("log")
    plt.xlabel("omega M")
    plt.ylabel("Transmission probability T_l(omega)")
    plt.title(f"G103 greybody transmission diagnostic, ell={args.ell}")
    plt.grid(alpha=0.3)
    plt.legend()
    p1 = plots/"G103_transmission_curves.png"
    plt.tight_layout(); plt.savefig(p1,dpi=200); plt.close()

    plt.figure(figsize=(10,6))
    plt.plot(omega_arr, Tpr/Tgr, "o-", label="proxy / GR")
    plt.plot(omega_arr, Tac/Tgr, "o-", label="action / GR")
    plt.axhline(1,color="black",lw=0.8)
    plt.xlabel("omega M")
    plt.ylabel("Transmission ratio")
    plt.title(f"G103 STAM/GR greybody transmission ratio, ell={args.ell}")
    plt.grid(alpha=0.3)
    plt.legend()
    p2 = plots/"G103_transmission_ratios.png"
    plt.tight_layout(); plt.savefig(p2,dpi=200); plt.close()

    plt.figure(figsize=(10,6))
    plt.plot(x, V_gr, label="GR")
    plt.plot(x, V_proxy, label="STAM proxy")
    plt.plot(x, V_action, label="STAM action-aware")
    plt.xlim(-60, 80)
    plt.xlabel("r*")
    plt.ylabel("V(r*)")
    plt.title("G103 potentials in tortoise coordinate")
    plt.grid(alpha=0.3)
    plt.legend()
    p3 = plots/"G103_potentials.png"
    plt.tight_layout(); plt.savefig(p3,dpi=200); plt.close()

    # Summary metrics around representative frequencies
    def nearest_row(w):
        return min(rows, key=lambda r: abs(r["omega"]-w))
    reps = [nearest_row(w) for w in [0.1,0.2,0.4,0.8,1.0]]

    summary = results/"G103_greybody_scattering_summary.md"
    md = []
    md.append("# G103 — Greybody / scattering diagnostic from STAM throat\n\n")
    md.append("## Purpose\n\n")
    md.append("G103 computes real-frequency axial-wave scattering for GR, STAM metric-only proxy, and STAM action-aware potentials. It probes the near-horizon tortoise-distance enhancement from G90 without relying on complex-frequency QNM extraction.\n\n")
    md.append("## Method\n\n")
    md.append("Solve:\n\n")
    md.append("```text\nψ''(r*) + [ω² - V(r*)]ψ = 0\n```\n\n")
    md.append("with a unit transmitted ingoing wave at the horizon and incident/reflected decomposition at infinity. Transmission is `T = 1/|A_in|²`.\n\n")
    md.append("## Representative results\n\n")
    md.append("| omega M | T_GR | T_proxy | T_action | proxy/GR | action/GR | T+R action |\n")
    md.append("|---:|---:|---:|---:|---:|---:|---:|\n")
    for r in reps:
        md.append(f"| {r['omega']:.3f} | {r['T_GR']:.4e} | {r['T_proxy']:.4e} | {r['T_action']:.4e} | "
                  f"{r['T_proxy_over_GR']:.3f} | {r['T_action_over_GR']:.3f} | {r['unitarity_action']:.5f} |\n")
    md.append("\n## Interpretation\n\n")
    md.append("If the action-aware transmission differs systematically from GR while unitarity remains close to one, this is a cleaner observable-channel signature of the STAM throat than the early low-l QNM attempts. The result should still be treated as a diagnostic until the full constrained-Σ axial perturbation equation is derived.\n\n")
    md.append("## Files\n\n")
    md.append(f"- CSV: `{csv_path}`\n")
    md.append(f"- Plot: `{p1}`\n")
    md.append(f"- Plot: `{p2}`\n")
    md.append(f"- Plot: `{p3}`\n")
    summary.write_text("".join(md), encoding="utf-8")

    print(f"CSV written: {csv_path}")
    print(f"Summary written: {summary}")
    print(f"Plot written: {p1}")
    print(f"Plot written: {p2}")
    print(f"Plot written: {p3}")


if __name__ == "__main__":
    main()
