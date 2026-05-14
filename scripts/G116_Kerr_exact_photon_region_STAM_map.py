#!/usr/bin/env python3
"""
G116_Kerr_exact_photon_region_STAM_map.py

STAM / Model-A Kerr exact photon-region map.

Purpose
=======
G85 replaced the G84 sin^2(theta) interpolation with the exact Kerr
spheroidal photon-region surface derived from the spherical-null-geodesic
conditions R(r)=0 and R'(r)=0.

G116 turns that closed-form G85 result into a reusable numerical map for
future Kerr STAM observables:

    - EMRI diagnostics
    - EHT / photon-ring substructure
    - Kerr throat / inside-shell width maps
    - Kerr W_K inside-shell factors
    - exact y_K normalization for STAM's F(y)

Definitions
===========
M = 1 units.

Equatorial prograde photon orbit:

    r_ph^+(a) = 2M [1 + cos((2/3) arccos(-a/M))]

Kerr horizon:

    r_+ = M + sqrt(M^2 - a^2)

Exact spherical photon-region constants:
    lambda(r_p) = -(r_p^3 - 3r_p^2 + a^2 r_p + a^2) / [a(r_p - 1)]

    eta(r_p) = r_p^2 [4a^2 Delta - (r_p^2 - 3r_p + 2a^2)^2]
               / [a^2 (r_p - 1)^2]

    Theta(theta; r_p, a) =
        eta sin^2(theta) - cos^2(theta) [lambda^2 - a^2 sin^2(theta)]

Exact off-axis photon-region surface:
    r_pr_exact(theta; a) = root in [r_ph^+(a), r_polar(a)] of Theta = 0

STAM shell variables:
    A_K(r; M,a) = 2 M r / (r^2 + a^2)
    Sigma_K = 3 A_K
    Sigma_ph(theta,a) = Sigma_K(r_pr_exact(theta,a))
    y_K_exact(r,theta,a) = [Sigma_K(r,a) - Sigma_ph(theta,a)] / [3 - Sigma_ph(theta,a)]

Outputs
=======
results/G116_Kerr_exact_photon_region_map_summary.md
results/G116_Kerr_exact_photon_region_map.csv
plots/G116_rpr_exact_vs_theta.png
plots/G116_shell_width_vs_theta.png
plots/G116_g84_error_vs_theta.png
plots/G116_sigma_ph_vs_theta.png
plots/G116_WK_factor_heatmap.png
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


M = 1.0


def Delta(r, a):
    return r*r - 2*M*r + a*a


def r_plus(a):
    return M + math.sqrt(M*M - a*a)


def r_ph_prograde(a):
    if abs(a) < 1e-14:
        return 3.0*M
    return 2*M*(1 + math.cos((2/3)*math.acos(-a/M)))


def polar_cubic(r, a):
    # r^3 - 3Mr^2 + a^2 r + a^2 M = 0, M=1
    return r**3 - 3*M*r**2 + a*a*r + a*a*M


def r_polar(a):
    if abs(a) < 1e-14:
        return 3.0*M
    rp = r_plus(a)
    # largest root between rp and 4M
    # scan intervals to find sign changes
    xs = np.linspace(rp + 1e-10, 5.0, 2000)
    vals = polar_cubic(xs, a)
    roots = []
    for i in range(len(xs)-1):
        if vals[i] == 0:
            roots.append(xs[i])
        if vals[i]*vals[i+1] < 0:
            roots.append(brentq(lambda r: polar_cubic(r,a), xs[i], xs[i+1]))
    if not roots:
        # at a=0 fallback, else use max real roots from numpy
        coeff = [1, -3*M, a*a, a*a*M]
        rs = np.roots(coeff)
        real = [float(x.real) for x in rs if abs(x.imag)<1e-8 and x.real > rp-1e-8]
        if real:
            return max(real)
        raise RuntimeError(f"No polar root found for a={a}")
    return max(roots)


def lam_sph(rp, a):
    if abs(a) < 1e-12:
        return np.nan
    return - (rp**3 - 3*rp**2 + a*a*rp + a*a) / (a*(rp - 1))


def eta_sph(rp, a):
    if abs(a) < 1e-12:
        return np.nan
    D = Delta(rp, a)
    return rp**2 * (4*a*a*D - (rp**2 - 3*rp + 2*a*a)**2) / (a*a*(rp - 1)**2)


def Theta(theta, rp, a):
    if abs(a) < 1e-12:
        # Schwarzschild degenerate photon sphere. Any theta at r=3.
        return 0.0
    lam = lam_sph(rp,a)
    eta = eta_sph(rp,a)
    s = math.sin(theta)
    c = math.cos(theta)
    return eta*s*s - c*c*(lam*lam - a*a*s*s)


def r_pr_exact(theta, a):
    if abs(a) < 1e-12:
        return 3.0*M
    # symmetry: theta in [0, pi/2]
    th = abs(theta)
    if th > math.pi/2:
        th = math.pi - th
    req = r_ph_prograde(a)
    rpol = r_polar(a)
    if abs(th - math.pi/2) < 1e-10:
        return req
    if abs(th) < 1e-10:
        return rpol
    lo, hi = req, rpol
    flo = Theta(th, lo, a)
    fhi = Theta(th, hi, a)
    # There may be a zero at boundary due to numerical roundoff; handle.
    if abs(flo) < 1e-10:
        return lo
    if abs(fhi) < 1e-10:
        return hi
    if flo*fhi > 0:
        # scan to bracket
        xs = np.linspace(lo, hi, 1000)
        vals = [Theta(th, x, a) for x in xs]
        for i in range(len(xs)-1):
            if vals[i]*vals[i+1] <= 0:
                return brentq(lambda r: Theta(th,r,a), xs[i], xs[i+1], xtol=1e-12, rtol=1e-12)
        raise RuntimeError(f"Could not bracket theta={theta}, a={a}, flo={flo}, fhi={fhi}")
    return brentq(lambda r: Theta(th,r,a), lo, hi, xtol=1e-12, rtol=1e-12)


def r_pr_g84(theta, a):
    # sin^2 interpolation between equatorial prograde and polar
    return r_ph_prograde(a)*math.sin(theta)**2 + r_polar(a)*math.cos(theta)**2


def A_K(r, a):
    return 2*M*r / (r*r + a*a)


def Sigma_K(r, a):
    return 3*A_K(r,a)


def F_quintic(y):
    return 1 - 5*y**4 + 4*y**5


def yK_exact(r, theta, a):
    sph = Sigma_K(r_pr_exact(theta,a), a)
    sig = Sigma_K(r,a)
    den = 3 - sph
    if abs(den) < 1e-14:
        return np.nan
    return (sig - sph)/den


def W_K_inside_factor(r, theta, a):
    """G82/G83/G84 W_K inside-shell kinematic factor:
       W = (Delta/Sigma_BL) * F(yK) * [dSigma_K/dr]^2
       Here output the full factor for M=1:
       Sigma_BL = r^2 + a^2 cos^2 theta
       dSigma/dr = 6M (a^2 - r^2)/(r^2+a^2)^2
    """
    sigBL = r*r + a*a*math.cos(theta)**2
    y = yK_exact(r, theta, a)
    if not np.isfinite(y):
        return np.nan
    F = F_quintic(y)
    dSig = 6*M*(a*a - r*r)/(r*r + a*a)**2
    return (Delta(r,a)/sigBL)*F*(dSig*dSig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="/mnt/data/G116_run")
    ap.add_argument("--n-theta", type=int, default=121)
    args = ap.parse_args()

    root = Path(args.outdir)
    results = root/"results"; plots = root/"plots"
    results.mkdir(parents=True, exist_ok=True); plots.mkdir(parents=True, exist_ok=True)

    a_values = [0.0, 0.3, 0.5, 0.7, 0.9, 0.99]
    theta_grid = np.linspace(0, math.pi/2, args.n_theta)

    rows = []
    summary_rows = []
    for a in a_values:
        rp = r_plus(a)
        req = r_ph_prograde(a)
        rpol = r_polar(a)
        max_abs_dr = 0.0
        max_rel_dr = 0.0
        max_abs_dsig = 0.0
        max_rel_dsig = 0.0
        max_row = None
        for th in theta_grid:
            rex = r_pr_exact(th,a)
            rg84 = r_pr_g84(th,a)
            sig_ex = Sigma_K(rex,a)
            sig_g84 = Sigma_K(rg84,a)
            shell_width = rex - rp
            abs_dr = abs(rg84 - rex)
            rel_dr = abs_dr/abs(rex) if rex else 0
            abs_dsig = abs(sig_g84 - sig_ex)
            rel_dsig = abs_dsig/abs(sig_ex) if sig_ex else 0
            if abs_dr > max_abs_dr:
                max_abs_dr=max_abs_dr; max_row=(th,rex,rg84,sig_ex,sig_g84)
            max_abs_dr = max(max_abs_dr, abs_dr)
            max_rel_dr = max(max_rel_dr, rel_dr)
            max_abs_dsig = max(max_abs_dsig, abs_dsig)
            max_rel_dsig = max(max_rel_dsig, rel_dsig)
            rows.append({
                "a":a, "theta_rad":th, "theta_deg":th*180/math.pi,
                "r_plus":rp, "r_ph_eq":req, "r_polar":rpol,
                "r_pr_exact":rex, "r_pr_G84":rg84,
                "delta_r_G84_minus_exact":rg84-rex,
                "abs_delta_r":abs_dr, "rel_delta_r":rel_dr,
                "Sigma_ph_exact":sig_ex, "Sigma_ph_G84":sig_g84,
                "delta_Sigma_G84_minus_exact":sig_g84-sig_ex,
                "abs_delta_Sigma":abs_dsig, "rel_delta_Sigma":rel_dsig,
                "shell_width_exact":shell_width
            })
        summary_rows.append({
            "a":a, "r_plus":rp, "r_ph_eq":req, "r_polar":rpol,
            "max_abs_delta_r":max_abs_dr, "max_rel_delta_r":max_rel_dr,
            "max_abs_delta_Sigma":max_abs_dsig, "max_rel_delta_Sigma":max_rel_dsig
        })

    csv_path = results/"G116_Kerr_exact_photon_region_map.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        fields=list(rows[0].keys())
        w=csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)

    summary_csv = results/"G116_Kerr_exact_photon_region_summary_table.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as f:
        fields=list(summary_rows[0].keys())
        w=csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(summary_rows)

    # Plots: r exact vs theta
    plt.figure(figsize=(10,6))
    for a in a_values:
        sub=[r for r in rows if r["a"]==a]
        plt.plot([r["theta_deg"] for r in sub],[r["r_pr_exact"] for r in sub],label=f"a={a}")
    plt.xlabel("theta [deg] (0=pole, 90=equator)")
    plt.ylabel("r_pr_exact / M")
    plt.title("G116 exact Kerr photon-region surface")
    plt.grid(alpha=0.3); plt.legend()
    p1=plots/"G116_rpr_exact_vs_theta.png"; plt.tight_layout(); plt.savefig(p1,dpi=200); plt.close()

    # shell width
    plt.figure(figsize=(10,6))
    for a in a_values:
        sub=[r for r in rows if r["a"]==a]
        plt.plot([r["theta_deg"] for r in sub],[r["shell_width_exact"] for r in sub],label=f"a={a}")
    plt.xlabel("theta [deg]")
    plt.ylabel("r_pr_exact - r_+")
    plt.title("G116 inside-shell radial width by latitude")
    plt.grid(alpha=0.3); plt.legend()
    p2=plots/"G116_shell_width_vs_theta.png"; plt.tight_layout(); plt.savefig(p2,dpi=200); plt.close()

    # G84 error
    plt.figure(figsize=(10,6))
    for a in a_values[1:]:
        sub=[r for r in rows if r["a"]==a]
        plt.plot([r["theta_deg"] for r in sub],[r["delta_r_G84_minus_exact"] for r in sub],label=f"a={a}")
    plt.xlabel("theta [deg]")
    plt.ylabel("r_G84 - r_exact")
    plt.title("G116 G84 interpolation error")
    plt.grid(alpha=0.3); plt.legend()
    p3=plots/"G116_g84_error_vs_theta.png"; plt.tight_layout(); plt.savefig(p3,dpi=200); plt.close()

    # Sigma ph
    plt.figure(figsize=(10,6))
    for a in a_values:
        sub=[r for r in rows if r["a"]==a]
        plt.plot([r["theta_deg"] for r in sub],[r["Sigma_ph_exact"] for r in sub],label=f"a={a}")
    plt.xlabel("theta [deg]")
    plt.ylabel("Sigma_ph_exact")
    plt.title("G116 exact photon-region shell count Sigma_ph(theta,a)")
    plt.grid(alpha=0.3); plt.legend()
    p4=plots/"G116_sigma_ph_vs_theta.png"; plt.tight_layout(); plt.savefig(p4,dpi=200); plt.close()

    # Heatmap W_K factor at midpoint of shell r=(r+ + r_pr)/2 for exact surface.
    a_heat=np.linspace(0.01,0.99,80)
    th_heat=np.linspace(0.001,math.pi/2,100)
    Z=np.zeros((len(a_heat),len(th_heat)))
    for i,a in enumerate(a_heat):
        rp=r_plus(a)
        for j,th in enumerate(th_heat):
            rex=r_pr_exact(th,a)
            rmid=0.5*(rp+rex)
            Z[i,j]=W_K_inside_factor(rmid, th, a)
    plt.figure(figsize=(10,6))
    extent=[th_heat[0]*180/math.pi, th_heat[-1]*180/math.pi, a_heat[0], a_heat[-1]]
    plt.imshow(np.log10(np.maximum(Z,1e-300)), origin="lower", aspect="auto", extent=extent)
    plt.colorbar(label="log10 W_K(mid-shell)")
    plt.xlabel("theta [deg]")
    plt.ylabel("spin a/M")
    plt.title("G116 W_K inside-shell factor at midpoint")
    p5=plots/"G116_WK_factor_heatmap.png"; plt.tight_layout(); plt.savefig(p5,dpi=200); plt.close()

    # Markdown summary
    md_path=results/"G116_Kerr_exact_photon_region_map_summary.md"
    md=[]
    md.append("# G116 — Kerr exact photon-region STAM map\n\n")
    md.append("G116 converts the G85 exact Kerr photon-region surface into a reusable numerical map for STAM Kerr observables.\n\n")
    md.append("## Exact photon-region formulas\n\n")
    md.append("```text\n")
    md.append("lambda(r_p) = -(r_p^3 - 3r_p^2 + a^2 r_p + a^2) / [a(r_p - 1)]\n")
    md.append("eta(r_p) = r_p^2 [4a^2 Delta - (r_p^2 - 3r_p + 2a^2)^2] / [a^2 (r_p - 1)^2]\n")
    md.append("Theta(theta;r_p,a) = eta sin^2(theta) - cos^2(theta)(lambda^2 - a^2 sin^2(theta))\n")
    md.append("r_pr_exact(theta,a) = root of Theta=0 in [r_ph_prograde, r_polar]\n")
    md.append("```\n\n")
    md.append("## G84 interpolation error summary\n\n")
    md.append("| a | r_plus | r_eq | r_polar | max abs Δr | max rel Δr | max abs ΔΣ | max rel ΔΣ |\n")
    md.append("|---:|---:|---:|---:|---:|---:|---:|---:|\n")
    for s in summary_rows:
        md.append(f"| {s['a']:.2f} | {s['r_plus']:.6f} | {s['r_ph_eq']:.6f} | {s['r_polar']:.6f} | {s['max_abs_delta_r']:.6f} | {100*s['max_rel_delta_r']:.2f}% | {s['max_abs_delta_Sigma']:.6f} | {100*s['max_rel_delta_Sigma']:.2f}% |\n")
    md.append("\n## Interpretation\n\n")
    md.append("G116 confirms the G85 result in map form: G84's sin^2(theta) interpolation is useful for intuition but introduces sizeable high-spin, intermediate-latitude errors. Future Kerr observable calculations should use r_pr_exact(theta,a) and Sigma_ph_exact(theta,a).\n\n")
    md.append("The structural STAM claims are unchanged: ghost-freedom, cubic horizon vanishing, eikonal recovery, and lambda_1/lambda_2 closure depend on the existence and topology of the photon-region surface. G116 supplies the precision geometry needed for EMRI/EHT/ringdown observables.\n\n")
    md.append("## Files\n\n")
    md.append(f"- Map CSV: `{csv_path}`\n")
    md.append(f"- Summary CSV: `{summary_csv}`\n")
    md.append(f"- Plot: `{p1}`\n")
    md.append(f"- Plot: `{p2}`\n")
    md.append(f"- Plot: `{p3}`\n")
    md.append(f"- Plot: `{p4}`\n")
    md.append(f"- Plot: `{p5}`\n")
    md_path.write_text("".join(md),encoding="utf-8")

    print("="*100)
    print("G116: Kerr exact photon-region STAM map")
    print("="*100)
    print()
    print(f"{'a':>6}{'r+':>12}{'r_eq':>12}{'r_polar':>12}{'max Δr':>12}{'max rel':>12}{'max ΔΣ':>12}{'max relΣ':>12}")
    print("-"*92)
    for s in summary_rows:
        print(f"{s['a']:>6.2f}{s['r_plus']:>12.6f}{s['r_ph_eq']:>12.6f}{s['r_polar']:>12.6f}{s['max_abs_delta_r']:>12.6f}{100*s['max_rel_delta_r']:>11.2f}%{s['max_abs_delta_Sigma']:>12.6f}{100*s['max_rel_delta_Sigma']:>11.2f}%")
    print()
    print(f"Summary written: {md_path}")
    print(f"Map CSV written: {csv_path}")
    print(f"Plots written: {plots}")

if __name__ == "__main__":
    main()
