#!/usr/bin/env python3
"""
G107_identify_Lmu_Dmu_decoherence.py

STAM / Model-A quantum-side worked example:
identify concrete L_mu and D_mu for environmental decoherence.

Purpose
=======
The framework has the structural resolution-rate form:

    Gamma_res = sum_mu <L_mu^dagger L_mu> D_mu

where:
    L_mu = physical interaction / write channel
    D_mu = distinguishability created by that channel

This script makes that concrete for a spatial superposition interacting with
an environment through scattering.

Core STAM interpretation
========================
Standard open quantum systems:
    environmental scattering suppresses off-diagonal coherence.

STAM interpretation:
    each physical interaction channel that creates distinguishability is a
    resolution/write channel.
    Gamma_res is the local proper-time SU-write intensity.
    Each write resolves 1 SU = A0 = 1/(12*pi).

For one scattering channel with rate gamma and wavelength lambda:

    L_k = sqrt(gamma_k) exp(i k x_hat)

    D(k, Delta_x) = 1 - sinc(k Delta_x)

Then:

    Gamma_res = gamma * D

    P_resolved(t) = 1 - exp(-Gamma_res t)

    N_write(t) = Gamma_res t

    A_ledger(t) = A0 * N_write(t)

For many channels:

    Gamma_res = sum_i gamma_i D_i

Cases implemented
=================
1. Single photon-scattering channel.
2. Thermal-ish photon bath integrated over a wavelength distribution.
3. Gas-collision decoherence using momentum transfer q.
4. Detector pixel absorption as high-distinguishability projectors.

Outputs
=======
results/G107_Lmu_Dmu_decoherence_summary.md
results/G107_Lmu_Dmu_decoherence_grid.csv
plots/G107_distinguishability_vs_separation.png
plots/G107_resolution_probability.png
plots/G107_write_count_vs_time.png
plots/G107_gamma_vs_wavelength.png
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


A0 = 1.0 / (12.0 * math.pi)


def sinc(x):
    """sin(x)/x, vectorized, with sinc(0)=1."""
    x = np.asarray(x, dtype=float)
    return np.where(np.abs(x) < 1e-12, 1.0, np.sin(x)/x)


def D_photon(delta_x_m, wavelength_m):
    """Approximate spatial distinguishability for isotropic photon scattering.

    k = 2*pi/lambda.
    D = 1 - sinc(k * Delta_x).
    """
    k = 2.0 * math.pi / wavelength_m
    return np.clip(1.0 - sinc(k * delta_x_m), 0.0, 2.0)


def D_gas(delta_x_m, q_kg_m_s):
    """Gas collision distinguishability from momentum transfer q.

    D = 1 - cos(q * Delta_x / hbar), direction-averaged version replaced by
    a bounded scalar channel for this diagnostic.
    """
    hbar = 1.054_571_817e-34
    phase = q_kg_m_s * delta_x_m / hbar
    return np.clip(1.0 - np.cos(phase), 0.0, 2.0)


def gamma_res_single(gamma_s, D):
    return gamma_s * D


def P_resolved(gamma_res_s, t_s):
    return 1.0 - np.exp(-gamma_res_s * t_s)


def N_write(gamma_res_s, t_s):
    return gamma_res_s * t_s


def A_ledger(gamma_res_s, t_s):
    return A0 * gamma_res_s * t_s


def thermal_photon_distribution(wavelengths_m, peak_wavelength_m, gamma_total_s):
    """Simple log-normal wavelength-rate distribution normalized to gamma_total."""
    lnlam = np.log(wavelengths_m)
    mu = np.log(peak_wavelength_m)
    sigma = 0.55
    weights = np.exp(-0.5*((lnlam-mu)/sigma)**2)
    weights = weights / np.trapz(weights, wavelengths_m)
    # convert density to channel rates on discrete grid
    rates_density = gamma_total_s * weights
    return rates_density


def integrate_gamma_res_photon_bath(delta_x_m, wavelengths_m, rates_density):
    Ds = D_photon(delta_x_m, wavelengths_m)
    return float(np.trapz(rates_density * Ds, wavelengths_m))


def main():
    p = argparse.ArgumentParser(description="G107 identify L_mu and D_mu decoherence worked example")
    p.add_argument("--outdir", default="/mnt/data/G107_decoherence_run")
    p.add_argument("--delta-x", type=float, default=1e-6, help="Superposition separation in meters")
    p.add_argument("--lambda", dest="lambda_m", type=float, default=500e-9, help="Single photon wavelength in meters")
    p.add_argument("--gamma", type=float, default=1e3, help="Single channel scattering rate in s^-1")
    p.add_argument("--time", type=float, default=0.01, help="Exposure time in seconds")
    args = p.parse_args()

    root = Path(args.outdir)
    results = root / "results"
    plots = root / "plots"
    results.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)

    dx = args.delta_x
    lam = args.lambda_m
    gamma = args.gamma
    t = args.time

    # Single photon channel
    D1 = float(D_photon(dx, lam))
    Gamma1 = gamma_res_single(gamma, D1)
    P1 = P_resolved(Gamma1, t)
    N1 = N_write(Gamma1, t)
    A1 = A_ledger(Gamma1, t)

    # Thermal-ish photon bath
    wavelengths = np.geomspace(100e-9, 100e-6, 600)
    rates_density = thermal_photon_distribution(wavelengths, peak_wavelength_m=10e-6, gamma_total_s=gamma)
    Gamma_bath = integrate_gamma_res_photon_bath(dx, wavelengths, rates_density)
    P_bath = P_resolved(Gamma_bath, t)
    N_bath = N_write(Gamma_bath, t)
    A_bath = A_ledger(Gamma_bath, t)

    # Gas collision channel
    # Example q from helium molecule thermal-ish velocity: q ~ m v with m ~ 6.6e-27 kg, v~1000 m/s.
    q = 6.6e-24
    Dgas = float(D_gas(dx, q))
    gamma_gas = gamma
    Gamma_gas = gamma_res_single(gamma_gas, Dgas)
    P_gas = P_resolved(Gamma_gas, t)
    N_gas = N_write(Gamma_gas, t)
    A_gas = A_ledger(Gamma_gas, t)

    # Detector pixel absorption
    # Projective channel: D≈1, rate gamma
    D_det = 1.0
    Gamma_det = gamma
    P_det = P_resolved(Gamma_det, t)
    N_det = N_write(Gamma_det, t)
    A_det = A_ledger(Gamma_det, t)

    rows = [
        {"case": "single_photon_scattering", "L_mu": "sqrt(gamma_k) exp(i k x_hat)", "D_mu": "1 - sinc(k Delta_x)", "gamma_s": gamma, "D": D1, "Gamma_res_s": Gamma1, "t_s": t, "P_resolved": P1, "N_write": N1, "A_ledger": A1},
        {"case": "thermal_photon_bath", "L_mu": "sqrt(gamma(lambda)) exp(i k(lambda) x_hat)", "D_mu": "integral[1 - sinc(k Delta_x)]", "gamma_s": gamma, "D": Gamma_bath/gamma if gamma else 0, "Gamma_res_s": Gamma_bath, "t_s": t, "P_resolved": P_bath, "N_write": N_bath, "A_ledger": A_bath},
        {"case": "gas_collision", "L_mu": "sqrt(Gamma_q) exp(i q x_hat / hbar)", "D_mu": "1 - cos(q Delta_x / hbar)", "gamma_s": gamma_gas, "D": Dgas, "Gamma_res_s": Gamma_gas, "t_s": t, "P_resolved": P_gas, "N_write": N_gas, "A_ledger": A_gas},
        {"case": "detector_absorption", "L_mu": "sqrt(gamma_i) |i><i|", "D_mu": "approximately 1", "gamma_s": gamma, "D": D_det, "Gamma_res_s": Gamma_det, "t_s": t, "P_resolved": P_det, "N_write": N_det, "A_ledger": A_det},
    ]

    csv_path = results / "G107_Lmu_Dmu_decoherence_grid.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        fields = list(rows[0].keys())
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # Plots
    dxs = np.geomspace(1e-10, 1e-3, 600)
    plt.figure(figsize=(10, 6))
    for lam_i in [100e-9, 500e-9, 1e-6, 10e-6, 100e-6]:
        plt.semilogx(dxs, D_photon(dxs, lam_i), label=f"lambda={lam_i:g} m")
    plt.axvline(dx, color="black", linestyle=":", label="chosen Δx")
    plt.xlabel("superposition separation Δx [m]")
    plt.ylabel("distinguishability D")
    plt.title("G107 distinguishability for photon scattering channels")
    plt.grid(alpha=0.3)
    plt.legend(fontsize=8)
    p1 = plots / "G107_distinguishability_vs_separation.png"
    plt.tight_layout()
    plt.savefig(p1, dpi=200)
    plt.close()

    ts = np.geomspace(1e-9, 10, 600)
    plt.figure(figsize=(10, 6))
    for r in rows:
        plt.semilogx(ts, P_resolved(r["Gamma_res_s"], ts), label=r["case"])
    plt.axvline(t, color="black", linestyle=":", label="chosen t")
    plt.xlabel("time [s]")
    plt.ylabel("P_resolved = 1 - exp(-Gamma_res t)")
    plt.title("G107 resolution probability")
    plt.grid(alpha=0.3)
    plt.legend(fontsize=8)
    p2 = plots / "G107_resolution_probability.png"
    plt.tight_layout()
    plt.savefig(p2, dpi=200)
    plt.close()

    plt.figure(figsize=(10, 6))
    for r in rows:
        plt.loglog(ts, np.maximum(N_write(r["Gamma_res_s"], ts), 1e-300), label=r["case"])
    plt.axvline(t, color="black", linestyle=":", label="chosen t")
    plt.xlabel("time [s]")
    plt.ylabel("N_write = Gamma_res t")
    plt.title("G107 expected SU write count")
    plt.grid(alpha=0.3)
    plt.legend(fontsize=8)
    p3 = plots / "G107_write_count_vs_time.png"
    plt.tight_layout()
    plt.savefig(p3, dpi=200)
    plt.close()

    plt.figure(figsize=(10, 6))
    gamma_res_lambda = rates_density * D_photon(dx, wavelengths)
    plt.loglog(wavelengths, rates_density, label="rate density gamma(lambda)")
    plt.loglog(wavelengths, gamma_res_lambda, label="resolution contribution gamma(lambda) D(lambda)")
    plt.axvline(lam, color="black", linestyle=":", label="single-channel lambda")
    plt.xlabel("wavelength [m]")
    plt.ylabel("rate density [s^-1 m^-1]")
    plt.title("G107 photon-bath channel contributions")
    plt.grid(alpha=0.3)
    plt.legend(fontsize=8)
    p4 = plots / "G107_gamma_vs_wavelength.png"
    plt.tight_layout()
    plt.savefig(p4, dpi=200)
    plt.close()

    # Summary
    summary = results / "G107_Lmu_Dmu_decoherence_summary.md"
    md = []
    md.append("# G107 — Identify L_mu and D_mu for decoherence channels\n\n")
    md.append("## Structural formula\n\n")
    md.append("```text\nGamma_res = sum_mu <L_mu^dagger L_mu> D_mu\n```\n\n")
    md.append("STAM interpretation:\n\n")
    md.append("- `L_mu` is a physical interaction / write channel.\n")
    md.append("- `D_mu` is the distinguishability created by that channel.\n")
    md.append("- `Gamma_res` is the proper-time SU-write intensity.\n")
    md.append("- each write resolves `1 SU = A0 = 1/(12π)`.\n\n")
    md.append(f"`A0 = {A0:.12f}`\n\n")
    md.append("## Chosen example\n\n")
    md.append(f"- separation Δx = `{dx:g}` m\n")
    md.append(f"- single photon wavelength λ = `{lam:g}` m\n")
    md.append(f"- single-channel rate γ = `{gamma:g}` s^-1\n")
    md.append(f"- exposure time t = `{t:g}` s\n\n")
    md.append("## Results\n\n")
    md.append("| case | L_mu | D_mu | gamma | D | Gamma_res | P_resolved | N_write | A_ledger |\n")
    md.append("|---|---|---|---:|---:|---:|---:|---:|---:|\n")
    for r in rows:
        md.append(f"| {r['case']} | `{r['L_mu']}` | `{r['D_mu']}` | {r['gamma_s']:.4e} | {r['D']:.4e} | {r['Gamma_res_s']:.4e} | {r['P_resolved']:.4e} | {r['N_write']:.4e} | {r['A_ledger']:.4e} |\n")
    md.append("\n## Interpretation\n\n")
    md.append("This is the first concrete identification of `L_mu` and `D_mu` for ordinary environmental decoherence. Photon scattering uses `L_k = sqrt(gamma_k) exp(i k x_hat)` and `D_k = 1 - sinc(k Delta_x)`. Gas collision uses `L_q = sqrt(Gamma_q) exp(i q x_hat / hbar)` and a momentum-transfer distinguishability. Detector absorption is the high-D limit with macroscopically distinguishable records.\n\n")
    md.append("In STAM language, standard decoherence rates become SU-write rates: `dN_write = Gamma_res dt`, `dA_ledger = A0 Gamma_res dt`.\n\n")
    md.append("## Files\n\n")
    md.append(f"- CSV: `{csv_path}`\n")
    md.append(f"- Plot: `{p1}`\n")
    md.append(f"- Plot: `{p2}`\n")
    md.append(f"- Plot: `{p3}`\n")
    md.append(f"- Plot: `{p4}`\n")
    summary.write_text("".join(md), encoding="utf-8")

    print("=" * 100)
    print("G107: identify L_mu and D_mu for decoherence")
    print("=" * 100)
    print(f"A0 = {A0:.12f}")
    print(f"Δx = {dx:g} m, lambda = {lam:g} m, gamma = {gamma:g} s^-1, t = {t:g} s")
    print()
    print(f"{'case':<26}{'D':>12}{'Gamma_res':>14}{'P_resolved':>14}{'N_write':>14}{'A_ledger':>14}")
    print("-" * 96)
    for r in rows:
        print(f"{r['case']:<26}{r['D']:>12.4e}{r['Gamma_res_s']:>14.4e}{r['P_resolved']:>14.4e}{r['N_write']:>14.4e}{r['A_ledger']:>14.4e}")
    print()
    print(f"Summary written: {summary}")
    print(f"CSV written: {csv_path}")
    print(f"Plots written: {plots}")


if __name__ == "__main__":
    main()
