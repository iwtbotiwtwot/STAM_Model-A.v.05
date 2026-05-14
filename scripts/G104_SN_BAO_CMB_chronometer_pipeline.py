#!/usr/bin/env python3
"""
G104_SN_BAO_CMB_chronometer_pipeline.py

Robust STAM / Model-A real cosmology pipeline scaffold:
SN + BAO + CMB compressed + chronometers.

Purpose
=======
This script is built to run in two modes:

1. TEMPLATE / SYNTHETIC MODE
   If data files are missing, it writes template CSVs and runs a synthetic
   sanity check. This verifies that the pipeline, model, likelihoods, and
   output machinery work.

2. REAL-DATA MODE
   If data files are present, it loads them and fits:
      - baseline flat LCDM
      - STAM distance-bias model

Core STAM convention from G98/G99
=================================
Photon-A traversal is a comoving/path correction first:

    D_C_obs(z) = D_C_intrinsic(z) + f_los(z) * b * shape(z)
    D_L_obs(z) = (1+z) * D_C_obs(z)

with:

    A0 = 1/(12π)
    b  = A0 * c / H0

Direct expansion probes are NOT photon-distance biased:

    H_model(z) = H0 * E(z)
    D_H(z) = c / H(z)

This is the key separation:
    SN, BAO transverse distances, CMB distances -> may receive photon-A bias
    BAO radial distance and chronometers -> direct expansion pressure tests

Important honesty
=================
This is a robust pipeline scaffold, not a final publication likelihood.
Full publication work still requires:
    - real covariance matrices for SN/BAO/CMB where available
    - nuisance parameters and priors appropriate to each dataset
    - proper treatment of r_d and CMB acoustic physics
    - MCMC / nested sampling for final uncertainties

CSV inputs
==========
Default data directory: ./data relative to run directory or --data-dir.

SN CSV: sn.csv
Required columns:
    z, mu
Optional:
    sigma_mu

BAO CSV: bao.csv
Required columns:
    z, observable, value, sigma
where observable is one of:
    DM_over_rd   transverse comoving distance / r_d
    DH_over_rd   radial Hubble distance / r_d = c/[H(z) r_d]
    DV_over_rd   volume distance / r_d
    DL_over_rd   luminosity distance / r_d   (rare; included for completeness)

CMB CSV: cmb.csv
Required columns:
    observable, value, sigma
where observable is one of:
    R_shift       sqrt(Ωm) H0 D_M(z*) / c
    lA            acoustic angular scale proxy π D_M(z*) / r_s
    DC_over_rd    D_C(z*) / r_d
Optional:
    zstar          default 1089

Chronometer CSV: chronometers.csv
Required columns:
    z, H, sigma_H

Models
======
baseline_lcdm:
    parameters: H0, Omega_m, optional SN offset Mcal, optional r_d

stam_bias:
    parameters: H0, Omega_m, f0, f1, optional SN offset Mcal, optional r_d
    f_los(z) = f0 + f1 z/(1+z)    (default)
    shape(z) = linear / stam_su / saturating / log

Outputs
=======
results/G104_fit_summary.md
results/G104_best_fit_params.csv
results/G104_probe_chi2_breakdown.csv
results/G104_model_comparison.csv
plots/G104_SN_hubble_residuals.png
plots/G104_BAO_residuals.png
plots/G104_chronometer_Hz.png
plots/G104_flos_shape.png
plots/G104_chi2_breakdown.png
data/template_*.csv if data missing

Examples
========
Synthetic sanity check:
    python G104_SN_BAO_CMB_chronometer_pipeline.py --outdir G104_run

Real data:
    python G104_SN_BAO_CMB_chronometer_pipeline.py --data-dir data --outdir G104_run --real-only

Fit with r_d free:
    python G104_SN_BAO_CMB_chronometer_pipeline.py --free-rd
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, List, Tuple, Any

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, minimize


# ---------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------

c_km_s = 299_792.458
A0 = 1.0 / (12.0 * math.pi)
MPC_TO_MLY = 3.261563776


# ---------------------------------------------------------------------
# Core cosmology
# ---------------------------------------------------------------------

def E_flat_lcdm(z, omega_m):
    z = np.asarray(z, dtype=float)
    omega_m = float(omega_m)
    return np.sqrt(omega_m * (1.0 + z)**3 + (1.0 - omega_m))


def H_z(z, H0, omega_m):
    return float(H0) * E_flat_lcdm(z, omega_m)


def comoving_distance_interp(z_values, H0, omega_m):
    """Return D_C(z) for arbitrary z array using a dense trapezoid grid.

    This avoids a slow per-point quad and is robust enough for fitting.
    """
    z_values = np.asarray(z_values, dtype=float)
    if z_values.size == 0:
        return np.array([])
    zmax = float(np.max(z_values))
    if zmax <= 0:
        return np.zeros_like(z_values)

    # Dense enough for z up to last scattering. Use log spacing plus linear
    # anchors to handle low-z and high-z smoothly.
    if zmax > 50:
        # include many points at low z and log high z
        z_low = np.linspace(0, min(5.0, zmax), 4000)
        z_high = np.geomspace(max(5.0, 1e-4), zmax, 4000)
        grid = np.unique(np.concatenate([z_low, z_high]))
    else:
        grid = np.linspace(0, zmax, max(3000, int(2000 * max(zmax, 1))))

    integrand = c_km_s / H_z(grid, H0, omega_m)
    dz = np.diff(grid)
    dc_grid = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * dz)])
    return np.interp(z_values, grid, dc_grid)


def D_L_from_D_C(z, D_C):
    return (1.0 + np.asarray(z, dtype=float)) * np.asarray(D_C, dtype=float)


def D_V(z, D_M, H):
    """Volume-averaged BAO distance D_V = [z D_M^2 c/H]^(1/3)."""
    z = np.asarray(z, dtype=float)
    return (z * D_M**2 * c_km_s / H)**(1.0 / 3.0)


def distance_modulus(D_L_mpc):
    D_L_mpc = np.asarray(D_L_mpc, dtype=float)
    return 5.0 * np.log10(D_L_mpc) + 25.0


def bridge_b_mpc(H0):
    return A0 * c_km_s / float(H0)


def bias_shape(z, model):
    z = np.asarray(z, dtype=float)
    if model == "linear":
        return z
    if model == "stam_su":
        return z * (1.0 + 3.0 * z / 20.0)
    if model == "saturating":
        return z / (1.0 + z)
    if model == "log":
        return np.log1p(z)
    raise ValueError(f"Unknown bias shape: {model}")


def flos_z(z, f0, f1, model="one_plus_z"):
    z = np.asarray(z, dtype=float)
    if model == "constant":
        return f0 + 0.0 * z
    if model == "one_plus_z":
        return f0 + f1 * z / (1.0 + z)
    if model == "log":
        return f0 + f1 * np.log1p(z)
    raise ValueError(f"Unknown f_los model: {model}")


def model_distances(z, params, model_kind, bias_shape_name, flos_model):
    """Return dictionary of distances for z.

    params keys:
        H0, omega_m, f0, f1, r_d, Mcal
    """
    H0 = params["H0"]
    om = params["omega_m"]
    z = np.asarray(z, dtype=float)
    D_C = comoving_distance_interp(z, H0, om)

    if model_kind == "stam":
        f0 = params.get("f0", 0.0)
        f1 = params.get("f1", 0.0)
        bias = bridge_b_mpc(H0) * flos_z(z, f0, f1, flos_model) * bias_shape(z, bias_shape_name)
        D_C_obs = D_C + bias
    else:
        D_C_obs = D_C

    H = H_z(z, H0, om)
    D_L = D_L_from_D_C(z, D_C_obs)
    D_V_val = D_V(z, D_C_obs, H)

    return {
        "D_C": D_C,
        "D_M_obs": D_C_obs,   # flat universe transverse comoving distance
        "D_L_obs": D_L,
        "D_V_obs": D_V_val,
        "H": H,
        "D_H": c_km_s / H,
    }


# ---------------------------------------------------------------------
# Data loading and templates
# ---------------------------------------------------------------------

def ensure_templates(data_dir: Path):
    data_dir.mkdir(parents=True, exist_ok=True)

    templates = {}

    sn_path = data_dir / "sn.csv"
    if not sn_path.exists():
        sn_path.write_text(
            "z,mu,sigma_mu\n"
            "0.01,33.18,0.12\n"
            "0.05,36.73,0.12\n"
            "0.10,38.31,0.12\n",
            encoding="utf-8"
        )
        templates["sn"] = sn_path

    bao_path = data_dir / "bao.csv"
    if not bao_path.exists():
        bao_path.write_text(
            "z,observable,value,sigma\n"
            "0.38,DM_over_rd,10.23,0.20\n"
            "0.38,DH_over_rd,24.9,0.75\n"
            "0.51,DM_over_rd,13.36,0.25\n"
            "0.51,DH_over_rd,22.3,0.67\n"
            "0.61,DM_over_rd,15.45,0.30\n"
            "0.61,DH_over_rd,20.1,0.60\n",
            encoding="utf-8"
        )
        templates["bao"] = bao_path

    cmb_path = data_dir / "cmb.csv"
    if not cmb_path.exists():
        cmb_path.write_text(
            "observable,value,sigma,zstar\n"
            "R_shift,1.75,0.02,1089\n"
            "DC_over_rd,94.0,0.5,1089\n",
            encoding="utf-8"
        )
        templates["cmb"] = cmb_path

    hz_path = data_dir / "chronometers.csv"
    if not hz_path.exists():
        hz_path.write_text(
            "z,H,sigma_H\n"
            "0.10,69.0,12.0\n"
            "0.20,72.0,10.0\n"
            "0.35,83.0,14.0\n"
            "0.60,97.0,15.0\n",
            encoding="utf-8"
        )
        templates["hz"] = hz_path

    return templates


def read_csv(path: Path):
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_sn(path: Path):
    rows = read_csv(path)
    out = []
    for r in rows:
        if not r.get("z") or not r.get("mu"):
            continue
        out.append({
            "z": float(r["z"]),
            "mu": float(r["mu"]),
            "sigma_mu": float(r.get("sigma_mu") or 0.15),
        })
    return out


def load_bao(path: Path):
    rows = read_csv(path)
    out = []
    for r in rows:
        if not r.get("z") or not r.get("observable") or not r.get("value") or not r.get("sigma"):
            continue
        out.append({
            "z": float(r["z"]),
            "observable": r["observable"].strip(),
            "value": float(r["value"]),
            "sigma": float(r["sigma"]),
        })
    return out


def load_cmb(path: Path):
    rows = read_csv(path)
    out = []
    for r in rows:
        if not r.get("observable") or not r.get("value") or not r.get("sigma"):
            continue
        out.append({
            "observable": r["observable"].strip(),
            "value": float(r["value"]),
            "sigma": float(r["sigma"]),
            "zstar": float(r.get("zstar") or 1089.0),
        })
    return out


def load_hz(path: Path):
    rows = read_csv(path)
    out = []
    for r in rows:
        if not r.get("z") or not r.get("H") or not r.get("sigma_H"):
            continue
        out.append({
            "z": float(r["z"]),
            "H": float(r["H"]),
            "sigma_H": float(r["sigma_H"]),
        })
    return out


def generate_synthetic_data(H0_ref=67.4, omega_m=0.30, r_d=147.1, seed=42):
    """Generate synthetic reference data from no-bias LCDM lower-H0 model."""
    rng = np.random.default_rng(seed)
    sn_z = np.concatenate([np.linspace(0.01, 0.15, 12), np.linspace(0.2, 1.6, 20)])
    dist = model_distances(sn_z, {"H0": H0_ref, "omega_m": omega_m}, "lcdm", "linear", "one_plus_z")
    sn = []
    for z, mu0 in zip(sn_z, distance_modulus(dist["D_L_obs"])):
        sig = 0.12
        sn.append({"z": float(z), "mu": float(mu0 + rng.normal(0, sig*0.15)), "sigma_mu": sig})

    bao_z = np.array([0.38, 0.51, 0.61, 1.48, 2.33])
    dist_b = model_distances(bao_z, {"H0": H0_ref, "omega_m": omega_m}, "lcdm", "linear", "one_plus_z")
    bao = []
    for z, dm, dh, dv in zip(bao_z, dist_b["D_M_obs"], dist_b["D_H"], dist_b["D_V_obs"]):
        bao.append({"z": float(z), "observable": "DM_over_rd", "value": float(dm/r_d), "sigma": float(0.02*dm/r_d)})
        bao.append({"z": float(z), "observable": "DH_over_rd", "value": float(dh/r_d), "sigma": float(0.03*dh/r_d)})
        bao.append({"z": float(z), "observable": "DV_over_rd", "value": float(dv/r_d), "sigma": float(0.025*dv/r_d)})

    zstar = 1089.0
    dist_c = model_distances(np.array([zstar]), {"H0": H0_ref, "omega_m": omega_m}, "lcdm", "linear", "one_plus_z")
    DCstar = dist_c["D_M_obs"][0]
    R_shift = math.sqrt(omega_m) * H0_ref * DCstar / c_km_s
    cmb = [
        {"observable": "R_shift", "value": R_shift, "sigma": 0.005*R_shift, "zstar": zstar},
        {"observable": "DC_over_rd", "value": DCstar/r_d, "sigma": 0.005*DCstar/r_d, "zstar": zstar},
    ]

    hz_z = np.array([0.1,0.2,0.35,0.5,0.7,1.0,1.5,2.0])
    hz = []
    for z, Hval in zip(hz_z, H_z(hz_z, H0_ref, omega_m)):
        sig = 0.07*Hval
        hz.append({"z": float(z), "H": float(Hval + rng.normal(0, sig*0.1)), "sigma_H": float(sig)})

    return sn, bao, cmb, hz


# ---------------------------------------------------------------------
# Likelihood
# ---------------------------------------------------------------------

def unpack_params(x, model_kind, free_rd, fit_Mcal):
    if model_kind == "lcdm":
        H0, om = x[0], x[1]
        idx = 2
        f0 = 0.0
        f1 = 0.0
    else:
        H0, om, f0, f1 = x[0], x[1], x[2], x[3]
        idx = 4

    r_d = x[idx] if free_rd else 147.1
    idx += 1 if free_rd else 0

    Mcal = x[idx] if fit_Mcal else 0.0

    return {
        "H0": float(H0),
        "omega_m": float(om),
        "f0": float(f0),
        "f1": float(f1),
        "r_d": float(r_d),
        "Mcal": float(Mcal),
    }


def bounds_for(model_kind, free_rd, fit_Mcal):
    if model_kind == "lcdm":
        b = [(50, 90), (0.1, 0.5)]
    else:
        b = [(50, 90), (0.1, 0.5), (-5, 10), (-20, 20)]
    if free_rd:
        b.append((120, 170))
    if fit_Mcal:
        b.append((-2, 2))
    return b


def chi2_for_dataset(params, model_kind, bias_shape_name, flos_model, sn, bao, cmb, hz, use):
    chi = {}
    n = {}

    if "SN" in use and sn:
        z = np.array([r["z"] for r in sn])
        mu_obs = np.array([r["mu"] for r in sn])
        sig = np.array([r["sigma_mu"] for r in sn])
        dist = model_distances(z, params, "stam" if model_kind == "stam" else "lcdm", bias_shape_name, flos_model)
        mu_model = distance_modulus(dist["D_L_obs"]) + params.get("Mcal", 0.0)
        chi["SN"] = float(np.sum(((mu_model - mu_obs)/sig)**2))
        n["SN"] = len(z)
    else:
        chi["SN"] = 0.0
        n["SN"] = 0

    if "BAO" in use and bao:
        csum = 0.0
        for r in bao:
            z = np.array([r["z"]])
            dist = model_distances(z, params, "stam" if model_kind == "stam" else "lcdm", bias_shape_name, flos_model)
            rd = params["r_d"]
            obs = r["observable"]
            if obs == "DM_over_rd":
                pred = dist["D_M_obs"][0] / rd
            elif obs == "DH_over_rd":
                pred = dist["D_H"][0] / rd
            elif obs == "DV_over_rd":
                pred = dist["D_V_obs"][0] / rd
            elif obs == "DL_over_rd":
                pred = dist["D_L_obs"][0] / rd
            else:
                continue
            csum += ((pred - r["value"])/r["sigma"])**2
        chi["BAO"] = float(csum)
        n["BAO"] = len(bao)
    else:
        chi["BAO"] = 0.0
        n["BAO"] = 0

    if "CMB" in use and cmb:
        csum = 0.0
        for r in cmb:
            zstar = np.array([r.get("zstar", 1089.0)])
            dist = model_distances(zstar, params, "stam" if model_kind == "stam" else "lcdm", bias_shape_name, flos_model)
            rd = params["r_d"]
            obs = r["observable"]
            if obs == "R_shift":
                pred = math.sqrt(params["omega_m"]) * params["H0"] * dist["D_M_obs"][0] / c_km_s
            elif obs == "lA":
                pred = math.pi * dist["D_M_obs"][0] / rd
            elif obs == "DC_over_rd":
                pred = dist["D_M_obs"][0] / rd
            else:
                continue
            csum += ((pred - r["value"])/r["sigma"])**2
        chi["CMB"] = float(csum)
        n["CMB"] = len(cmb)
    else:
        chi["CMB"] = 0.0
        n["CMB"] = 0

    if "HZ" in use and hz:
        z = np.array([r["z"] for r in hz])
        Hobs = np.array([r["H"] for r in hz])
        sig = np.array([r["sigma_H"] for r in hz])
        Hmod = H_z(z, params["H0"], params["omega_m"])
        chi["HZ"] = float(np.sum(((Hmod - Hobs)/sig)**2))
        n["HZ"] = len(hz)
    else:
        chi["HZ"] = 0.0
        n["HZ"] = 0

    chi["TOTAL"] = sum(chi.values())
    n["TOTAL"] = sum(n.values())
    return chi, n


def fit_model(model_kind, bias_shape_name, flos_model, sn, bao, cmb, hz, use, free_rd, fit_Mcal):
    bounds = bounds_for(model_kind, free_rd, fit_Mcal)

    def obj(x):
        p = unpack_params(x, model_kind, free_rd, fit_Mcal)
        # guard invalid
        if p["omega_m"] <= 0 or p["omega_m"] >= 1 or p["H0"] <= 0:
            return 1e99
        chi, n = chi2_for_dataset(p, model_kind, bias_shape_name, flos_model, sn, bao, cmb, hz, use)
        # mild physical penalty against strongly negative f_los over data range
        if model_kind == "stam":
            zcheck = np.array([0, 0.5, 2, 10, 1089.0])
            fvals = flos_z(zcheck, p["f0"], p["f1"], flos_model)
            if np.min(fvals) < -0.5:
                return chi["TOTAL"] + 1e6 * abs(np.min(fvals)+0.5)
        return chi["TOTAL"]

    de = differential_evolution(obj, bounds=bounds, seed=11, polish=False, tol=1e-7)
    local = minimize(obj, de.x, method="Nelder-Mead", options={"maxiter": 8000, "xatol": 1e-8, "fatol": 1e-8})
    xbest = local.x if local.fun <= de.fun else de.x
    params = unpack_params(xbest, model_kind, free_rd, fit_Mcal)
    chi, n = chi2_for_dataset(params, model_kind, bias_shape_name, flos_model, sn, bao, cmb, hz, use)
    k = len(xbest)
    aic = chi["TOTAL"] + 2*k
    bic = chi["TOTAL"] + k * math.log(max(n["TOTAL"], 1))
    return params, chi, n, aic, bic


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------

def plot_outputs(outdir, params_lcdm, params_stam, sn, bao, hz, bias_shape_name, flos_model):
    plots = outdir / "plots"
    plots.mkdir(parents=True, exist_ok=True)

    # SN residuals
    if sn:
        z = np.array([r["z"] for r in sn])
        mu_obs = np.array([r["mu"] for r in sn])
        sig = np.array([r["sigma_mu"] for r in sn])
        order = np.argsort(z)
        z = z[order]; mu_obs = mu_obs[order]; sig = sig[order]

        mu_lcdm = distance_modulus(model_distances(z, params_lcdm, "lcdm", bias_shape_name, flos_model)["D_L_obs"]) + params_lcdm.get("Mcal", 0.0)
        mu_stam = distance_modulus(model_distances(z, params_stam, "stam", bias_shape_name, flos_model)["D_L_obs"]) + params_stam.get("Mcal", 0.0)

        plt.figure(figsize=(10,6))
        plt.errorbar(z, mu_obs - mu_lcdm, yerr=sig, fmt="o", ms=3, alpha=0.6, label="data - LCDM")
        plt.plot(z, mu_stam - mu_lcdm, label="STAM - LCDM", lw=2)
        plt.axhline(0, color="black", lw=0.8)
        plt.xlabel("z")
        plt.ylabel("distance modulus residual")
        plt.title("G104 SN residuals")
        plt.grid(alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(plots / "G104_SN_hubble_residuals.png", dpi=200)
        plt.close()

    # H(z)
    if hz:
        z = np.array([r["z"] for r in hz])
        Hobs = np.array([r["H"] for r in hz])
        sig = np.array([r["sigma_H"] for r in hz])
        zgrid = np.linspace(0, max(2.0, np.max(z)*1.05), 300)
        plt.figure(figsize=(10,6))
        plt.errorbar(z, Hobs, yerr=sig, fmt="o", label="chronometers")
        plt.plot(zgrid, H_z(zgrid, params_lcdm["H0"], params_lcdm["omega_m"]), label="LCDM")
        plt.plot(zgrid, H_z(zgrid, params_stam["H0"], params_stam["omega_m"]), label="STAM intrinsic H(z)")
        plt.xlabel("z")
        plt.ylabel("H(z) km/s/Mpc")
        plt.title("G104 chronometer H(z)")
        plt.grid(alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(plots / "G104_Hz_chronometer_fit.png", dpi=200)
        plt.close()

    # f_los
    zgrid = np.linspace(0, 5, 300)
    plt.figure(figsize=(10,6))
    plt.plot(zgrid, flos_z(zgrid, params_stam.get("f0",0), params_stam.get("f1",0), flos_model), label="f_los(z)")
    plt.axhline(0, color="black", lw=0.8)
    plt.xlabel("z")
    plt.ylabel("f_los")
    plt.title("G104 best-fit STAM line-of-sight amplification")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plots / "G104_flos_shape.png", dpi=200)
    plt.close()


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="/mnt/data/G104_run")
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--synthetic", action="store_true", help="Force synthetic sanity data.")
    ap.add_argument("--real-only", action="store_true", help="Require real CSVs; fail if missing.")
    ap.add_argument("--bias-shape", choices=["linear", "stam_su", "saturating", "log"], default="linear")
    ap.add_argument("--flos-model", choices=["constant", "one_plus_z", "log"], default="one_plus_z")
    ap.add_argument("--free-rd", action="store_true")
    ap.add_argument("--no-fit-Mcal", action="store_true")
    ap.add_argument("--use", default="SN,BAO,CMB,HZ", help="Comma list: SN,BAO,CMB,HZ")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    results = outdir / "results"
    plots = outdir / "plots"
    data_dir = Path(args.data_dir) if args.data_dir else outdir / "data"
    results.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    templates = ensure_templates(data_dir)

    # Decide real vs synthetic.
    paths = {
        "sn": data_dir / "sn.csv",
        "bao": data_dir / "bao.csv",
        "cmb": data_dir / "cmb.csv",
        "hz": data_dir / "chronometers.csv",
    }
    using_templates = any(str(p) in [str(v) for v in templates.values()] for p in paths.values())

    if args.real_only and using_templates:
        raise FileNotFoundError(f"Real-only requested but missing data files. Templates written in {data_dir}")

    if args.synthetic or using_templates:
        sn, bao, cmb, hz = generate_synthetic_data()
        data_mode = "synthetic"
    else:
        sn = load_sn(paths["sn"])
        bao = load_bao(paths["bao"])
        cmb = load_cmb(paths["cmb"])
        hz = load_hz(paths["hz"])
        data_mode = "real_csv"

    use = [x.strip().upper() for x in args.use.split(",") if x.strip()]
    fit_Mcal = not args.no_fit_Mcal

    # Fit baseline and STAM.
    params_lcdm, chi_lcdm, n_lcdm, aic_lcdm, bic_lcdm = fit_model(
        "lcdm", args.bias_shape, args.flos_model, sn, bao, cmb, hz,
        use, args.free_rd, fit_Mcal
    )
    params_stam, chi_stam, n_stam, aic_stam, bic_stam = fit_model(
        "stam", args.bias_shape, args.flos_model, sn, bao, cmb, hz,
        use, args.free_rd, fit_Mcal
    )

    # Write results.
    param_csv = results / "G104_best_fit_params.csv"
    with param_csv.open("w", newline="", encoding="utf-8") as f:
        fields = ["model", "H0", "omega_m", "f0", "f1", "r_d", "Mcal", "chi2", "N", "chi2_per_N", "AIC", "BIC"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for name, p, chi, n, aic, bic in [
            ("lcdm", params_lcdm, chi_lcdm, n_lcdm, aic_lcdm, bic_lcdm),
            ("stam", params_stam, chi_stam, n_stam, aic_stam, bic_stam),
        ]:
            w.writerow({
                "model": name,
                "H0": p["H0"],
                "omega_m": p["omega_m"],
                "f0": p.get("f0", 0.0),
                "f1": p.get("f1", 0.0),
                "r_d": p.get("r_d", 147.1),
                "Mcal": p.get("Mcal", 0.0),
                "chi2": chi["TOTAL"],
                "N": n["TOTAL"],
                "chi2_per_N": chi["TOTAL"]/max(n["TOTAL"],1),
                "AIC": aic,
                "BIC": bic,
            })

    breakdown_csv = results / "G104_probe_chi2_breakdown.csv"
    with breakdown_csv.open("w", newline="", encoding="utf-8") as f:
        fields = ["model", "probe", "chi2", "N", "chi2_per_N"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for name, chi, n in [("lcdm", chi_lcdm, n_lcdm), ("stam", chi_stam, n_stam)]:
            for probe in ["SN", "BAO", "CMB", "HZ", "TOTAL"]:
                w.writerow({
                    "model": name,
                    "probe": probe,
                    "chi2": chi.get(probe, 0.0),
                    "N": n.get(probe, 0),
                    "chi2_per_N": chi.get(probe, 0.0)/max(n.get(probe, 1),1),
                })

    comparison_csv = results / "G104_model_comparison.csv"
    with comparison_csv.open("w", newline="", encoding="utf-8") as f:
        fields = ["quantity", "value"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerow({"quantity": "delta_chi2_LCDM_minus_STAM", "value": chi_lcdm["TOTAL"] - chi_stam["TOTAL"]})
        w.writerow({"quantity": "delta_AIC_LCDM_minus_STAM", "value": aic_lcdm - aic_stam})
        w.writerow({"quantity": "delta_BIC_LCDM_minus_STAM", "value": bic_lcdm - bic_stam})
        w.writerow({"quantity": "data_mode", "value": data_mode})

    plot_outputs(outdir, params_lcdm, params_stam, sn, bao, hz, args.bias_shape, args.flos_model)

    summary = results / "G104_fit_summary.md"
    md = []
    md.append("# G104 — SN + BAO + CMB + chronometer pipeline\n\n")
    md.append(f"- data mode: `{data_mode}`\n")
    md.append(f"- bias shape: `{args.bias_shape}`\n")
    md.append(f"- f_los model: `{args.flos_model}`\n")
    md.append(f"- free r_d: `{args.free_rd}`\n")
    md.append(f"- fit SN Mcal offset: `{fit_Mcal}`\n")
    md.append(f"- probes used: `{','.join(use)}`\n\n")
    if using_templates and not args.synthetic:
        md.append("Templates were created because one or more data files were missing. The run used synthetic sanity data. Replace templates in the data directory and rerun with `--real-only` for real-data mode.\n\n")

    md.append("## Best-fit parameters\n\n")
    md.append("| Model | H0 | Omega_m | f0 | f1 | r_d | Mcal | chi2 | N | chi2/N | AIC | BIC |\n")
    md.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
    for name, p, chi, n, aic, bic in [
        ("lcdm", params_lcdm, chi_lcdm, n_lcdm, aic_lcdm, bic_lcdm),
        ("stam", params_stam, chi_stam, n_stam, aic_stam, bic_stam),
    ]:
        md.append(f"| {name} | {p['H0']:.4f} | {p['omega_m']:.4f} | {p.get('f0',0):.4f} | {p.get('f1',0):.4f} | {p.get('r_d',147.1):.3f} | {p.get('Mcal',0):.4f} | {chi['TOTAL']:.3f} | {n['TOTAL']} | {chi['TOTAL']/max(n['TOTAL'],1):.3f} | {aic:.3f} | {bic:.3f} |\n")

    md.append("\n## Probe chi2 breakdown\n\n")
    md.append("| Model | SN | BAO | CMB | HZ | Total |\n")
    md.append("|---|---:|---:|---:|---:|---:|\n")
    for name, chi, n in [("lcdm", chi_lcdm, n_lcdm), ("stam", chi_stam, n_stam)]:
        md.append(f"| {name} | {chi.get('SN',0):.3f} | {chi.get('BAO',0):.3f} | {chi.get('CMB',0):.3f} | {chi.get('HZ',0):.3f} | {chi.get('TOTAL',0):.3f} |\n")

    md.append("\n## Interpretation\n\n")
    md.append("This pipeline tests whether the STAM photon-A distance layer can improve distance probes while direct expansion probes remain pressure tests. In synthetic mode, it should recover the lower-H0 reference model unless the STAM bias is strongly preferred by the generated data. In real-data mode, compare Δχ2, AIC, BIC, and probe-level residuals.\n\n")
    md.append("A positive STAM result requires more than better SN residuals: it must avoid wrecking BAO radial distances, chronometers, and CMB compressed distances under one locked convention.\n\n")
    md.append("## Files\n\n")
    md.append(f"- Best-fit params: `{param_csv}`\n")
    md.append(f"- Probe breakdown: `{breakdown_csv}`\n")
    md.append(f"- Model comparison: `{comparison_csv}`\n")
    md.append("- Plots in `plots/`\n")
    summary.write_text("".join(md), encoding="utf-8")

    print("=" * 100)
    print("G104: SN + BAO + CMB + chronometer pipeline")
    print("=" * 100)
    print(f"data mode: {data_mode}")
    print(f"data dir: {data_dir}")
    print()
    print("Best fits:")
    for name, p, chi, n, aic, bic in [
        ("LCDM", params_lcdm, chi_lcdm, n_lcdm, aic_lcdm, bic_lcdm),
        ("STAM", params_stam, chi_stam, n_stam, aic_stam, bic_stam),
    ]:
        print(f"{name}: H0={p['H0']:.4f}, Om={p['omega_m']:.4f}, "
              f"f0={p.get('f0',0):.4f}, f1={p.get('f1',0):.4f}, "
              f"chi2={chi['TOTAL']:.3f}, N={n['TOTAL']}, chi2/N={chi['TOTAL']/max(n['TOTAL'],1):.3f}, "
              f"AIC={aic:.3f}, BIC={bic:.3f}")
    print()
    print(f"Delta chi2 LCDM-STAM = {chi_lcdm['TOTAL'] - chi_stam['TOTAL']:.3f}")
    print(f"Delta AIC  LCDM-STAM = {aic_lcdm - aic_stam:.3f}")
    print(f"Delta BIC  LCDM-STAM = {bic_lcdm - bic_stam:.3f}")
    print()
    print(f"Summary written: {summary}")
    print(f"Templates/data directory: {data_dir}")
    print(f"Results directory: {results}")


if __name__ == "__main__":
    main()
