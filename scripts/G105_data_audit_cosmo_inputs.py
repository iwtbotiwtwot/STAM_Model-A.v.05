#!/usr/bin/env python3
"""
G105_data_audit_cosmo_inputs.py

Cosmology input-data audit for the STAM / Model-A G104 pipeline.

Purpose
=======
G104 can run in "real_csv" mode even when the files are tiny template-like CSVs.
G105 audits the actual input files before interpreting the cosmology fit.

It checks:

    data/sn.csv
    data/bao.csv
    data/cmb.csv
    data/chronometers.csv

and reports:

    - file presence
    - row counts
    - required-column checks
    - redshift ranges
    - observable counts
    - value/sigma ranges
    - likely-template warnings
    - per-row residuals for a baseline LCDM model and a STAM model if requested

This is a diagnostic utility, not a fit.

Usage
=====
Audit only:

    python G105_data_audit_cosmo_inputs.py --data-dir data --outdir G105_run

Audit against G104 best-fit params if available:

    python G105_data_audit_cosmo_inputs.py \
        --data-dir G104_run/data \
        --g104-params G104_run/results/G104_best_fit_params.csv \
        --outdir G105_run

Outputs
=======
results/G105_data_audit_summary.md
results/G105_file_audit.csv
results/G105_row_residuals.csv
plots/G105_data_coverage.png
plots/G105_residuals_by_probe.png
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Constants and model helpers — same conventions as G104
# ---------------------------------------------------------------------

c_km_s = 299_792.458
A0 = 1.0 / (12.0 * math.pi)


def E_flat_lcdm(z, omega_m):
    z = np.asarray(z, dtype=float)
    return np.sqrt(omega_m * (1.0 + z)**3 + (1.0 - omega_m))


def H_z(z, H0, omega_m):
    return float(H0) * E_flat_lcdm(z, omega_m)


def comoving_distance_interp(z_values, H0, omega_m):
    z_values = np.asarray(z_values, dtype=float)
    if z_values.size == 0:
        return np.array([])
    zmax = float(np.max(z_values))
    if zmax <= 0:
        return np.zeros_like(z_values)

    if zmax > 50:
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
    z = np.asarray(z, dtype=float)
    return (z * D_M**2 * c_km_s / H)**(1.0 / 3.0)


def distance_modulus(D_L_mpc):
    return 5.0 * np.log10(np.asarray(D_L_mpc, dtype=float)) + 25.0


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
    return z


def flos_z(z, f0, f1, model="one_plus_z"):
    z = np.asarray(z, dtype=float)
    if model == "constant":
        return f0 + 0.0 * z
    if model == "one_plus_z":
        return f0 + f1 * z / (1.0 + z)
    if model == "log":
        return f0 + f1 * np.log1p(z)
    return f0 + f1 * z / (1.0 + z)


def model_distances(z, params, model_kind, bias_shape_name="linear", flos_model="one_plus_z"):
    H0 = params["H0"]
    om = params["omega_m"]
    z = np.asarray(z, dtype=float)
    D_C = comoving_distance_interp(z, H0, om)

    if model_kind == "stam":
        f0 = params.get("f0", 0.0)
        f1 = params.get("f1", 0.0)
        bias = bridge_b_mpc(H0) * flos_z(z, f0, f1, flos_model) * bias_shape(z, bias_shape_name)
        D_M_obs = D_C + bias
    else:
        D_M_obs = D_C

    H = H_z(z, H0, om)
    D_L = D_L_from_D_C(z, D_M_obs)
    return {
        "D_M_obs": D_M_obs,
        "D_L_obs": D_L,
        "D_V_obs": D_V(z, D_M_obs, H),
        "H": H,
        "D_H": c_km_s / H,
    }


# ---------------------------------------------------------------------
# CSV utilities
# ---------------------------------------------------------------------

def read_csv(path: Path):
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def numeric_or_none(x):
    try:
        if x is None or str(x).strip() == "":
            return None
        return float(x)
    except Exception:
        return None


def minmax(vals):
    vals = [v for v in vals if v is not None and np.isfinite(v)]
    if not vals:
        return "", ""
    return min(vals), max(vals)


def is_template_like(file_kind, rows):
    # Heuristics based on the templates G104 writes.
    n = len(rows)
    if file_kind == "sn" and n <= 5:
        return True
    if file_kind == "bao" and n <= 6:
        return True
    if file_kind == "cmb" and n <= 2:
        return True
    if file_kind == "chronometers" and n <= 4:
        return True
    return False


def validate_file(kind, rows):
    if kind == "sn":
        required = ["z", "mu"]
    elif kind == "bao":
        required = ["z", "observable", "value", "sigma"]
    elif kind == "cmb":
        required = ["observable", "value", "sigma"]
    elif kind == "chronometers":
        required = ["z", "H", "sigma_H"]
    else:
        required = []

    missing = []
    if rows:
        cols = set(rows[0].keys())
        for c in required:
            if c not in cols:
                missing.append(c)
    else:
        missing = required

    return missing


def load_params_from_g104(path: Path):
    """Read G104_best_fit_params.csv if present."""
    if not path or not path.exists():
        return None, None
    rows = read_csv(path)
    lcdm = None
    stam = None
    for r in rows:
        p = {
            "H0": float(r["H0"]),
            "omega_m": float(r["omega_m"]),
            "f0": float(r.get("f0", 0.0) or 0.0),
            "f1": float(r.get("f1", 0.0) or 0.0),
            "r_d": float(r.get("r_d", 147.1) or 147.1),
            "Mcal": float(r.get("Mcal", 0.0) or 0.0),
        }
        if r.get("model") == "lcdm":
            lcdm = p
        elif r.get("model") == "stam":
            stam = p
    return lcdm, stam


# ---------------------------------------------------------------------
# Residual audit
# ---------------------------------------------------------------------

def compute_residuals(sn, bao, cmb, hz, lcdm_params, stam_params, bias_shape_name, flos_model):
    residuals = []

    if lcdm_params is None:
        lcdm_params = {"H0": 67.4, "omega_m": 0.30, "r_d": 147.1, "Mcal": 0.0}
    if stam_params is None:
        stam_params = {"H0": 73.04, "omega_m": 0.30, "f0": 0.0, "f1": 0.0, "r_d": 147.1, "Mcal": 0.0}

    # SN residuals
    for i, r in enumerate(sn):
        z = numeric_or_none(r.get("z"))
        obs = numeric_or_none(r.get("mu"))
        sig = numeric_or_none(r.get("sigma_mu")) or 0.15
        if z is None or obs is None:
            continue
        zarr = np.array([z])
        pred_lcdm = distance_modulus(model_distances(zarr, lcdm_params, "lcdm", bias_shape_name, flos_model)["D_L_obs"])[0] + lcdm_params.get("Mcal", 0.0)
        pred_stam = distance_modulus(model_distances(zarr, stam_params, "stam", bias_shape_name, flos_model)["D_L_obs"])[0] + stam_params.get("Mcal", 0.0)
        residuals.append({
            "probe": "SN", "row": i, "z": z, "observable": "mu", "value": obs, "sigma": sig,
            "lcdm_pred": pred_lcdm, "stam_pred": pred_stam,
            "lcdm_pull": (pred_lcdm - obs)/sig, "stam_pull": (pred_stam - obs)/sig
        })

    # BAO residuals
    for i, r in enumerate(bao):
        z = numeric_or_none(r.get("z"))
        obsname = str(r.get("observable", "")).strip()
        obs = numeric_or_none(r.get("value"))
        sig = numeric_or_none(r.get("sigma"))
        if z is None or obs is None or sig is None:
            continue
        zarr = np.array([z])
        dist_l = model_distances(zarr, lcdm_params, "lcdm", bias_shape_name, flos_model)
        dist_s = model_distances(zarr, stam_params, "stam", bias_shape_name, flos_model)
        rd_l = lcdm_params.get("r_d", 147.1)
        rd_s = stam_params.get("r_d", 147.1)

        def pred(dist, rd):
            if obsname == "DM_over_rd":
                return dist["D_M_obs"][0]/rd
            if obsname == "DH_over_rd":
                return dist["D_H"][0]/rd
            if obsname == "DV_over_rd":
                return dist["D_V_obs"][0]/rd
            if obsname == "DL_over_rd":
                return dist["D_L_obs"][0]/rd
            return None

        pl = pred(dist_l, rd_l)
        ps = pred(dist_s, rd_s)
        if pl is None or ps is None:
            continue

        residuals.append({
            "probe": "BAO", "row": i, "z": z, "observable": obsname, "value": obs, "sigma": sig,
            "lcdm_pred": pl, "stam_pred": ps,
            "lcdm_pull": (pl - obs)/sig, "stam_pull": (ps - obs)/sig
        })

    # CMB residuals
    for i, r in enumerate(cmb):
        obsname = str(r.get("observable", "")).strip()
        obs = numeric_or_none(r.get("value"))
        sig = numeric_or_none(r.get("sigma"))
        zstar = numeric_or_none(r.get("zstar")) or 1089.0
        if obs is None or sig is None:
            continue
        zarr = np.array([zstar])
        dist_l = model_distances(zarr, lcdm_params, "lcdm", bias_shape_name, flos_model)
        dist_s = model_distances(zarr, stam_params, "stam", bias_shape_name, flos_model)

        def pred(dist, p):
            rd = p.get("r_d", 147.1)
            if obsname == "R_shift":
                return math.sqrt(p["omega_m"]) * p["H0"] * dist["D_M_obs"][0] / c_km_s
            if obsname == "lA":
                return math.pi * dist["D_M_obs"][0] / rd
            if obsname == "DC_over_rd":
                return dist["D_M_obs"][0] / rd
            return None

        pl = pred(dist_l, lcdm_params)
        ps = pred(dist_s, stam_params)
        if pl is None or ps is None:
            continue
        residuals.append({
            "probe": "CMB", "row": i, "z": zstar, "observable": obsname, "value": obs, "sigma": sig,
            "lcdm_pred": pl, "stam_pred": ps,
            "lcdm_pull": (pl - obs)/sig, "stam_pull": (ps - obs)/sig
        })

    # Chronometer residuals
    for i, r in enumerate(hz):
        z = numeric_or_none(r.get("z"))
        obs = numeric_or_none(r.get("H"))
        sig = numeric_or_none(r.get("sigma_H"))
        if z is None or obs is None or sig is None:
            continue
        pl = H_z(np.array([z]), lcdm_params["H0"], lcdm_params["omega_m"])[0]
        ps = H_z(np.array([z]), stam_params["H0"], stam_params["omega_m"])[0]
        residuals.append({
            "probe": "HZ", "row": i, "z": z, "observable": "H", "value": obs, "sigma": sig,
            "lcdm_pred": pl, "stam_pred": ps,
            "lcdm_pull": (pl - obs)/sig, "stam_pull": (ps - obs)/sig
        })

    return residuals


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="G105 cosmology input-data audit")
    ap.add_argument("--data-dir", default="/mnt/data/G104_run/data")
    ap.add_argument("--outdir", default="/mnt/data/G105_run")
    ap.add_argument("--g104-params", default="/mnt/data/G104_run/results/G104_best_fit_params.csv")
    ap.add_argument("--bias-shape", default="linear")
    ap.add_argument("--flos-model", default="one_plus_z")
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    outdir = Path(args.outdir)
    results_dir = outdir / "results"
    plots_dir = outdir / "plots"
    results_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    files = {
        "sn": data_dir / "sn.csv",
        "bao": data_dir / "bao.csv",
        "cmb": data_dir / "cmb.csv",
        "chronometers": data_dir / "chronometers.csv",
    }

    data = {k: read_csv(p) for k, p in files.items()}
    lcdm_params, stam_params = load_params_from_g104(Path(args.g104_params))

    # File audit rows
    audit_rows = []
    for kind, path in files.items():
        rows = data[kind]
        missing = validate_file(kind, rows)
        zvals = []
        sigvals = []
        valvals = []

        for r in rows:
            if "z" in r:
                zvals.append(numeric_or_none(r.get("z")))
            if kind == "sn":
                valvals.append(numeric_or_none(r.get("mu")))
                sigvals.append(numeric_or_none(r.get("sigma_mu")))
            elif kind == "bao":
                valvals.append(numeric_or_none(r.get("value")))
                sigvals.append(numeric_or_none(r.get("sigma")))
            elif kind == "cmb":
                valvals.append(numeric_or_none(r.get("value")))
                sigvals.append(numeric_or_none(r.get("sigma")))
                zvals.append(numeric_or_none(r.get("zstar")))
            elif kind == "chronometers":
                valvals.append(numeric_or_none(r.get("H")))
                sigvals.append(numeric_or_none(r.get("sigma_H")))

        zmin, zmax = minmax(zvals)
        vmin, vmax = minmax(valvals)
        smin, smax = minmax(sigvals)

        observable_counts = ""
        if kind in ("bao", "cmb") and rows:
            observable_counts = dict(Counter(r.get("observable", "").strip() for r in rows))

        audit_rows.append({
            "file_kind": kind,
            "path": str(path),
            "exists": path.exists(),
            "row_count": len(rows),
            "missing_required_columns": ";".join(missing),
            "z_min": zmin,
            "z_max": zmax,
            "value_min": vmin,
            "value_max": vmax,
            "sigma_min": smin,
            "sigma_max": smax,
            "observable_counts": str(observable_counts),
            "template_like_warning": is_template_like(kind, rows),
        })

    # Residuals
    residuals = compute_residuals(
        data["sn"], data["bao"], data["cmb"], data["chronometers"],
        lcdm_params, stam_params,
        args.bias_shape, args.flos_model
    )

    # Write CSVs
    audit_csv = results_dir / "G105_file_audit.csv"
    with audit_csv.open("w", newline="", encoding="utf-8") as f:
        fields = list(audit_rows[0].keys())
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(audit_rows)

    residual_csv = results_dir / "G105_row_residuals.csv"
    with residual_csv.open("w", newline="", encoding="utf-8") as f:
        if residuals:
            fields = list(residuals[0].keys())
        else:
            fields = ["probe", "row", "z", "observable", "value", "sigma", "lcdm_pred", "stam_pred", "lcdm_pull", "stam_pull"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in residuals:
            w.writerow(r)

    # Probe summary
    probe_summary = defaultdict(lambda: {"N": 0, "lcdm_chi2": 0.0, "stam_chi2": 0.0})
    for r in residuals:
        p = r["probe"]
        probe_summary[p]["N"] += 1
        probe_summary[p]["lcdm_chi2"] += r["lcdm_pull"]**2
        probe_summary[p]["stam_chi2"] += r["stam_pull"]**2

    # Plots
    # Data coverage
    plt.figure(figsize=(9, 5))
    kinds = [r["file_kind"] for r in audit_rows]
    counts = [r["row_count"] for r in audit_rows]
    colors = ["tab:red" if r["template_like_warning"] else "tab:blue" for r in audit_rows]
    plt.bar(kinds, counts, color=colors)
    plt.ylabel("row count")
    plt.title("G105 cosmology input-data coverage")
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    coverage_plot = plots_dir / "G105_data_coverage.png"
    plt.savefig(coverage_plot, dpi=200)
    plt.close()

    # Pulls by probe
    if residuals:
        probes = sorted(set(r["probe"] for r in residuals))
        x = np.arange(len(probes))
        lcdm_chi = [probe_summary[p]["lcdm_chi2"]/max(probe_summary[p]["N"],1) for p in probes]
        stam_chi = [probe_summary[p]["stam_chi2"]/max(probe_summary[p]["N"],1) for p in probes]
        width = 0.35
        plt.figure(figsize=(9, 5))
        plt.bar(x-width/2, lcdm_chi, width, label="LCDM")
        plt.bar(x+width/2, stam_chi, width, label="STAM")
        plt.xticks(x, probes)
        plt.ylabel("chi2 / N")
        plt.title("G105 residuals by probe")
        plt.grid(axis="y", alpha=0.3)
        plt.legend()
        plt.tight_layout()
        residual_plot = plots_dir / "G105_residuals_by_probe.png"
        plt.savefig(residual_plot, dpi=200)
        plt.close()
    else:
        residual_plot = None

    # Console output
    print("=" * 100)
    print("G105: cosmology input-data audit")
    print("=" * 100)
    print(f"data dir: {data_dir}")
    print()
    print(f"{'file':<16}{'exists':>8}{'rows':>8}{'z range':>20}{'template?':>12}{'missing columns':>24}")
    print("-" * 92)
    for r in audit_rows:
        zrange = f"{r['z_min']}..{r['z_max']}"
        print(f"{r['file_kind']:<16}{str(r['exists']):>8}{r['row_count']:>8}{zrange:>20}{str(r['template_like_warning']):>12}{r['missing_required_columns']:>24}")
        if r["observable_counts"]:
            print(f"  observables: {r['observable_counts']}")

    print()
    print("Probe residual summary:")
    print(f"{'probe':<8}{'N':>6}{'LCDM chi2/N':>16}{'STAM chi2/N':>16}")
    print("-" * 50)
    for p, s in sorted(probe_summary.items()):
        print(f"{p:<8}{s['N']:>6}{s['lcdm_chi2']/max(s['N'],1):>16.3f}{s['stam_chi2']/max(s['N'],1):>16.3f}")

    # Markdown summary
    summary = results_dir / "G105_data_audit_summary.md"
    md = []
    md.append("# G105 — Cosmology input-data audit\n\n")
    md.append(f"Data directory: `{data_dir}`\n\n")
    md.append("## File audit\n\n")
    md.append("| File | Exists | Rows | z range | Template-like warning | Missing columns | Observable counts |\n")
    md.append("|---|---:|---:|---:|---:|---|---|\n")
    for r in audit_rows:
        md.append(f"| {r['file_kind']} | {r['exists']} | {r['row_count']} | {r['z_min']}..{r['z_max']} | {r['template_like_warning']} | {r['missing_required_columns']} | {r['observable_counts']} |\n")

    md.append("\n## Probe residual summary\n\n")
    md.append("| Probe | N | LCDM chi2/N | STAM chi2/N |\n")
    md.append("|---|---:|---:|---:|\n")
    for p, s in sorted(probe_summary.items()):
        md.append(f"| {p} | {s['N']} | {s['lcdm_chi2']/max(s['N'],1):.3f} | {s['stam_chi2']/max(s['N'],1):.3f} |\n")

    md.append("\n## Interpretation guide\n\n")
    md.append("- If a file is flagged `template-like`, the corresponding G104 result is not a serious real-data test.\n")
    md.append("- `DH_over_rd` BAO and chronometer `H(z)` are direct expansion probes; they are expected pressure points for a pure photon-distance-bias model.\n")
    md.append("- A serious real-data cosmology run should use many SN rows, modern BAO compressed values with covariance where possible, real CMB compressed constraints, and a full chronometer table.\n\n")
    md.append("## Files\n\n")
    md.append(f"- File audit CSV: `{audit_csv}`\n")
    md.append(f"- Row residual CSV: `{residual_csv}`\n")
    md.append(f"- Coverage plot: `{coverage_plot}`\n")
    if residual_plot:
        md.append(f"- Residual plot: `{residual_plot}`\n")
    summary.write_text("".join(md), encoding="utf-8")

    print()
    print(f"Summary written: {summary}")
    print(f"File audit CSV: {audit_csv}")
    print(f"Row residual CSV: {residual_csv}")
    print(f"Coverage plot: {coverage_plot}")
    if residual_plot:
        print(f"Residual plot: {residual_plot}")


if __name__ == "__main__":
    main()
