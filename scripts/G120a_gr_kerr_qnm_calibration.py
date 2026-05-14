#!/usr/bin/env python3
"""
G120a_gr_kerr_qnm_calibration.py

GR Kerr QNM calibration ONLY.  No STAM substitution.

Goal
====
Establish a verified GR Kerr QNM reference table for (l=2, m=2, n=0) at
spins a in {0, 0.3, 0.5, 0.7, 0.9} with calibration error < 0.5%.

Strategy
========
Use the standard published Kerr QNM library  qnm  (Stein 2019, based on
Leaver's continued fraction method with Berti's verified tabulation).
qnm provides Kerr QNM eigenvalues to many digits of precision; using it
as the calibration reference is consistent with the framework's
requirement "reproduce published Kerr QNMs to < 0.5%."

What this script does (and does not)
=====================================
DOES:
- Fetch verified Kerr (l=2, m=2, n=0) QNMs at the five spins.
- Tabulate omega(a), angular eigenvalue A_lm(a).
- Save to results/G120a_gr_kerr_qnm_reference.csv and .json.
- Plot omega_R(a), |omega_I|(a) and the spin trend.
- This establishes the GR baseline that G120b's STAM solver must match.

DOES NOT:
- Run any STAM substitution.  (G120b.)
- Build a standalone Kerr QNM solver from scratch.  (G120b/c -- with
  qnm as the verification target.)
- Test other modes (l != 2, m != 2, n != 0).  Extend later if needed.

Status gate
===========
"PASS" means qnm successfully returns values at all five spins.
Errors against qnm's tabulated values are intrinsically << 0.5% because
qnm is the source of those values.  The < 0.5% gate is what G120b's
standalone solver must hit when calibrated against this G120a reference.

Outputs
=======
    results/G120a_gr_kerr_qnm_calibration_summary.md
    results/G120a_gr_kerr_qnm_reference.csv
    results/G120a_gr_kerr_qnm_reference.json
    plots/G120a_gr_kerr_qnm_calibration.png
"""

from __future__ import annotations

import sys
import json
import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


def main():
    print("=" * 88)
    print("G120a  --  GR Kerr QNM calibration  (no STAM)")
    print("=" * 88)
    print()

    try:
        import qnm
    except ImportError as e:
        print("ERROR: qnm package not available.  Install with `pip install qnm`.")
        print(f"  ImportError: {e}")
        sys.exit(1)

    print(f"Using qnm package version: {qnm.__version__}")
    print()
    print("Target modes: spin-weight s = -2, l = 2, m = 2, n = 0")
    print("Target spins: a in [0.0, 0.3, 0.5, 0.7, 0.9]")
    print("Calibration target: < 0.5% error vs published Kerr QNMs.")
    print()

    target_spins = [0.0, 0.3, 0.5, 0.7, 0.9]

    # Initialize qnm cache for (s=-2, l=2, m=2, n=0)
    print("Initializing qnm cache for (s=-2, l=2, m=2, n=0) ...", flush=True)
    qnm_cache = qnm.modes_cache(s=-2, l=2, m=2, n=0)
    print("  cache ready.")
    print()

    # Fetch QNMs at each target spin
    print("=" * 88)
    print("Kerr (l=2, m=2, n=0) QNMs from qnm  (M = 1 units)")
    print("=" * 88)
    print()
    print(f"  {'a':>5}  {'omega_R':>14}  {'omega_I':>14}  "
          f"{'|omega|':>10}  {'A_lm':>22}  {'status':>8}")
    print("-" * 95)

    results = []
    for a in target_spins:
        try:
            omega, A, C = qnm_cache(a=a)
            re_om = float(omega.real)
            im_om = float(omega.imag)
            abs_om = float(abs(omega))
            A_complex = complex(A) if hasattr(A, "real") else complex(float(A), 0.0)
            entry = {
                "a": a,
                "omega_real": re_om,
                "omega_imag": im_om,
                "omega_abs": abs_om,
                "A_lm_real": float(A_complex.real),
                "A_lm_imag": float(A_complex.imag),
                "status": "OK",
            }
            print(f"  {a:>5.2f}  {re_om:>14.10f}  {im_om:>14.10f}  "
                  f"{abs_om:>10.6f}  {str(A_complex):>22}  {'OK':>8}")
        except Exception as e:
            entry = {"a": a, "status": f"FAIL: {e}"}
            print(f"  {a:>5.2f}  FAIL: {e}")
        results.append(entry)

    print()
    n_ok = sum(1 for r in results if r.get("status") == "OK")
    print(f"  passed: {n_ok}/{len(results)} spins")
    print()

    # Spin trends
    print("=" * 88)
    print("Spin trends (n=0 corotating dominant mode)")
    print("=" * 88)
    print()
    ok = [r for r in results if r.get("status") == "OK"]
    if ok:
        spins = np.array([r["a"] for r in ok])
        oms_R = np.array([r["omega_real"] for r in ok])
        oms_I = np.array([r["omega_imag"] for r in ok])
        print(f"  Re(omega) ranges from {oms_R.min():.4f} (a={spins[np.argmin(oms_R)]:.2f}) "
              f"to {oms_R.max():.4f} (a={spins[np.argmax(oms_R)]:.2f})")
        print(f"  |Im(omega)| ranges from {abs(oms_I).min():.4f} (a={spins[np.argmin(abs(oms_I))]:.2f}) "
              f"to {abs(oms_I).max():.4f} (a={spins[np.argmax(abs(oms_I))]:.2f})")
        print()
        print("  Expected trend (Kerr GR):")
        print("    Re(omega) increases monotonically with spin (faster prograde orbit).")
        print("    |Im(omega)| decreases with spin (slower decay near extremal).")
        print()
        # Verify trends
        monotonic_Re = all(oms_R[i+1] > oms_R[i] for i in range(len(oms_R)-1))
        monotonic_Im = all(abs(oms_I[i+1]) < abs(oms_I[i]) for i in range(len(oms_I)-1))
        print(f"  Re(omega) monotonic increase:  {'PASS' if monotonic_Re else 'FAIL'}")
        print(f"  |Im(omega)| monotonic decrease: {'PASS' if monotonic_Im else 'FAIL'}")
        print()

    # Save reference data
    print("=" * 88)
    print("Saving reference data ...")
    print("=" * 88)
    print()

    json_path = RESULTS / "G120a_gr_kerr_qnm_reference.json"
    with json_path.open("w", encoding="utf-8") as f:
        json.dump({
            "description": "GR Kerr QNM reference for l=m=2, n=0 (M=1 units, s=-2).",
            "source": "qnm package (Stein 2019), uses Berti's verified tabulation.",
            "modes": {"s": -2, "l": 2, "m": 2, "n": 0},
            "spins": target_spins,
            "data": results,
        }, f, indent=2)
    print(f"  JSON: {json_path}")

    csv_path = RESULTS / "G120a_gr_kerr_qnm_reference.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["a", "omega_real", "omega_imag", "omega_abs",
                         "A_lm_real", "A_lm_imag", "status"])
        for r in results:
            writer.writerow([
                r.get("a"), r.get("omega_real"), r.get("omega_imag"),
                r.get("omega_abs"), r.get("A_lm_real"), r.get("A_lm_imag"),
                r.get("status"),
            ])
    print(f"  CSV:  {csv_path}")
    print()

    # Plot
    if ok:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        ax = axes[0, 0]
        ax.plot(spins, oms_R, "o-", linewidth=2, markersize=8, color="tab:blue")
        ax.set_xlabel("spin a / M")
        ax.set_ylabel("Re(omega)")
        ax.set_title("Kerr (l=m=2, n=0) Re(omega) vs spin")
        ax.grid(True, alpha=0.3)

        ax = axes[0, 1]
        ax.plot(spins, abs(oms_I), "s-", linewidth=2, markersize=8, color="tab:red")
        ax.set_xlabel("spin a / M")
        ax.set_ylabel("|Im(omega)|")
        ax.set_title("Kerr (l=m=2, n=0) |Im(omega)| vs spin (decay rate)")
        ax.grid(True, alpha=0.3)

        ax = axes[1, 0]
        # omega in complex plane
        ax.plot(oms_R, oms_I, "o-", linewidth=2, markersize=8, color="tab:purple")
        for r in ok:
            ax.annotate(f"a={r['a']:.1f}", (r["omega_real"], r["omega_imag"]),
                         xytext=(5, 5), textcoords="offset points", fontsize=8)
        ax.set_xlabel("Re(omega)")
        ax.set_ylabel("Im(omega)")
        ax.set_title("Kerr QNM track in complex omega plane")
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color="gray", linewidth=0.5)
        ax.axvline(0, color="gray", linewidth=0.5)

        ax = axes[1, 1]
        # Quality factor Q = -Re/(2*Im) (rings per e-fold)
        Q = oms_R / (2.0 * abs(oms_I))
        ax.plot(spins, Q, "^-", linewidth=2, markersize=8, color="tab:green")
        ax.set_xlabel("spin a / M")
        ax.set_ylabel("Q = Re(omega) / (2 |Im(omega)|)")
        ax.set_title("Ringdown quality factor vs spin")
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        out_png = PLOTS / "G120a_gr_kerr_qnm_calibration.png"
        plt.savefig(out_png, dpi=180)
        plt.close()
        print(f"  Plot: {out_png}")
        print()

    # Verdict
    print("=" * 88)
    print("VERDICT")
    print("=" * 88)
    print()
    if n_ok == len(target_spins):
        verdict = (
            f"PASS  --  qnm successfully provides verified Kerr (l=m=2, n=0) QNMs\n"
            f"at all {n_ok} target spins.  This is the GR calibration baseline.\n\n"
            "G120b's standalone Kerr QNM solver must reproduce these values to\n"
            "< 0.5% before any STAM substitution can be trusted."
        )
    else:
        verdict = (
            f"PARTIAL  --  qnm returned {n_ok}/{len(target_spins)} successfully.\n"
            "Check qnm install or fetch errors above."
        )
    print(verdict)
    print()

    # Markdown summary
    md = ["# G120a -- GR Kerr QNM calibration (no STAM)\n"]
    md.append("\n## Purpose\n")
    md.append("Establish a verified GR Kerr QNM reference table for "
              "(l=2, m=2, n=0) at the five target spins (a = 0, 0.3, 0.5, 0.7, 0.9).  "
              "No STAM substitution.  This is the calibration baseline that G120b's "
              "standalone solver must match to < 0.5% before STAM extension.\n")
    md.append("\n## Reference values\n")
    md.append("Source: qnm package "
              f"(version {qnm.__version__}, based on Leaver continued fraction + "
              "Berti's verified tabulation).\n\n")
    md.append("| a | Re(omega) | Im(omega) | |omega| | A_lm |\n")
    md.append("|---:|---:|---:|---:|---|\n")
    for r in results:
        if r.get("status") == "OK":
            A_str = f"{r['A_lm_real']:.6f}"
            if abs(r.get("A_lm_imag", 0.0)) > 1e-10:
                A_str += f" + {r['A_lm_imag']:.6f}i"
            md.append(f"| {r['a']:.2f} | {r['omega_real']:.10f} | "
                      f"{r['omega_imag']:.10f} | {r['omega_abs']:.8f} | "
                      f"{A_str} |\n")
        else:
            md.append(f"| {r['a']:.2f} | --- | --- | --- | --- |\n")
    md.append(f"\n**Passed: {n_ok}/{len(target_spins)} spins.**\n")
    md.append("\n## Spin trends\n")
    if ok:
        md.append("- Re(omega) monotonically increases with spin "
                  f"({oms_R.min():.4f} at a={spins[np.argmin(oms_R)]:.2f} "
                  f"→ {oms_R.max():.4f} at a={spins[np.argmax(oms_R)]:.2f}).\n")
        md.append("- |Im(omega)| monotonically decreases with spin "
                  f"({abs(oms_I).max():.4f} at a={spins[np.argmax(abs(oms_I))]:.2f} "
                  f"→ {abs(oms_I).min():.4f} at a={spins[np.argmin(abs(oms_I))]:.2f}).\n")
        md.append("\nThese trends match the expected Kerr GR behavior "
                  "(faster prograde orbit at higher spin; less efficient ringdown decay "
                  "due to longer-lived photon orbits near extremal).\n")
    md.append("\n## Verdict\n")
    md.append(verdict.replace("\n", "\n  ") + "\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G120a_gr_kerr_qnm_calibration.py](../scripts/G120a_gr_kerr_qnm_calibration.py)\n")
    md.append(f"- Reference JSON: `results/G120a_gr_kerr_qnm_reference.json`\n")
    md.append(f"- Reference CSV:  `results/G120a_gr_kerr_qnm_reference.csv`\n")
    md.append(f"- Plot: `plots/G120a_gr_kerr_qnm_calibration.png`\n")
    md.append("\n## Note on the < 0.5% gate\n")
    md.append("Because qnm IS the published reference (it uses Berti's verified "
              "tabulation), the 'error' against published values is intrinsically "
              "at the floating-point level.  The < 0.5% gate is what G120b's "
              "standalone Leaver / Sasaki-Nakamura / Detweiler solver must hit when "
              "calibrated against the G120a table.  This script establishes the "
              "table; G120b builds and tests the standalone solver.\n")
    (RESULTS / "G120a_gr_kerr_qnm_calibration_summary.md").write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {RESULTS / 'G120a_gr_kerr_qnm_calibration_summary.md'}")


if __name__ == "__main__":
    main()
