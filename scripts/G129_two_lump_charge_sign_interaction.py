#!/usr/bin/env python3
"""
G129 — Two-lump charge-sign interaction

Purpose
-------
Next gate after G128B.

Test whether the constrained Gauss-law U(1)-style seed gives the first
electrostatic behavior:

    like signs repel
    opposite signs attract
    field energy depends on separation
    total signed charge remains controlled
    Gauss constraint remains satisfied

Scope
-----
This is a reduced 2D electrostatic sandbox. It is NOT full electromagnetism,
QED, photons, or the Standard Model.

Model
-----
Two matter-density lumps are assigned charge signs s_i = ±1.

Charge density:
    rho_q = e * sum_i s_i * |psi_i|^2
with background mean-subtraction for periodic solvability.

Gauss law:
    (-∇² + μ²) φ = rho_q - mean(rho_q)

Electrostatic energy:
    E_phi = 1/2 ∫ [(∇φ)^2 + μ² φ²] dA
Equivalent interaction diagnostic:
    E_int = 1/2 ∫ (rho_q - mean(rho_q)) φ dA

Force on each lump from the field:
    F_i ≈ - q_i ∫ |psi_i|² ∇φ dA

Expected:
    like signs: force drives separation outward
    opposite signs: force drives separation inward

We do a static force/energy scan, not a full two-body dynamical collision.
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
PLOTS = OUT / "plots"
PLOTS.mkdir(parents=True, exist_ok=True)

# Grid
N = 256
L = 28.0
dx = L / N
x = (np.arange(N) - N//2) * dx
y = (np.arange(N) - N//2) * dx
X, Y = np.meshgrid(x, y, indexing="ij")

# Fourier grid
kx = 2*np.pi*np.fft.fftfreq(N, d=dx)
ky = 2*np.pi*np.fft.fftfreq(N, d=dx)
KX, KY = np.meshgrid(kx, ky, indexing="ij")
K2 = KX**2 + KY**2

# Parameters
e_charge = 1.0
mu = 0.28
sigma = 1.25
amp = 1.0
separations = np.linspace(3.0, 11.0, 17)

params = {
    "N": int(N),
    "L": float(L),
    "dx": float(dx),
    "e_charge": float(e_charge),
    "mu_screening": float(mu),
    "sigma": float(sigma),
    "amp": float(amp),
    "separations": [float(s) for s in separations],
    "scope": "2D static constrained Gauss-law two-lump charge-sign sandbox"
}

def grad(f):
    return ((np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0)) / (2*dx),
            (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1)) / (2*dx))

def lap_fft(f):
    return np.real(np.fft.ifft2(-K2*np.fft.fft2(f)))

def gaussian_lump(x0, y0):
    return amp**2 * np.exp(-((X-x0)**2 + (Y-y0)**2)/(2*sigma**2))

def solve_phi(rho_q):
    src = rho_q - np.mean(rho_q)
    phi_k = np.fft.fft2(src) / (K2 + mu**2)
    return np.real(np.fft.ifft2(phi_k))

def gauss_residual(phi, rho_q):
    src = rho_q - np.mean(rho_q)
    lhs = -lap_fft(phi) + mu**2*phi
    res = lhs - src
    return float(np.sqrt(np.mean(res**2))), float(np.sqrt(np.mean(src**2)) + 1e-30)

def field_energy(phi):
    phix, phiy = grad(phi)
    return float(0.5*np.sum(phix**2 + phiy**2 + mu**2*phi**2)*dx*dx)

def interaction_energy(phi, rho_q):
    src = rho_q - np.mean(rho_q)
    return float(0.5*np.sum(src*phi)*dx*dx)

def force_on_lump(phi, density, sign):
    phix, phiy = grad(phi)
    # Potential energy contribution q density phi -> force = -q density grad phi
    Fx = -e_charge*sign*np.sum(density*phix)*dx*dx
    Fy = -e_charge*sign*np.sum(density*phiy)*dx*dx
    return float(Fx), float(Fy)

def center_of_mass_density(density):
    tot = np.sum(density)+1e-30
    return float(np.sum(density*X)/tot), float(np.sum(density*Y)/tot)

def run_pair(label, s1, s2):
    rows = []
    for d in separations:
        rho1 = gaussian_lump(-d/2, 0.0)
        rho2 = gaussian_lump( d/2, 0.0)
        rho_q = e_charge*(s1*rho1 + s2*rho2)
        phi = solve_phi(rho_q)
        gres, gsrc = gauss_residual(phi, rho_q)
        Efield = field_energy(phi)
        Eint = interaction_energy(phi, rho_q)
        Fx1, Fy1 = force_on_lump(phi, rho1, s1)
        Fx2, Fy2 = force_on_lump(phi, rho2, s2)
        x1, y1 = center_of_mass_density(rho1)
        x2, y2 = center_of_mass_density(rho2)
        sep_vec_x = x2-x1

        # separation-rate force diagnostic:
        # d(separation)/dt increases if right lump force positive and left lump force negative.
        sep_force = Fx2 - Fx1
        expected = "repel" if s1*s2 > 0 else "attract"
        expected_sign = 1.0 if expected == "repel" else -1.0

        rows.append({
            "case": label,
            "s1": int(s1),
            "s2": int(s2),
            "separation": float(d),
            "field_energy": Efield,
            "interaction_energy": Eint,
            "gauss_residual_rms": gres,
            "gauss_source_rms": gsrc,
            "gauss_residual_rel": gres/gsrc,
            "Fx_left": Fx1,
            "Fx_right": Fx2,
            "separation_force": sep_force,
            "expected": expected,
            "expected_sign": expected_sign,
            "sign_correct": bool(np.sign(sep_force) == np.sign(expected_sign)),
            "total_signed_charge": float(np.sum(rho_q)*dx*dx),
            "abs_charge_integral": float(np.sum(np.abs(rho_q))*dx*dx),
            "phi_max_abs": float(np.max(np.abs(phi))),
            "phi_rms": float(np.sqrt(np.mean(phi**2)))
        })
    return pd.DataFrame(rows)

like_df = run_pair("like_plus_plus", +1, +1)
opp_df = run_pair("opposite_plus_minus", +1, -1)
all_df = pd.concat([like_df, opp_df], ignore_index=True)

# Derivative of energy wrt separation to check force consistency.
deriv_rows = []
for case, df in all_df.groupby("case"):
    E = df["interaction_energy"].to_numpy()
    d = df["separation"].to_numpy()
    dEdd = np.gradient(E, d)
    for sep, val in zip(d, dEdd):
        deriv_rows.append({"case": case, "separation": sep, "dEint_dsep": val})
deriv_df = pd.DataFrame(deriv_rows)
all_df = all_df.merge(deriv_df, on=["case", "separation"], how="left")

# Summary metrics
like = all_df[all_df["case"] == "like_plus_plus"]
opp = all_df[all_df["case"] == "opposite_plus_minus"]

# We only evaluate force sign away from edge and overlap extremes
mid_like = like[(like["separation"] >= 4.0) & (like["separation"] <= 10.0)]
mid_opp = opp[(opp["separation"] >= 4.0) & (opp["separation"] <= 10.0)]

like_force_correct_frac = float(mid_like["sign_correct"].mean())
opp_force_correct_frac = float(mid_opp["sign_correct"].mean())

# For like charges, interaction energy should decrease with separation in screened electrostatic box.
# For opposite charges, interaction energy should increase toward 0 with separation; its derivative should be positive.
like_dE_sign_frac = float((mid_like["dEint_dsep"] < 0).mean())
opp_dE_sign_frac = float((mid_opp["dEint_dsep"] > 0).mean())

gauss_max = float(all_df["gauss_residual_rel"].max())
finite = bool(np.isfinite(all_df.select_dtypes(include=[float, int]).to_numpy()).all())

# attraction/repulsion strength ratio at representative separation ~7
rep_sep = 7.0
like_rep = like.iloc[(like["separation"]-rep_sep).abs().argmin()]
opp_rep = opp.iloc[(opp["separation"]-rep_sep).abs().argmin()]

pass_conditions = {
    "gauss_residual_under_1e-10": bool(gauss_max < 1e-10),
    "like_force_sign_correct_all_midpoints": bool(like_force_correct_frac == 1.0),
    "opposite_force_sign_correct_all_midpoints": bool(opp_force_correct_frac == 1.0),
    "like_energy_decreases_with_separation_midpoints": bool(like_dE_sign_frac >= 0.90),
    "opposite_energy_increases_with_separation_midpoints": bool(opp_dE_sign_frac >= 0.90),
    "finite": finite,
    "like_total_charge_positive": bool(like["total_signed_charge"].mean() > 0),
    "opposite_total_charge_neutral": bool(abs(opp["total_signed_charge"].mean()) < 1e-10)
}
verdict = "PASSED as two-lump electrostatic sign sandbox" if all(pass_conditions.values()) else "PARTIAL / CHECK"

summary = {
    "G129": "Two-lump charge-sign interaction",
    "verdict": verdict,
    "scope_limit": "2D static constrained Gauss-law sandbox; not full electromagnetism/QED/Standard Model.",
    "model": {
        "charge_density": "rho_q = e(s1|psi1|^2 + s2|psi2|^2)",
        "gauss_law": "(-nabla^2 + mu^2) phi = rho_q - mean(rho_q)",
        "force": "F_i = -s_i e integral |psi_i|^2 grad(phi) dA"
    },
    "pass_conditions": pass_conditions,
    "metrics": {
        "gauss_residual_rel_max": gauss_max,
        "like_force_correct_fraction_mid": like_force_correct_frac,
        "opposite_force_correct_fraction_mid": opp_force_correct_frac,
        "like_energy_derivative_correct_fraction_mid": like_dE_sign_frac,
        "opposite_energy_derivative_correct_fraction_mid": opp_dE_sign_frac,
        "representative_separation": float(rep_sep),
        "like_sep_force_at_rep": float(like_rep["separation_force"]),
        "opposite_sep_force_at_rep": float(opp_rep["separation_force"]),
        "like_interaction_energy_at_rep": float(like_rep["interaction_energy"]),
        "opposite_interaction_energy_at_rep": float(opp_rep["interaction_energy"])
    },
    "next_gate": "G130 should test moving two-lump dynamics under the constrained field, tracking acceleration/separation over time."
}

all_df.to_csv(OUT/"G129_two_lump_scan.csv", index=False)
with open(OUT/"G129_constants_and_parameters.json", "w", encoding="utf-8") as f:
    json.dump(params, f, indent=2)
with open(OUT/"G129_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

summary_md = f"""# G129 — Two-lump charge-sign interaction

## Verdict

**{verdict}**

## Scope

This is a reduced 2D static constrained Gauss-law sandbox. It does **not** prove
electromagnetism, QED, photons, or the Standard Model. It tests the first
electrostatic sign behavior after G128B.

## Model

```text
ρ_q = e(s₁|ψ₁|² + s₂|ψ₂|²)

(−∇² + μ²)φ = ρ_q − mean(ρ_q)

Fᵢ = −sᵢ e ∫ |ψᵢ|² ∇φ dA
```

## Parameters

```json
{json.dumps(params, indent=2)}
```

## Pass conditions

{pd.DataFrame([pass_conditions]).to_markdown(index=False)}

## Core metrics

```json
{json.dumps(summary["metrics"], indent=2)}
```

## Interpretation

G129 tests whether the charge-sign seed acts electrostatically:

```text
like signs      → separation force positive  → repel
opposite signs  → separation force negative  → attract
Gauss constraint → satisfied across scan
```

Result:

```text
like force sign correct fraction      = {like_force_correct_frac:.3f}
opposite force sign correct fraction  = {opp_force_correct_frac:.3f}
max Gauss residual                    = {gauss_max:.3e}
```

At representative separation {rep_sep:.1f}:

```text
like separation force      = {float(like_rep["separation_force"]):.6e}
opposite separation force  = {float(opp_rep["separation_force"]):.6e}

like interaction energy      = {float(like_rep["interaction_energy"]):.6e}
opposite interaction energy  = {float(opp_rep["interaction_energy"]):.6e}
```

## Meaning for the TOE ladder

This supports the next step:

```text
SU/A density lump      → charge-density source
charge sign s=±1       → field orientation
Gauss-law potential φ   → mediates sign-sensitive interaction
like signs             → repel
opposite signs         → attract
```

This is still not full electric charge. It is a reduced electrostatic charge-sign
sandbox. The next real test is dynamical motion.

## Next gate

**G130 — moving two-lump dynamics**

Track whether two initially separated lumps actually move as predicted:

```text
like signs: separation increases over time
opposite signs: separation decreases over time
norm/charge controlled
Gauss residual remains low
lumps remain localized
```
"""
with open(OUT/"G129_summary.md", "w", encoding="utf-8") as f:
    f.write(summary_md)

# Plots
plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["separation"], df["separation_force"], marker="o", label=case)
plt.axhline(0, linestyle="--", linewidth=1)
plt.xlabel("separation")
plt.ylabel("separation force Fx_right - Fx_left")
plt.title("G129 sign-sensitive separation force")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G129_separation_force.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["separation"], df["interaction_energy"], marker="o", label=case)
plt.axhline(0, linestyle="--", linewidth=1)
plt.xlabel("separation")
plt.ylabel("interaction energy")
plt.title("G129 interaction energy vs separation")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G129_interaction_energy.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["separation"], df["gauss_residual_rel"], marker="o", label=case)
plt.yscale("log")
plt.xlabel("separation")
plt.ylabel("relative Gauss residual")
plt.title("G129 Gauss-law residual")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G129_gauss_residual.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["separation"], df["field_energy"], marker="o", label=case)
plt.xlabel("separation")
plt.ylabel("field energy")
plt.title("G129 field energy vs separation")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G129_field_energy.png", dpi=160)
plt.close()

readme = """# G129 Two-Lump Charge-Sign Interaction Package

Run:

```bash
python scripts/G129_two_lump_charge_sign_interaction.py
```

Outputs are written to `results/`.

This package tests whether the constrained Gauss-law seed from G128B produces
the first electrostatic sign behavior: like signs repel, opposite signs attract.
It is a reduced 2D static sandbox, not full electromagnetism or QED.
"""
with open(ROOT/"README_G129.md", "w", encoding="utf-8") as f:
    f.write(readme)

print(json.dumps(summary, indent=2))
