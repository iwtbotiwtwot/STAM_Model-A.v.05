#!/usr/bin/env python3
"""
G130 — Moving two-lump dynamics

Purpose
-------
Follow-up to G129. G129 showed static sign-sensitive force:
    like signs repel
    opposite signs attract

G130 lets two localized SU/A density lumps evolve under a constrained
Gauss-law scalar potential and tracks whether their centers actually move
as predicted.

Scope
-----
Reduced 2D electrostatic sandbox. This is NOT full electromagnetism, QED,
photons, or the Standard Model.

Important modeling choice
-------------------------
Full NLSE two-soliton collisions can be dominated by interference and internal
breathing. For this test, we use a reduced collective-coordinate dynamics model:

    each lump is represented by a Gaussian density packet with center x_i(t)
    charge sign s_i = ±1
    Gauss law solves phi from the two charge densities each time step
    force on each center is computed from -s_i e ∫ rho_i grad(phi) dA
    centers evolve by damped Newtonian dynamics

This isolates the electrostatic sign behavior from packet-shape numerics.
The next harder test would be a full wavefunction evolution with separated
component labels and re-localization tracking.

Model
-----
rho_i(x,y) = Q exp(-((x-x_i)^2 + (y-y_i)^2)/(2 sigma^2))

rho_q = e (s1 rho_1 + s2 rho_2)
(-nabla^2 + mu^2) phi = rho_q - mean(rho_q)

F_i = -s_i e ∫ rho_i grad(phi) dA

Dynamics:
M_eff x_i'' + gamma x_i' = F_i

Pass condition
--------------
- like signs: separation increases monotonically or net positive final change
- opposite signs: separation decreases monotonically or net negative final change
- Gauss residual small
- total signed charge controlled
- forces remain equal/opposite to tolerance
- no boundary collision
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
L = 34.0
dx = L / N
x = (np.arange(N) - N//2) * dx
y = (np.arange(N) - N//2) * dx
X, Y = np.meshgrid(x, y, indexing="ij")

# Fourier
kx = 2*np.pi*np.fft.fftfreq(N, d=dx)
ky = 2*np.pi*np.fft.fftfreq(N, d=dx)
KX, KY = np.meshgrid(kx, ky, indexing="ij")
K2 = KX**2 + KY**2

# Parameters
e_charge = 1.0
mu = 0.28
sigma = 1.25
Qamp = 1.0
M_eff = 22.0
gamma = 0.35
dt = 0.012
steps = 1800
sample_every = 5
initial_sep = 7.0
max_center_abs = 12.5

params = {
    "N": int(N),
    "L": float(L),
    "dx": float(dx),
    "e_charge": float(e_charge),
    "mu_screening": float(mu),
    "sigma": float(sigma),
    "Qamp": float(Qamp),
    "M_eff": float(M_eff),
    "gamma": float(gamma),
    "dt": float(dt),
    "steps": int(steps),
    "sample_every": int(sample_every),
    "initial_sep": float(initial_sep),
    "max_center_abs_allowed": float(max_center_abs),
    "scope": "2D collective-coordinate constrained Gauss-law two-lump motion sandbox"
}

def grad(f):
    return ((np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0)) / (2*dx),
            (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1)) / (2*dx))

def lap_fft(f):
    return np.real(np.fft.ifft2(-K2*np.fft.fft2(f)))

def gaussian_density(x0, y0):
    return Qamp * np.exp(-((X-x0)**2 + (Y-y0)**2)/(2*sigma**2))

def solve_phi(rho_q):
    src = rho_q - np.mean(rho_q)
    phi_k = np.fft.fft2(src)/(K2 + mu**2)
    return np.real(np.fft.ifft2(phi_k))

def gauss_residual(phi, rho_q):
    src = rho_q - np.mean(rho_q)
    lhs = -lap_fft(phi) + mu**2*phi
    res = lhs - src
    return float(np.sqrt(np.mean(res**2))), float(np.sqrt(np.mean(src**2)) + 1e-30)

def force_on_lump(phi, rho_i, sign):
    phix, phiy = grad(phi)
    Fx = -e_charge*sign*np.sum(rho_i*phix)*dx*dx
    Fy = -e_charge*sign*np.sum(rho_i*phiy)*dx*dx
    return np.array([float(Fx), float(Fy)])

def field_energy(phi):
    phix, phiy = grad(phi)
    return float(0.5*np.sum(phix**2 + phiy**2 + mu**2*phi**2)*dx*dx)

def interaction_energy(phi, rho_q):
    src = rho_q - np.mean(rho_q)
    return float(0.5*np.sum(src*phi)*dx*dx)

def run_case(label, s1, s2):
    r1 = np.array([-initial_sep/2, 0.0], dtype=float)
    r2 = np.array([ initial_sep/2, 0.0], dtype=float)
    v1 = np.array([0.0, 0.0], dtype=float)
    v2 = np.array([0.0, 0.0], dtype=float)
    rows = []

    sep0 = float(np.linalg.norm(r2-r1))

    for step in range(steps+1):
        rho1 = gaussian_density(r1[0], r1[1])
        rho2 = gaussian_density(r2[0], r2[1])
        rho_q = e_charge*(s1*rho1 + s2*rho2)
        phi = solve_phi(rho_q)
        F1 = force_on_lump(phi, rho1, s1)
        F2 = force_on_lump(phi, rho2, s2)
        gres, gsrc = gauss_residual(phi, rho_q)

        sep_vec = r2-r1
        sep = float(np.linalg.norm(sep_vec))
        sep_rate = float(np.dot(v2-v1, sep_vec)/(sep+1e-30))
        equal_opposite_err = float(np.linalg.norm(F1+F2)/(np.linalg.norm(F1)+np.linalg.norm(F2)+1e-30))

        if step % sample_every == 0:
            rows.append({
                "case": label,
                "step": int(step),
                "t": float(step*dt),
                "s1": int(s1),
                "s2": int(s2),
                "x1": float(r1[0]),
                "y1": float(r1[1]),
                "x2": float(r2[0]),
                "y2": float(r2[1]),
                "v1x": float(v1[0]),
                "v2x": float(v2[0]),
                "separation": sep,
                "separation_change": sep-sep0,
                "separation_rate": sep_rate,
                "F1x": float(F1[0]),
                "F2x": float(F2[0]),
                "separation_force": float(F2[0]-F1[0]),
                "equal_opposite_force_error": equal_opposite_err,
                "field_energy": field_energy(phi),
                "interaction_energy": interaction_energy(phi, rho_q),
                "gauss_residual_rms": gres,
                "gauss_source_rms": gsrc,
                "gauss_residual_rel": gres/gsrc,
                "total_signed_charge": float(np.sum(rho_q)*dx*dx),
                "abs_charge_integral": float(np.sum(np.abs(rho_q))*dx*dx),
                "max_center_abs": float(max(abs(r1[0]), abs(r2[0]), abs(r1[1]), abs(r2[1]))),
                "phi_rms": float(np.sqrt(np.mean(phi**2)))
            })

        if step == steps:
            break

        # Damped velocity-Verlet-ish Euler-Cromer.
        # Damping prevents endless oscillation and lets sign of interaction show cleanly.
        a1 = (F1 - gamma*v1)/M_eff
        a2 = (F2 - gamma*v2)/M_eff
        v1 = v1 + dt*a1
        v2 = v2 + dt*a2
        r1 = r1 + dt*v1
        r2 = r2 + dt*v2

    df = pd.DataFrame(rows)
    final = df.iloc[-1]
    # monotonic fraction: most sampled separation changes in expected direction
    dsep = np.diff(df["separation"].to_numpy())
    if s1*s2 > 0:
        expected = "repel"
        direction_ok_fraction = float((dsep >= -1e-4).mean())
        final_direction_ok = bool(final["separation_change"] > 0)
    else:
        expected = "attract"
        direction_ok_fraction = float((dsep <= 1e-4).mean())
        final_direction_ok = bool(final["separation_change"] < 0)

    summary = {
        "case": label,
        "s1": int(s1),
        "s2": int(s2),
        "expected": expected,
        "initial_separation": sep0,
        "final_separation": float(final["separation"]),
        "final_separation_change": float(final["separation_change"]),
        "direction_ok_fraction": direction_ok_fraction,
        "final_direction_ok": final_direction_ok,
        "max_gauss_residual_rel": float(df["gauss_residual_rel"].max()),
        "mean_equal_opposite_force_error": float(df["equal_opposite_force_error"].mean()),
        "max_center_abs": float(df["max_center_abs"].max()),
        "field_energy_initial": float(df.iloc[0]["field_energy"]),
        "field_energy_final": float(final["field_energy"]),
        "interaction_energy_initial": float(df.iloc[0]["interaction_energy"]),
        "interaction_energy_final": float(final["interaction_energy"]),
        "total_signed_charge_mean": float(df["total_signed_charge"].mean()),
        "finite": bool(np.isfinite(df.select_dtypes(include=[float, int]).to_numpy()).all())
    }
    return df, summary

like_df, like_summary = run_case("like_plus_plus_motion", +1, +1)
opp_df, opp_summary = run_case("opposite_plus_minus_motion", +1, -1)
all_df = pd.concat([like_df, opp_df], ignore_index=True)
summaries = pd.DataFrame([like_summary, opp_summary])

pass_conditions = {
    "like_final_separation_increases": bool(like_summary["final_separation_change"] > 0),
    "opposite_final_separation_decreases": bool(opp_summary["final_separation_change"] < 0),
    "like_direction_ok_fraction_over_0p90": bool(like_summary["direction_ok_fraction"] > 0.90),
    "opposite_direction_ok_fraction_over_0p90": bool(opp_summary["direction_ok_fraction"] > 0.90),
    "gauss_residual_under_1e-10": bool(max(like_summary["max_gauss_residual_rel"], opp_summary["max_gauss_residual_rel"]) < 1e-10),
    "equal_opposite_force_error_under_1e-2": bool(max(like_summary["mean_equal_opposite_force_error"], opp_summary["mean_equal_opposite_force_error"]) < 1e-2),
    "centers_remain_inside_box": bool(max(like_summary["max_center_abs"], opp_summary["max_center_abs"]) < max_center_abs),
    "finite": bool(like_summary["finite"] and opp_summary["finite"])
}
verdict = "PASSED as moving two-lump sign-dynamics sandbox" if all(pass_conditions.values()) else "PARTIAL / CHECK"

summary = {
    "G130": "Moving two-lump dynamics",
    "verdict": verdict,
    "scope_limit": "2D collective-coordinate constrained Gauss-law sandbox; not full EM/QED/Standard Model.",
    "model": {
        "charge_density": "rho_q = e(s1 rho1 + s2 rho2)",
        "gauss_law": "(-nabla^2 + mu^2) phi = rho_q - mean(rho_q)",
        "force": "F_i = -s_i e integral rho_i grad(phi) dA",
        "motion": "M_eff x_i'' + gamma x_i' = F_i"
    },
    "pass_conditions": pass_conditions,
    "like_summary": like_summary,
    "opposite_summary": opp_summary,
    "next_gate": "G131 should test full wavepacket two-component dynamics rather than collective-coordinate motion."
}

all_df.to_csv(OUT/"G130_motion_timeseries.csv", index=False)
summaries.to_csv(OUT/"G130_case_summaries.csv", index=False)
with open(OUT/"G130_constants_and_parameters.json", "w", encoding="utf-8") as f:
    json.dump(params, f, indent=2)
with open(OUT/"G130_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

summary_md = f"""# G130 — Moving two-lump dynamics

## Verdict

**{verdict}**

## Scope

This is a reduced 2D collective-coordinate constrained Gauss-law sandbox. It does
**not** prove electromagnetism, QED, photons, or the Standard Model.

This test isolates the motion of two localized density packets under the
Gauss-law potential validated in G128B/G129.

## Model

```text
ρ_q = e(s₁ρ₁ + s₂ρ₂)

(−∇² + μ²)φ = ρ_q − mean(ρ_q)

Fᵢ = −sᵢ e ∫ ρᵢ ∇φ dA

M_eff xᵢ'' + γxᵢ' = Fᵢ
```

## Parameters

```json
{json.dumps(params, indent=2)}
```

## Case summaries

{summaries.to_markdown(index=False)}

## Pass conditions

{pd.DataFrame([pass_conditions]).to_markdown(index=False)}

## Core result

```text
Like signs:
  initial separation = {like_summary["initial_separation"]:.6f}
  final separation   = {like_summary["final_separation"]:.6f}
  change             = {like_summary["final_separation_change"]:.6f}
  direction fraction = {like_summary["direction_ok_fraction"]:.3f}

Opposite signs:
  initial separation = {opp_summary["initial_separation"]:.6f}
  final separation   = {opp_summary["final_separation"]:.6f}
  change             = {opp_summary["final_separation_change"]:.6f}
  direction fraction = {opp_summary["direction_ok_fraction"]:.3f}
```

## Interpretation

G130 supports the dynamical version of G129:

```text
like signs      → separation increases over time
opposite signs  → separation decreases over time
Gauss constraint → remains satisfied
field forces     → equal/opposite to numerical tolerance
```

This is still a reduced model. It uses collective-coordinate Gaussian packets,
not full wavefunction packet dynamics. But it confirms that the constrained
charge-sign field does not merely give the right static force; it drives the
expected motion.

## Next gate

**G131 — full two-component wavepacket dynamics**

The next harder test should evolve two labeled wavefunction components:

```text
ψ₁, ψ₂ evolved dynamically
ρ_q = e(s₁|ψ₁|² + s₂|ψ₂|²)
φ solved by Gauss law each step
track packet centroids, norm, localization, and charge conservation
```
"""
with open(OUT/"G130_summary.md", "w", encoding="utf-8") as f:
    f.write(summary_md)

# Plots
plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["t"], df["separation"], label=case)
plt.xlabel("time")
plt.ylabel("separation")
plt.title("G130 separation over time")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G130_separation_over_time.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["t"], df["separation_rate"], label=case)
plt.axhline(0, linestyle="--", linewidth=1)
plt.xlabel("time")
plt.ylabel("separation rate")
plt.title("G130 separation rate")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G130_separation_rate.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["t"], df["gauss_residual_rel"], label=case)
plt.yscale("log")
plt.xlabel("time")
plt.ylabel("relative Gauss residual")
plt.title("G130 Gauss-law residual")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G130_gauss_residual.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["t"], df["interaction_energy"], label=case)
plt.xlabel("time")
plt.ylabel("interaction energy")
plt.title("G130 interaction energy")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G130_interaction_energy.png", dpi=160)
plt.close()

readme = """# G130 Moving Two-Lump Dynamics Package

Run:

```bash
python scripts/G130_moving_two_lump_dynamics.py
```

Outputs are written to `results/`.

This package tests whether the static charge-sign behavior from G129 drives
actual motion in a reduced collective-coordinate model. It is not full
electromagnetism, QED, or Standard Model physics.
"""
with open(ROOT/"README_G130.md", "w", encoding="utf-8") as f:
    f.write(readme)

print(json.dumps(summary, indent=2))
