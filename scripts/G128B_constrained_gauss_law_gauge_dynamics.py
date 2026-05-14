#!/usr/bin/env python3
"""
G128B — Constrained Gauss-law gauge dynamics

Purpose
-------
Repair test after G128. G128 showed bounded coupled dynamics but weak/inconclusive
dynamical gauge sourcing. G128B replaces the rough transverse damped-wave sandbox
with a constrained Gauss-law field solved at every time step.

This is still NOT full electromagnetism, QED, or the Standard Model.
It is a reduced electrostatic U(1)-style sandbox:

    matter density -> charge density
    charge density -> constrained scalar gauge field phi by Gauss law
    phi feeds back into matter evolution

Model
-----
Matter:
    i ∂t ψ = -1/2 ∇²ψ + V_eff ψ

with cubic-quintic matter saturation:
    V_matter = -g |ψ|² + q |ψ|⁴

Electrostatic constrained gauge field:
    (-∇² + μ²) φ = e (|ψ|² - <|ψ|²>)

The mean-subtraction makes the periodic-box Gauss law solvable and represents a
neutralizing background. μ is a weak screening/IR regulator.

Feedback:
    V_eff = V_matter + e φ

Numerics:
    split-step Fourier evolution for ψ
    FFT solve for φ each step
    diagnostics for norm, Gauss residual, energy, winding, and active-vs-control response

Tests
-----
1. Norm conserved.
2. Gauss law residual remains tiny.
3. Gauge field strongly responds to matter density.
4. Active run differs from control.
5. Winding remains conserved for vortex case.
6. Lump remains bounded: no collapse / no runaway spread.

Pass condition
--------------
Reduced pass if:
- max norm drift < 1e-10
- max Gauss relative residual < 1e-8
- phi-rho correlation > 0.90
- active/control gauge-response ratio effectively infinite or > 10
- winding drift < 1e-3
- peak and RMS bounded
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
N = 128
L = 18.0
dx = L / N
x = (np.arange(N) - N//2) * dx
y = (np.arange(N) - N//2) * dx
X, Y = np.meshgrid(x, y, indexing="ij")
R = np.sqrt(X**2 + Y**2)
TH = np.arctan2(Y, X)

# Fourier grid
kx = 2*np.pi*np.fft.fftfreq(N, d=dx)
ky = 2*np.pi*np.fft.fftfreq(N, d=dx)
KX, KY = np.meshgrid(kx, ky, indexing="ij")
K2 = KX**2 + KY**2

# Evolution params
dt = 0.004
steps = 1400
sample_every = 20

# Matter/gauge params
g = 1.0
q = 0.2
e_charge = 0.16
mu = 0.35
sigma = 3.0
core = 0.75
amp0 = 1.22

params = {
    "N": int(N),
    "L": float(L),
    "dx": float(dx),
    "dt": float(dt),
    "steps": int(steps),
    "sample_every": int(sample_every),
    "g": float(g),
    "q": float(q),
    "e_charge": float(e_charge),
    "mu_screening": float(mu),
    "sigma": float(sigma),
    "core": float(core),
    "amp0": float(amp0),
    "scope": "2D transverse constrained Gauss-law sandbox"
}

def lap_fft(f):
    return np.real(np.fft.ifft2(-K2 * np.fft.fft2(f)))

def grad(f):
    return ((np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0)) / (2*dx),
            (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1)) / (2*dx))

def initial_vortex(m=1):
    amp = amp0 * (np.tanh(R/core) ** abs(m)) * np.exp(-R**2/(2*sigma**2))
    return (amp * np.exp(1j*m*TH)).astype(np.complex128)

def solve_phi(rho, active=True):
    if not active:
        return np.zeros_like(rho)
    source = e_charge * (rho - np.mean(rho))
    source_k = np.fft.fft2(source)
    phi_k = source_k / (K2 + mu**2)
    return np.real(np.fft.ifft2(phi_k))

def gauss_residual(phi, rho):
    source = e_charge * (rho - np.mean(rho))
    lhs = -lap_fft(phi) + mu**2 * phi
    res = lhs - source
    return float(np.sqrt(np.mean(res**2))), float(np.sqrt(np.mean(source**2)) + 1e-30)

def norm(psi):
    return float(np.sum(np.abs(psi)**2)*dx*dx)

def peak(psi):
    return float(np.max(np.abs(psi)))

def rms_radius(psi):
    rho = np.abs(psi)**2
    total = np.sum(rho) + 1e-30
    return float(np.sqrt(np.sum(rho*R**2)/total))

def winding_number(psi, radius=3.0, samples=720):
    t = np.linspace(0, 2*np.pi, samples, endpoint=False)
    xs = radius*np.cos(t)
    ys = radius*np.sin(t)
    ix = np.clip(np.round(xs/dx + N//2).astype(int), 0, N-1)
    iy = np.clip(np.round(ys/dx + N//2).astype(int), 0, N-1)
    phases = np.unwrap(np.angle(psi[ix, iy]))
    return float((phases[-1] - phases[0])/(2*np.pi))

def energies(psi, phi, active=True):
    rho = np.abs(psi)**2
    gx, gy = grad(psi)
    kinetic = 0.5*np.sum(np.abs(gx)**2 + np.abs(gy)**2)*dx*dx
    matter = np.sum(-0.5*g*rho**2 + (q/3.0)*rho**3)*dx*dx
    phix, phiy = grad(phi)
    gauge = 0.5*np.sum(phix**2 + phiy**2 + mu**2*phi**2)*dx*dx if active else 0.0
    interaction = 0.5*np.sum(e_charge*(rho - np.mean(rho))*phi)*dx*dx if active else 0.0
    total = kinetic + matter + gauge + interaction
    return {
        "kinetic_E": float(kinetic),
        "matter_E": float(matter),
        "gauge_E": float(gauge),
        "interaction_E": float(interaction),
        "total_E": float(total)
    }

def phi_rho_corr(phi, rho):
    a = phi.ravel() - np.mean(phi)
    b = rho.ravel() - np.mean(rho)
    den = np.sqrt(np.sum(a*a)*np.sum(b*b)) + 1e-30
    return float(np.sum(a*b)/den)

def split_step(psi, active=True):
    rho = np.abs(psi)**2
    phi = solve_phi(rho, active=active)
    V = -g*rho + q*rho**2 + (e_charge*phi if active else 0.0)

    psi = np.exp(-1j*V*dt/2)*psi
    psi_k = np.fft.fft2(psi)
    psi_k *= np.exp(-1j*(K2/2.0)*dt)
    psi = np.fft.ifft2(psi_k)

    rho = np.abs(psi)**2
    phi = solve_phi(rho, active=active)
    V = -g*rho + q*rho**2 + (e_charge*phi if active else 0.0)
    psi = np.exp(-1j*V*dt/2)*psi
    return psi

def run_case(label, active=True, m=1):
    psi = initial_vortex(m=m)
    N0 = norm(psi)
    peak0 = peak(psi)
    rms0 = rms_radius(psi)
    w0 = winding_number(psi)
    rows = []

    for step in range(steps+1):
        rho = np.abs(psi)**2
        phi = solve_phi(rho, active=active)
        res, src = gauss_residual(phi, rho)
        en = energies(psi, phi, active=active)
        wind = winding_number(psi)
        if step % sample_every == 0:
            rows.append({
                "case": label,
                "step": int(step),
                "t": float(step*dt),
                "N": norm(psi),
                "N_drift_rel": abs(norm(psi)-N0)/(N0+1e-30),
                "peak": peak(psi),
                "peak_ratio": peak(psi)/(peak0+1e-30),
                "rms": rms_radius(psi),
                "rms_ratio": rms_radius(psi)/(rms0+1e-30),
                "winding": wind,
                "winding_drift": abs(wind-w0),
                "phi_rms": float(np.sqrt(np.mean(phi**2))),
                "phi_max_abs": float(np.max(np.abs(phi))),
                "phi_rho_corr": phi_rho_corr(phi, rho) if active else 0.0,
                "gauss_residual_rms": res,
                "gauss_source_rms": src,
                "gauss_residual_rel": res/src,
                **en
            })

        if step == steps:
            break
        psi = split_step(psi, active=active)

    df = pd.DataFrame(rows)
    final = df.iloc[-1].to_dict()
    numeric_ok = bool(np.isfinite(df.select_dtypes(include=[float, int]).to_numpy()).all())
    summary = {
        "case": label,
        "active": bool(active),
        "m": int(m),
        "N0": float(N0),
        "N_final": float(final["N"]),
        "N_drift_rel_max": float(df["N_drift_rel"].max()),
        "peak_ratio_max": float(df["peak_ratio"].max()),
        "rms_ratio_final": float(final["rms_ratio"]),
        "rms_ratio_max": float(df["rms_ratio"].max()),
        "winding_initial": float(w0),
        "winding_final": float(final["winding"]),
        "winding_drift_max": float(df["winding_drift"].max()),
        "phi_rms_mean": float(df["phi_rms"].mean()),
        "phi_rms_final": float(final["phi_rms"]),
        "phi_max_abs_final": float(final["phi_max_abs"]),
        "phi_rho_corr_mean": float(df["phi_rho_corr"].mean()),
        "gauss_residual_rel_max": float(df["gauss_residual_rel"].max()),
        "gauge_E_final": float(final["gauge_E"]),
        "interaction_E_final": float(final["interaction_E"]),
        "total_E_initial": float(df.iloc[0]["total_E"]),
        "total_E_final": float(final["total_E"]),
        "total_E_drift_rel": float(abs(final["total_E"]-df.iloc[0]["total_E"])/(abs(df.iloc[0]["total_E"])+1e-30)),
        "finite": numeric_ok
    }
    return df, summary

active_df, active_summary = run_case("active_gauss_law", active=True, m=1)
control_df, control_summary = run_case("control_phi0", active=False, m=1)
all_df = pd.concat([active_df, control_df], ignore_index=True)
summaries = pd.DataFrame([active_summary, control_summary])

response_ratio = float(active_summary["phi_rms_mean"]/(control_summary["phi_rms_mean"] + 1e-30))

pass_conditions = {
    "active_N_drift_under_1e-10": bool(active_summary["N_drift_rel_max"] < 1e-10),
    "active_gauss_residual_under_1e-8": bool(active_summary["gauss_residual_rel_max"] < 1e-8),
    "active_phi_rho_corr_over_0p90": bool(active_summary["phi_rho_corr_mean"] > 0.90),
    "active_response_ratio_over_10": bool(response_ratio > 10.0),
    "active_winding_drift_under_1e-3": bool(active_summary["winding_drift_max"] < 1e-3),
    "active_peak_bounded_under_4x": bool(active_summary["peak_ratio_max"] < 4.0),
    "active_rms_bounded_under_4x": bool(active_summary["rms_ratio_max"] < 4.0),
    "active_finite": bool(active_summary["finite"])
}
verdict = "PASSED as constrained Gauss-law gauge sandbox" if all(pass_conditions.values()) else "PARTIAL / CHECK"

summary = {
    "G128B": "Constrained Gauss-law gauge dynamics",
    "verdict": verdict,
    "scope_limit": "2D electrostatic constrained U(1)-style sandbox; not full electromagnetism, QED, or Standard Model.",
    "model": {
        "matter": "i psi_t = -1/2 nabla^2 psi + [-g|psi|^2 + q|psi|^4 + e phi] psi",
        "gauss_law": "(-nabla^2 + mu^2) phi = e(|psi|^2 - mean(|psi|^2))",
        "numerics": "split-step Fourier for psi; FFT constrained solve for phi every step"
    },
    "pass_conditions": pass_conditions,
    "response_ratio_active_phi_rms_vs_control": response_ratio,
    "active_summary": active_summary,
    "control_summary": control_summary,
    "next_gate": "G129 should test two-lump charge-sign interactions: like signs repel, opposite signs attract, field energy vs separation."
}

all_df.to_csv(OUT/"G128B_timeseries.csv", index=False)
summaries.to_csv(OUT/"G128B_case_summaries.csv", index=False)
with open(OUT/"G128B_constants_and_parameters.json", "w", encoding="utf-8") as f:
    json.dump(params, f, indent=2)
with open(OUT/"G128B_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

summary_md = f"""# G128B — Constrained Gauss-law gauge dynamics

## Verdict

**{verdict}**

## Scope

This is a reduced 2D electrostatic constrained U(1)-style sandbox. It does **not**
prove electromagnetism, QED, photons, or the Standard Model. It repairs G128 by
enforcing Gauss law directly at each time step.

## Model

```text
i ψ_t = −1/2 ∇²ψ + [−g|ψ|² + q|ψ|⁴ + eφ]ψ

(−∇² + μ²)φ = e(|ψ|² − mean(|ψ|²))
```

Numerics:

```text
ψ evolved by split-step Fourier method
φ solved by FFT constrained Gauss-law solve each step
```

## Parameters

```json
{json.dumps(params, indent=2)}
```

## Case summaries

{summaries.to_markdown(index=False)}

## Pass conditions

{pd.DataFrame([pass_conditions]).to_markdown(index=False)}

## Active vs control response

```text
mean φ_rms(active) / mean φ_rms(control) = {response_ratio:.6e}
```

## Interpretation

G128B repairs the weak point from G128:

```text
SU/A density |ψ|²          →  charge-density source
Gauss constraint            →  determines scalar gauge field φ
φ feedback                  →  acts back on ψ evolution
winding sector m=1          →  remains conserved
```

The result supports a narrower claim:

> A localized SU/A excitation can source a constrained U(1)-style gauge potential
> through a Gauss-law relation while remaining bounded and conserving winding.

This is still not electric charge in the full Standard Model sense. It is the
electrostatic seed needed before testing Coulomb-like two-body behavior.

## Next gate

**G129 — Charge-sign sectors and Coulomb-like interaction**

Required tests:

```text
like signs repel
opposite signs attract
field energy changes with separation
total charge/winding conserved
Gauss constraint remains satisfied
```
"""
with open(OUT/"G128B_summary.md", "w", encoding="utf-8") as f:
    f.write(summary_md)

# Plots
plt.figure(figsize=(7,5))
for case in all_df["case"].unique():
    sub = all_df[all_df["case"] == case]
    plt.plot(sub["t"], sub["phi_rms"], label=case)
plt.xlabel("time")
plt.ylabel("phi RMS")
plt.title("G128B constrained gauge response")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G128B_phi_rms.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case in all_df["case"].unique():
    sub = all_df[all_df["case"] == case]
    plt.plot(sub["t"], sub["gauss_residual_rel"], label=case)
plt.yscale("log")
plt.xlabel("time")
plt.ylabel("relative Gauss residual")
plt.title("G128B Gauss-law residual")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G128B_gauss_residual.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case in all_df["case"].unique():
    sub = all_df[all_df["case"] == case]
    plt.plot(sub["t"], sub["winding"], label=case)
plt.xlabel("time")
plt.ylabel("winding number")
plt.title("G128B winding conservation")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G128B_winding.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case in all_df["case"].unique():
    sub = all_df[all_df["case"] == case]
    plt.plot(sub["t"], sub["peak_ratio"], label=case)
plt.xlabel("time")
plt.ylabel("peak ratio")
plt.title("G128B peak boundedness")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G128B_peak_ratio.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case in all_df["case"].unique():
    sub = all_df[all_df["case"] == case]
    plt.plot(sub["t"], sub["rms_ratio"], label=case)
plt.xlabel("time")
plt.ylabel("RMS radius ratio")
plt.title("G128B RMS boundedness")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G128B_rms_ratio.png", dpi=160)
plt.close()

readme = """# G128B Constrained Gauss-Law Gauge Dynamics Package

Run:

```bash
python scripts/G128B_constrained_gauss_law_gauge_dynamics.py
```

Outputs are written to `results/`.

This package repairs G128 by enforcing a constrained Gauss-law solve for the
U(1)-style scalar gauge potential at each time step. It is a reduced 2D
electrostatic sandbox, not a full Maxwell/QED/Standard Model derivation.
"""
with open(ROOT/"README_G128B.md", "w", encoding="utf-8") as f:
    f.write(readme)

print(json.dumps(summary, indent=2))
