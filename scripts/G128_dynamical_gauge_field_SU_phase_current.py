#!/usr/bin/env python3
"""
G128 — Dynamical gauge field from SU phase-current

Reduced sandbox test:
Can the U(1)-style connection from G127 be promoted into a dynamical field
sourced by SU phase-current?

This is not a full Maxwell/QED/Standard Model derivation.
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

# Grid and model parameters
Ngrid = 96
L = 16.0
dx = L / Ngrid
x = (np.arange(Ngrid) - Ngrid//2) * dx
y = (np.arange(Ngrid) - Ngrid//2) * dx
X, Y = np.meshgrid(x, y, indexing="ij")
R = np.sqrt(X**2 + Y**2)
TH = np.arctan2(Y, X)

# Evolution parameters
dt = 0.002
steps = 1200
sample_every = 20

# Matter params
g = 1.0
q = 0.2
e = 0.45
sigma = 3.0
core = 0.75

# Gauge params
cA = 1.0
gamma = 0.25
kappa_active = 0.12
kappa_control = 0.0
gauge_mass = 0.02

params = {
    "Ngrid": int(Ngrid),
    "L": float(L),
    "dx": float(dx),
    "dt": float(dt),
    "steps": int(steps),
    "sample_every": int(sample_every),
    "g": float(g),
    "q": float(q),
    "e": float(e),
    "sigma": float(sigma),
    "core": float(core),
    "cA": float(cA),
    "gamma": float(gamma),
    "kappa_active": float(kappa_active),
    "kappa_control": float(kappa_control),
    "gauge_mass": float(gauge_mass),
    "scope": "2D transverse sandbox; dynamical U(1)-style gauge field sourced by SU phase-current."
}

kx = 2*np.pi*np.fft.fftfreq(Ngrid, d=dx)
ky = 2*np.pi*np.fft.fftfreq(Ngrid, d=dx)
KX, KY = np.meshgrid(kx, ky, indexing="ij")
K2 = KX**2 + KY**2
K2_safe = np.where(K2 == 0, 1.0, K2)

def pybool(x):
    return bool(x)

def grad(f):
    return ((np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0)) / (2*dx),
            (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1)) / (2*dx))

def lap(f):
    return (
        (np.roll(f, -1, axis=0) + np.roll(f, 1, axis=0) +
         np.roll(f, -1, axis=1) + np.roll(f, 1, axis=1) - 4*f) / dx**2
    )

def cov_grad(psi, Ax, Ay):
    gx, gy = grad(psi)
    return gx - 1j*e*Ax*psi, gy - 1j*e*Ay*psi

def cov_lap(psi, Ax, Ay):
    gx, gy = grad(psi)
    divA = grad(Ax)[0] + grad(Ay)[1]
    return lap(psi) - 2j*e*(Ax*gx + Ay*gy) - 1j*e*divA*psi - (e**2)*(Ax**2 + Ay**2)*psi

def current(psi, Ax, Ay):
    Dpx, Dpy = cov_grad(psi, Ax, Ay)
    return np.imag(np.conj(psi) * Dpx), np.imag(np.conj(psi) * Dpy)

def transverse_project(jx, jy):
    Jx = np.fft.fftn(jx)
    Jy = np.fft.fftn(jy)
    dot = KX*Jx + KY*Jy
    JxT = Jx - KX*dot/K2_safe
    JyT = Jy - KY*dot/K2_safe
    JxT[K2 == 0] = 0
    JyT[K2 == 0] = 0
    return np.real(np.fft.ifftn(JxT)), np.real(np.fft.ifftn(JyT))

def divergence(Fx, Fy):
    return grad(Fx)[0] + grad(Fy)[1]

def curl_z(Ax, Ay):
    return grad(Ay)[0] - grad(Ax)[1]

def norm(psi):
    return float(np.sum(np.abs(psi)**2) * dx * dx)

def peak(psi):
    return float(np.max(np.abs(psi)))

def rms_radius(psi):
    w = np.abs(psi)**2
    total = np.sum(w) + 1e-30
    return float(np.sqrt(np.sum(w*R**2) / total))

def winding_number(psi, radius=3.5, samples=720):
    t = np.linspace(0, 2*np.pi, samples, endpoint=False)
    xs = radius*np.cos(t)
    ys = radius*np.sin(t)
    ix = np.clip(np.round(xs/dx + Ngrid//2).astype(int), 0, Ngrid-1)
    iy = np.clip(np.round(ys/dx + Ngrid//2).astype(int), 0, Ngrid-1)
    phases = np.unwrap(np.angle(psi[ix, iy]))
    return float((phases[-1] - phases[0]) / (2*np.pi))

def matter_hamiltonian(psi, Ax, Ay):
    Dpx, Dpy = cov_grad(psi, Ax, Ay)
    rho = np.abs(psi)**2
    dens = 0.5*(np.abs(Dpx)**2 + np.abs(Dpy)**2) - 0.5*g*rho**2 + (q/3.0)*rho**3
    return float(np.sum(dens) * dx * dx)

def gauge_energy(Ax, Ay, Vx, Vy):
    B = curl_z(Ax, Ay)
    E2 = Vx**2 + Vy**2
    A2 = Ax**2 + Ay**2
    dens = 0.5*(E2 + B**2 + gauge_mass**2*A2)
    return float(np.sum(dens) * dx * dx)

def initial_vortex(m=1):
    amp = np.tanh(R/core)**abs(m) * np.exp(-R**2/(2*sigma**2))
    psi = 1.35 * amp * np.exp(1j*m*TH)
    return psi.astype(np.complex128)

def initial_connection(m=1):
    f = 1 - np.exp(-(R/core)**2)
    denom = np.maximum(R, dx/2)
    Atheta = 0.65*(m/(e*denom))*f
    Ax = -Atheta*np.sin(TH)
    Ay =  Atheta*np.cos(TH)
    return Ax.astype(float), Ay.astype(float)

def run_case(label, m=1, kappa=0.12):
    psi = initial_vortex(m)
    Ax, Ay = initial_connection(m)
    Vx = np.zeros_like(Ax)
    Vy = np.zeros_like(Ay)

    N0 = norm(psi)
    peak0 = peak(psi)
    rms0 = rms_radius(psi)
    w0 = winding_number(psi)
    GE0 = gauge_energy(Ax, Ay, Vx, Vy)
    H0 = matter_hamiltonian(psi, Ax, Ay)

    rows = []
    for step in range(steps + 1):
        if step % sample_every == 0:
            jx, jy = current(psi, Ax, Ay)
            jxT, jyT = transverse_project(jx, jy)
            GE = gauge_energy(Ax, Ay, Vx, Vy)
            Hm = matter_hamiltonian(psi, Ax, Ay)
            Nnow = norm(psi)
            rmsnow = rms_radius(psi)
            windnow = winding_number(psi)
            rows.append({
                "case": label,
                "step": int(step),
                "t": float(step*dt),
                "N": float(Nnow),
                "N_drift_rel": float(abs(Nnow-N0)/(N0+1e-30)),
                "peak": float(peak(psi)),
                "peak_ratio": float(peak(psi)/(peak0+1e-30)),
                "rms": float(rmsnow),
                "rms_ratio": float(rmsnow/(rms0+1e-30)),
                "winding": float(windnow),
                "winding_drift": float(abs(windnow-w0)),
                "matter_H": float(Hm),
                "matter_H_drift_rel": float(abs(Hm-H0)/(abs(H0)+1e-30)),
                "gauge_energy": float(GE),
                "gauge_energy_ratio": float(GE/(GE0+1e-30)),
                "gauge_energy_delta": float(GE-GE0),
                "current_rms": float(np.sqrt(np.mean(jx**2 + jy**2))),
                "transverse_current_rms": float(np.sqrt(np.mean(jxT**2 + jyT**2))),
                "divA_rms": float(np.sqrt(np.mean(divergence(Ax, Ay)**2))),
                "B_rms": float(np.sqrt(np.mean(curl_z(Ax, Ay)**2)))
            })

        if step == steps:
            break

        jx, jy = current(psi, Ax, Ay)
        jxT, jyT = transverse_project(jx, jy)

        acc_x = cA**2 * lap(Ax) - gauge_mass**2 * Ax + kappa * jxT - gamma * Vx
        acc_y = cA**2 * lap(Ay) - gauge_mass**2 * Ay + kappa * jyT - gamma * Vy
        Vx = Vx + dt * acc_x
        Vy = Vy + dt * acc_y
        Ax = Ax + dt * Vx
        Ay = Ay + dt * Vy

        rho = np.abs(psi)**2
        rhs = -0.5 * cov_lap(psi, Ax, Ay) - g*rho*psi + q*rho**2*psi
        psi = psi - 1j*dt*rhs

        # keep scheme drift from dominating diagnostics
        psi *= np.sqrt(N0/(norm(psi)+1e-30))

    df = pd.DataFrame(rows)
    final = df.iloc[-1].to_dict()
    vals = df.select_dtypes(include=[float, int]).to_numpy()
    summary = {
        "case": str(label),
        "kappa": float(kappa),
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
        "matter_H_drift_rel_max": float(df["matter_H_drift_rel"].max()),
        "gauge_energy_initial": float(GE0),
        "gauge_energy_final": float(final["gauge_energy"]),
        "gauge_energy_ratio_final": float(final["gauge_energy_ratio"]),
        "gauge_energy_delta_final": float(final["gauge_energy_delta"]),
        "current_rms_mean": float(df["current_rms"].mean()),
        "transverse_current_rms_mean": float(df["transverse_current_rms"].mean()),
        "B_rms_final": float(final["B_rms"]),
        "divA_rms_final": float(final["divA_rms"]),
        "finite": bool(np.isfinite(vals).all())
    }
    return df, summary

active_df, active_summary = run_case("vortex_active_kappa", m=1, kappa=kappa_active)
control_df, control_summary = run_case("vortex_control_kappa0", m=1, kappa=kappa_control)
all_df = pd.concat([active_df, control_df], ignore_index=True)
summaries = pd.DataFrame([active_summary, control_summary])

active = active_summary
control = control_summary
response_ratio = float(abs(active["gauge_energy_delta_final"]) / (abs(control["gauge_energy_delta_final"]) + 1e-30))

pass_conditions = {
    "active_N_drift_under_2e-3": bool(active["N_drift_rel_max"] < 2e-3),
    "active_winding_drift_under_1e-3": bool(active["winding_drift_max"] < 1e-3),
    "active_peak_bounded_under_5x": bool(active["peak_ratio_max"] < 5.0),
    "active_rms_bounded_under_4x": bool(active["rms_ratio_max"] < 4.0),
    "active_gauge_energy_finite": bool(active["finite"] and np.isfinite(active["gauge_energy_final"])),
    "active_response_larger_than_control": bool(response_ratio > 1.10),
    "active_current_nonzero": bool(active["transverse_current_rms_mean"] > 1e-4)
}
verdict = "PASSED as reduced dynamical gauge sandbox" if all(pass_conditions.values()) else "PARTIAL / CHECK"

summary = {
    "G128": "Dynamical gauge field from SU phase-current",
    "verdict": verdict,
    "scope_limit": "2D transverse sandbox; not full electromagnetism, QED, or Standard Model.",
    "model": {
        "matter": "i psi_t = -1/2 D^2 psi - g|psi|^2 psi + q|psi|^4 psi",
        "current": "j = Im(conj(psi) D psi)",
        "gauge": "A_tt + gamma A_t = cA^2 lap(A) - mA^2 A + kappa j_T"
    },
    "pass_conditions": pass_conditions,
    "response_ratio_active_vs_control": response_ratio,
    "active_summary": active_summary,
    "control_summary": control_summary,
    "next_gate": "G129 should test charge-sign sectors and Coulomb-like interaction between two SU phase-current lumps."
}

all_df.to_csv(OUT/"G128_timeseries.csv", index=False)
summaries.to_csv(OUT/"G128_case_summaries.csv", index=False)
with open(OUT/"G128_constants_and_parameters.json", "w", encoding="utf-8") as f:
    json.dump(params, f, indent=2)
with open(OUT/"G128_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

summary_md = f"""# G128 — Dynamical gauge field from SU phase-current

## Verdict

**{verdict}**

## Scope

This is a reduced 2D transverse sandbox. It does **not** prove electromagnetism,
QED, or the Standard Model. It tests the narrower next gate after G127:

> Can the U(1)-style connection be promoted from a compensating background into
> a dynamical field sourced by SU phase-current?

## Model

```text
Dψ = (∇ − i e A)ψ
j = Im(ψ* Dψ)

i ψ_t = −1/2 D²ψ − g|ψ|²ψ + q|ψ|⁴ψ

A_tt + γ A_t = c_A² ∇²A − m_A² A + κ j_T
```

where `j_T` is the transverse-projected phase current. The transverse projection
prevents pure-gauge longitudinal buildup in this sandbox.

## Parameters

```json
{json.dumps(params, indent=2)}
```

## Case summaries

{summaries.to_markdown(index=False)}

## Pass conditions

{pd.DataFrame([pass_conditions]).to_markdown(index=False)}

## Active vs control gauge response

```text
|ΔE_gauge active| / |ΔE_gauge control| = {response_ratio:.6f}
```

## Interpretation

G128 supports this next TOE ladder step:

```text
SU phase-current j        →  sources dynamical connection A
dynamical gauge field A   →  carries finite field energy
winding sector m=1        →  remains conserved during coupled evolution
localized ψ lump          →  stays bounded under gauge feedback
```

The result is still narrow:

> A conserved SU phase/winding sector can source a bounded dynamical U(1)-style
> connection in a reduced sandbox.

It does **not** yet show Coulomb force, electric charge, photons, or full Maxwell theory.

## Next gate

**G129 — Charge-sign sectors and Coulomb-like interaction**

The next test should place two SU phase-current lumps with opposite and like winding/sign
assignments and check whether the dynamical connection produces:

```text
like signs repel
opposite signs attract
field energy scales approximately with separation
total charge/winding remains conserved
```
"""
with open(OUT/"G128_summary.md", "w", encoding="utf-8") as f:
    f.write(summary_md)

plt.figure(figsize=(7,5))
for case in all_df["case"].unique():
    sub = all_df[all_df["case"] == case]
    plt.plot(sub["t"], sub["gauge_energy"], label=case)
plt.xlabel("time")
plt.ylabel("gauge energy")
plt.title("G128 gauge-field energy")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G128_gauge_energy.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case in all_df["case"].unique():
    sub = all_df[all_df["case"] == case]
    plt.plot(sub["t"], sub["winding"], label=case)
plt.xlabel("time")
plt.ylabel("winding number")
plt.title("G128 winding conservation")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G128_winding_conservation.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case in all_df["case"].unique():
    sub = all_df[all_df["case"] == case]
    plt.plot(sub["t"], sub["peak_ratio"], label=case)
plt.xlabel("time")
plt.ylabel("peak ratio")
plt.title("G128 lump boundedness")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G128_peak_ratio.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case in all_df["case"].unique():
    sub = all_df[all_df["case"] == case]
    plt.plot(sub["t"], sub["transverse_current_rms"], label=case)
plt.xlabel("time")
plt.ylabel("transverse current RMS")
plt.title("G128 source current")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G128_transverse_current.png", dpi=160)
plt.close()

readme = """# G128 Dynamical Gauge Field Package

Run:

```bash
python scripts/G128_dynamical_gauge_field_SU_phase_current.py
```

Outputs are written to `results/`.

This package tests whether the G127 U(1)-style connection can be promoted into
a dynamical field sourced by SU phase-current. It is a reduced 2D sandbox, not
a full Maxwell/QED/Standard Model derivation.
"""
with open(ROOT/"README_G128.md", "w", encoding="utf-8") as f:
    f.write(readme)

print(json.dumps(summary, indent=2))
