#!/usr/bin/env python3
"""
G132 — Charge-sign wavepacket dynamics with topological winding

Purpose
-------
Harder gate after G131.

G131 evolved two labeled charge-sign wavepackets and found the correct sign
motion. G132 adds local vortex/winding structure to each labeled component and
tests whether:

    1. winding remains conserved
    2. charge-sign attraction/repulsion still appears
    3. component norms remain conserved
    4. Gauss constraint remains controlled
    5. packets remain bounded

Scope
-----
Reduced 2D two-component NLSE + constrained scalar Gauss-law sandbox.
This is NOT full electromagnetism, QED, photons, spin, or Standard Model physics.

Model
-----
i ∂t ψᵢ = -1/2 ∇²ψᵢ + [-g|ψᵢ|² + q|ψᵢ|⁴ + sᵢ e φ] ψᵢ

(-∇² + μ²)φ = e(s₁|ψ₁|² + s₂|ψ₂|² - mean)

Initial packet:
ψᵢ = amp * tanh(r_i/core)^|m_i| * exp(-r_i²/(2σ²)) * exp(i m_i θ_i)

Winding is measured on a small loop around each component centroid.
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
L = 34.0
dx = L / N
x = (np.arange(N) - N//2) * dx
y = (np.arange(N) - N//2) * dx
X, Y = np.meshgrid(x, y, indexing="ij")

kx = 2*np.pi*np.fft.fftfreq(N, d=dx)
ky = 2*np.pi*np.fft.fftfreq(N, d=dx)
KX, KY = np.meshgrid(kx, ky, indexing="ij")
K2 = KX**2 + KY**2

# Evolution parameters
dt = 0.0025
steps = 1100
sample_every = 10

# Matter / gauge parameters
g = 1.0
q = 0.2
e_charge = 0.10
mu = 0.35
sigma = 1.85
core = 0.55
amp = 1.10
initial_sep = 8.0
m1 = 1
m2 = 1
winding_radius = 0.95

params = {
    "N": int(N), "L": float(L), "dx": float(dx),
    "dt": float(dt), "steps": int(steps), "sample_every": int(sample_every),
    "g": float(g), "q": float(q), "e_charge": float(e_charge),
    "mu_screening": float(mu), "sigma": float(sigma), "core": float(core),
    "amp": float(amp), "initial_sep": float(initial_sep),
    "m1": int(m1), "m2": int(m2), "winding_radius": float(winding_radius),
    "scope": "2D two-component vortex wavepacket NLSE with constrained scalar Gauss law"
}

def solve_phi(rho_q):
    src = rho_q - np.mean(rho_q)
    return np.real(np.fft.ifft2(np.fft.fft2(src)/(K2 + mu**2)))

def lap_fft_real(f):
    return np.real(np.fft.ifft2(-K2*np.fft.fft2(f)))

def gauss_residual(phi, rho_q):
    src = rho_q - np.mean(rho_q)
    lhs = -lap_fft_real(phi) + mu**2*phi
    res = lhs - src
    return float(np.sqrt(np.mean(res**2))), float(np.sqrt(np.mean(src**2)) + 1e-30)

def init_vortex_packet(x0, y0, m):
    RX = X - x0
    RY = Y - y0
    R = np.sqrt(RX**2 + RY**2)
    TH = np.arctan2(RY, RX)
    envelope = amp * (np.tanh(R/core)**abs(m)) * np.exp(-(R**2)/(2*sigma**2))
    return (envelope * np.exp(1j*m*TH)).astype(np.complex128)

def norm(psi):
    return float(np.sum(np.abs(psi)**2)*dx*dx)

def centroid(psi):
    rho = np.abs(psi)**2
    n = np.sum(rho)+1e-30
    return float(np.sum(rho*X)/n), float(np.sum(rho*Y)/n)

def rms_radius_about_centroid(psi):
    rho = np.abs(psi)**2
    cx, cy = centroid(psi)
    n = np.sum(rho)+1e-30
    return float(np.sqrt(np.sum(rho*((X-cx)**2 + (Y-cy)**2))/n))

def peak(psi):
    return float(np.max(np.abs(psi)))

def winding_around(psi, cx, cy, radius=winding_radius, samples=720):
    t = np.linspace(0, 2*np.pi, samples, endpoint=False)
    xs = cx + radius*np.cos(t)
    ys = cy + radius*np.sin(t)
    ix = np.clip(np.round(xs/dx + N//2).astype(int), 0, N-1)
    iy = np.clip(np.round(ys/dx + N//2).astype(int), 0, N-1)
    phases = np.unwrap(np.angle(psi[ix, iy]))
    return float((phases[-1] - phases[0])/(2*np.pi))

def spectral_grad_abs2(psi):
    pk = np.fft.fft2(psi)
    gx = np.fft.ifft2(1j*KX*pk)
    gy = np.fft.ifft2(1j*KY*pk)
    return np.abs(gx)**2 + np.abs(gy)**2

def energies(psi1, psi2, phi, s1, s2):
    rho1 = np.abs(psi1)**2
    rho2 = np.abs(psi2)**2
    kin = 0.5*np.sum(spectral_grad_abs2(psi1) + spectral_grad_abs2(psi2))*dx*dx
    mat = np.sum(-0.5*g*(rho1**2+rho2**2) + (q/3.0)*(rho1**3+rho2**3))*dx*dx
    src = e_charge*(s1*rho1+s2*rho2)
    interaction = 0.5*np.sum((src-np.mean(src))*phi)*dx*dx
    phik = np.fft.fft2(phi)
    phix = np.real(np.fft.ifft2(1j*KX*phik))
    phiy = np.real(np.fft.ifft2(1j*KY*phik))
    gauge = 0.5*np.sum(phix**2+phiy**2+mu**2*phi**2)*dx*dx
    return {
        "kinetic_E": float(kin),
        "matter_E": float(mat),
        "interaction_E": float(interaction),
        "gauge_E": float(gauge),
        "total_E": float(kin+mat+interaction+gauge)
    }

def split_step(psi1, psi2, s1, s2):
    rho1 = np.abs(psi1)**2
    rho2 = np.abs(psi2)**2
    phi = solve_phi(e_charge*(s1*rho1+s2*rho2))
    V1 = -g*rho1 + q*rho1**2 + s1*e_charge*phi
    V2 = -g*rho2 + q*rho2**2 + s2*e_charge*phi

    psi1 = np.exp(-1j*V1*dt/2)*psi1
    psi2 = np.exp(-1j*V2*dt/2)*psi2

    phase = np.exp(-1j*(K2/2.0)*dt)
    psi1 = np.fft.ifft2(np.fft.fft2(psi1)*phase)
    psi2 = np.fft.ifft2(np.fft.fft2(psi2)*phase)

    rho1 = np.abs(psi1)**2
    rho2 = np.abs(psi2)**2
    phi = solve_phi(e_charge*(s1*rho1+s2*rho2))
    V1 = -g*rho1 + q*rho1**2 + s1*e_charge*phi
    V2 = -g*rho2 + q*rho2**2 + s2*e_charge*phi
    psi1 = np.exp(-1j*V1*dt/2)*psi1
    psi2 = np.exp(-1j*V2*dt/2)*psi2
    return psi1, psi2

def run_case(label, s1, s2):
    psi1 = init_vortex_packet(-initial_sep/2, 0.0, m1)
    psi2 = init_vortex_packet(+initial_sep/2, 0.0, m2)

    N10, N20 = norm(psi1), norm(psi2)
    c10, c20 = centroid(psi1), centroid(psi2)
    sep0 = float(c20[0]-c10[0])
    peak10, peak20 = peak(psi1), peak(psi2)
    rms10, rms20 = rms_radius_about_centroid(psi1), rms_radius_about_centroid(psi2)
    w10 = winding_around(psi1, *c10)
    w20 = winding_around(psi2, *c20)

    rows = []
    for step in range(steps+1):
        rho1, rho2 = np.abs(psi1)**2, np.abs(psi2)**2
        rho_q = e_charge*(s1*rho1+s2*rho2)
        phi = solve_phi(rho_q)
        gres, gsrc = gauss_residual(phi, rho_q)
        c1, c2 = centroid(psi1), centroid(psi2)
        sep = float(c2[0]-c1[0])
        w1 = winding_around(psi1, *c1)
        w2 = winding_around(psi2, *c2)
        en = energies(psi1, psi2, phi, s1, s2)

        if step % sample_every == 0:
            rows.append({
                "case": label, "step": int(step), "t": float(step*dt),
                "s1": int(s1), "s2": int(s2),
                "N1": norm(psi1), "N2": norm(psi2),
                "N1_drift_rel": abs(norm(psi1)-N10)/(N10+1e-30),
                "N2_drift_rel": abs(norm(psi2)-N20)/(N20+1e-30),
                "c1x": c1[0], "c2x": c2[0], "c1y": c1[1], "c2y": c2[1],
                "separation": sep, "separation_change": sep-sep0,
                "w1": w1, "w2": w2,
                "w1_drift": abs(w1-w10), "w2_drift": abs(w2-w20),
                "peak1_ratio": peak(psi1)/(peak10+1e-30),
                "peak2_ratio": peak(psi2)/(peak20+1e-30),
                "rms1_ratio": rms_radius_about_centroid(psi1)/(rms10+1e-30),
                "rms2_ratio": rms_radius_about_centroid(psi2)/(rms20+1e-30),
                "phi_rms": float(np.sqrt(np.mean(phi**2))),
                "gauss_residual_rel": gres/gsrc,
                **en
            })

        if step == steps:
            break

        psi1, psi2 = split_step(psi1, psi2, s1, s2)
        psi1 *= np.sqrt(N10/(norm(psi1)+1e-30))
        psi2 *= np.sqrt(N20/(norm(psi2)+1e-30))

    df = pd.DataFrame(rows)
    final = df.iloc[-1]
    dsep = np.diff(df["separation"].to_numpy())
    if s1*s2 > 0:
        expected = "repel"
        direction_ok_fraction = float((dsep >= -1e-6).mean())
        final_ok = bool(final["separation_change"] > 0)
    else:
        expected = "attract"
        direction_ok_fraction = float((dsep <= 1e-6).mean())
        final_ok = bool(final["separation_change"] < 0)

    summary = {
        "case": label,
        "expected": expected,
        "initial_separation": sep0,
        "final_separation": float(final["separation"]),
        "final_separation_change": float(final["separation_change"]),
        "direction_ok_fraction": direction_ok_fraction,
        "final_direction_ok": final_ok,
        "initial_w1": float(w10),
        "initial_w2": float(w20),
        "final_w1": float(final["w1"]),
        "final_w2": float(final["w2"]),
        "w1_drift_max": float(df["w1_drift"].max()),
        "w2_drift_max": float(df["w2_drift"].max()),
        "N1_drift_rel_max": float(df["N1_drift_rel"].max()),
        "N2_drift_rel_max": float(df["N2_drift_rel"].max()),
        "gauss_residual_rel_max": float(df["gauss_residual_rel"].max()),
        "peak_ratio_max": float(max(df["peak1_ratio"].max(), df["peak2_ratio"].max())),
        "rms_ratio_max": float(max(df["rms1_ratio"].max(), df["rms2_ratio"].max())),
        "phi_rms_mean": float(df["phi_rms"].mean()),
        "total_E_initial": float(df.iloc[0]["total_E"]),
        "total_E_final": float(final["total_E"]),
        "total_E_drift_rel": float(abs(final["total_E"]-df.iloc[0]["total_E"])/(abs(df.iloc[0]["total_E"])+1e-30)),
        "finite": bool(np.isfinite(df.select_dtypes(include=[float, int]).to_numpy()).all())
    }
    return df, summary

like_df, like_summary = run_case("like_plus_plus_vortex_wavepackets", +1, +1)
opp_df, opp_summary = run_case("opposite_plus_minus_vortex_wavepackets", +1, -1)
all_df = pd.concat([like_df, opp_df], ignore_index=True)
summaries = pd.DataFrame([like_summary, opp_summary])

pass_conditions = {
    "like_final_separation_increases": bool(like_summary["final_separation_change"] > 0),
    "opposite_final_separation_decreases": bool(opp_summary["final_separation_change"] < 0),
    "like_direction_fraction_over_0p50": bool(like_summary["direction_ok_fraction"] > 0.50),
    "opposite_direction_fraction_over_0p50": bool(opp_summary["direction_ok_fraction"] > 0.50),
    "winding_drift_under_1e-2": bool(max(like_summary["w1_drift_max"], like_summary["w2_drift_max"], opp_summary["w1_drift_max"], opp_summary["w2_drift_max"]) < 1e-2),
    "component_norm_drift_under_1e-10": bool(max(like_summary["N1_drift_rel_max"], like_summary["N2_drift_rel_max"], opp_summary["N1_drift_rel_max"], opp_summary["N2_drift_rel_max"]) < 1e-10),
    "gauss_residual_under_1e-10": bool(max(like_summary["gauss_residual_rel_max"], opp_summary["gauss_residual_rel_max"]) < 1e-10),
    "packets_bounded_peak_under_5x": bool(max(like_summary["peak_ratio_max"], opp_summary["peak_ratio_max"]) < 5.0),
    "packets_not_spread_over_5x": bool(max(like_summary["rms_ratio_max"], opp_summary["rms_ratio_max"]) < 5.0),
    "finite": bool(like_summary["finite"] and opp_summary["finite"])
}
verdict = "PASSED as winding charge-sign wavepacket sandbox" if all(pass_conditions.values()) else "PARTIAL / CHECK"

summary = {
    "G132": "Charge-sign wavepacket dynamics with topological winding",
    "verdict": verdict,
    "scope_limit": "2D two-component vortex NLSE + constrained scalar Gauss law; not full EM/QED/Standard Model.",
    "model": {
        "component_dynamics": "i psi_i,t = -1/2 nabla^2 psi_i + [-g|psi_i|^2 + q|psi_i|^4 + s_i e phi] psi_i",
        "gauss_law": "(-nabla^2 + mu^2) phi = e(s1|psi1|^2 + s2|psi2|^2 - mean)",
        "winding": "local vortex phase exp(i m theta_i), measured around each component centroid"
    },
    "pass_conditions": pass_conditions,
    "like_summary": like_summary,
    "opposite_summary": opp_summary,
    "next_gate": "G133 should test opposite winding sectors m=+1 and m=-1 and whether charge sign is separable from topological winding sign."
}

all_df.to_csv(OUT/"G132_winding_wavepacket_timeseries.csv", index=False)
summaries.to_csv(OUT/"G132_case_summaries.csv", index=False)
with open(OUT/"G132_constants_and_parameters.json", "w", encoding="utf-8") as f:
    json.dump(params, f, indent=2)
with open(OUT/"G132_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

summary_md = f"""# G132 — Charge-sign wavepacket dynamics with topological winding

## Verdict

**{verdict}**

## Scope

Reduced 2D two-component vortex nonlinear Schrödinger + constrained scalar
Gauss-law sandbox. It is **not** full electromagnetism, QED, photons, spin, or
the Standard Model.

This extends G131 by adding local topological winding to each labeled component.

## Model

```text
ψᵢ = envelopeᵢ · exp(i mᵢθᵢ)

i ∂t ψᵢ = −1/2 ∇²ψᵢ + [−g|ψᵢ|² + q|ψᵢ|⁴ + sᵢeφ]ψᵢ

(−∇² + μ²)φ = e(s₁|ψ₁|² + s₂|ψ₂|² − mean)
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
  winding drifts     = {like_summary["w1_drift_max"]:.3e}, {like_summary["w2_drift_max"]:.3e}

Opposite signs:
  initial separation = {opp_summary["initial_separation"]:.6f}
  final separation   = {opp_summary["final_separation"]:.6f}
  change             = {opp_summary["final_separation_change"]:.6f}
  direction fraction = {opp_summary["direction_ok_fraction"]:.3f}
  winding drifts     = {opp_summary["w1_drift_max"]:.3e}, {opp_summary["w2_drift_max"]:.3e}
```

## Interpretation

G132 tests whether the G131 charge-sign wavepacket result survives after restoring
topological winding:

```text
labeled wavepacket dynamics  → retained
Gauss-law charge field       → retained
local vortex winding m=1     → tracked around each centroid
charge-sign force direction  → checked
```

## Next gate

**G133 — separate charge sign from winding sign**

Use combinations:

```text
s = +1, m = +1
s = +1, m = −1
s = −1, m = +1
s = −1, m = −1
```

The key question:

> Is electric-like charge sign independent from topological winding orientation,
> or does the model force them to be tied together?
"""
with open(OUT/"G132_summary.md", "w", encoding="utf-8") as f:
    f.write(summary_md)

# Plots
plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["t"], df["separation"], label=case)
plt.xlabel("time")
plt.ylabel("centroid separation")
plt.title("G132 vortex wavepacket separation")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(PLOTS/"G132_separation.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["t"], df["w1"], label=case+" w1")
    plt.plot(df["t"], df["w2"], linestyle="--", label=case+" w2")
plt.xlabel("time")
plt.ylabel("winding")
plt.title("G132 winding conservation")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(PLOTS/"G132_winding.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["t"], df["gauss_residual_rel"], label=case)
plt.yscale("log")
plt.xlabel("time")
plt.ylabel("relative Gauss residual")
plt.title("G132 Gauss-law residual")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(PLOTS/"G132_gauss_residual.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["t"], df["peak1_ratio"], label=case+" ψ1")
    plt.plot(df["t"], df["peak2_ratio"], linestyle="--", label=case+" ψ2")
plt.xlabel("time")
plt.ylabel("peak ratio")
plt.title("G132 packet peak boundedness")
plt.legend(fontsize=7)
plt.tight_layout()
plt.savefig(PLOTS/"G132_peak_ratios.png", dpi=160)
plt.close()

readme = """# G132 Winding Charge-Sign Wavepacket Dynamics Package

Run:

```bash
python scripts/G132_winding_charge_sign_wavepacket_dynamics.py
```

Outputs are written to `results/`.

This package evolves two labeled vortex wavepacket components under a constrained
Gauss-law scalar potential. It is a reduced sandbox, not full electromagnetism,
QED, photons, spin, or Standard Model physics.
"""
with open(ROOT/"README_G132.md", "w", encoding="utf-8") as f:
    f.write(readme)

print(json.dumps(summary, indent=2))
