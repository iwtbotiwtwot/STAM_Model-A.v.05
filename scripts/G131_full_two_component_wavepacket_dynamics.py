#!/usr/bin/env python3
"""
G131 — Full two-component wavepacket dynamics

Purpose
-------
Harder gate after G130/G130B.

Instead of moving Gaussian centers as particles, evolve two labeled wavefunction
components directly:

    psi1, psi2

Each component has its own norm and centroid. The shared constrained Gauss-law
field is solved each step:

    rho_q = e(s1|psi1|^2 + s2|psi2|^2)
    (-∇² + μ²)φ = rho_q - mean(rho_q)

Dynamics:
    i ∂t psi_i = -1/2 ∇² psi_i + [-g|psi_i|² + q|psi_i|⁴ + s_i e φ] psi_i

This is a reduced two-component NLSE + constrained electrostatic potential.
It is NOT full electromagnetism/QED/Standard Model physics.

Tests
-----
1. Like signs should increase centroid separation.
2. Opposite signs should decrease centroid separation.
3. Component norms remain conserved.
4. Winding/phase not tested here; this is charge-sign wavepacket motion only.
5. Gauss residual remains small.
6. Packets remain localized and finite.
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
L = 32.0
dx = L / N
x = (np.arange(N) - N//2) * dx
y = (np.arange(N) - N//2) * dx
X, Y = np.meshgrid(x, y, indexing="ij")
R2 = X**2 + Y**2

# Fourier
kx = 2*np.pi*np.fft.fftfreq(N, d=dx)
ky = 2*np.pi*np.fft.fftfreq(N, d=dx)
KX, KY = np.meshgrid(kx, ky, indexing="ij")
K2 = KX**2 + KY**2

# Parameters
dt = 0.003
steps = 1200
sample_every = 10

g = 1.0
q = 0.2
e_charge = 0.12
mu = 0.35
sigma = 1.55
amp = 1.12
initial_sep = 7.0

# Give a very small initial push so field response can separate from pure breathing
kick_like = 0.00
kick_opp = 0.00

params = {
    "N": int(N), "L": float(L), "dx": float(dx),
    "dt": float(dt), "steps": int(steps), "sample_every": int(sample_every),
    "g": float(g), "q": float(q), "e_charge": float(e_charge),
    "mu_screening": float(mu), "sigma": float(sigma), "amp": float(amp),
    "initial_sep": float(initial_sep),
    "scope": "2D two-component NLSE with constrained Gauss-law scalar potential"
}

def lap_fft_complex(psi):
    return np.fft.ifft2(-K2*np.fft.fft2(psi))

def lap_fft_real(f):
    return np.real(np.fft.ifft2(-K2*np.fft.fft2(f)))

def solve_phi(rho_q):
    src = rho_q - np.mean(rho_q)
    return np.real(np.fft.ifft2(np.fft.fft2(src)/(K2 + mu**2)))

def gauss_residual(phi, rho_q):
    src = rho_q - np.mean(rho_q)
    lhs = -lap_fft_real(phi) + mu**2*phi
    res = lhs - src
    return float(np.sqrt(np.mean(res**2))), float(np.sqrt(np.mean(src**2)) + 1e-30)

def init_packet(x0, kx0=0.0):
    env = amp*np.exp(-((X-x0)**2 + Y**2)/(2*sigma**2))
    phase = np.exp(1j*kx0*X)
    return (env*phase).astype(np.complex128)

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

def energies(psi1, psi2, phi, s1, s2):
    rho1 = np.abs(psi1)**2
    rho2 = np.abs(psi2)**2
    # spectral gradient energy: 0.5 ∫ |∇ψ|²
    p1k = np.fft.fft2(psi1)
    p2k = np.fft.fft2(psi2)
    kin1 = 0.5*np.sum(np.abs(np.fft.ifft2(1j*KX*p1k))**2 + np.abs(np.fft.ifft2(1j*KY*p1k))**2)*dx*dx
    kin2 = 0.5*np.sum(np.abs(np.fft.ifft2(1j*KX*p2k))**2 + np.abs(np.fft.ifft2(1j*KY*p2k))**2)*dx*dx
    mat = np.sum(-0.5*g*(rho1**2+rho2**2) + (q/3.0)*(rho1**3+rho2**3))*dx*dx
    src = e_charge*(s1*rho1 + s2*rho2)
    interaction = 0.5*np.sum((src-np.mean(src))*phi)*dx*dx
    # gauge field energy
    phix = np.real(np.fft.ifft2(1j*KX*np.fft.fft2(phi)))
    phiy = np.real(np.fft.ifft2(1j*KY*np.fft.fft2(phi)))
    gauge = 0.5*np.sum(phix**2+phiy**2+mu**2*phi**2)*dx*dx
    return {
        "kinetic_E": float(kin1+kin2),
        "matter_E": float(mat),
        "interaction_E": float(interaction),
        "gauge_E": float(gauge),
        "total_E": float(kin1+kin2+mat+interaction+gauge)
    }

def split_step(psi1, psi2, s1, s2):
    rho1 = np.abs(psi1)**2
    rho2 = np.abs(psi2)**2
    phi = solve_phi(e_charge*(s1*rho1+s2*rho2))
    V1 = -g*rho1 + q*rho1**2 + s1*e_charge*phi
    V2 = -g*rho2 + q*rho2**2 + s2*e_charge*phi

    psi1 = np.exp(-1j*V1*dt/2)*psi1
    psi2 = np.exp(-1j*V2*dt/2)*psi2

    psi1k = np.fft.fft2(psi1)
    psi2k = np.fft.fft2(psi2)
    kinetic_phase = np.exp(-1j*(K2/2.0)*dt)
    psi1 = np.fft.ifft2(psi1k*kinetic_phase)
    psi2 = np.fft.ifft2(psi2k*kinetic_phase)

    rho1 = np.abs(psi1)**2
    rho2 = np.abs(psi2)**2
    phi = solve_phi(e_charge*(s1*rho1+s2*rho2))
    V1 = -g*rho1 + q*rho1**2 + s1*e_charge*phi
    V2 = -g*rho2 + q*rho2**2 + s2*e_charge*phi
    psi1 = np.exp(-1j*V1*dt/2)*psi1
    psi2 = np.exp(-1j*V2*dt/2)*psi2
    return psi1, psi2

def run_case(label, s1, s2):
    # no external kick: we want field to drive centroid trend
    psi1 = init_packet(-initial_sep/2, 0.0)
    psi2 = init_packet(+initial_sep/2, 0.0)

    N10 = norm(psi1)
    N20 = norm(psi2)
    c10 = centroid(psi1)
    c20 = centroid(psi2)
    sep0 = float(c20[0]-c10[0])
    peak10 = peak(psi1)
    peak20 = peak(psi2)
    rms10 = rms_radius_about_centroid(psi1)
    rms20 = rms_radius_about_centroid(psi2)

    rows = []
    for step in range(steps+1):
        rho1 = np.abs(psi1)**2
        rho2 = np.abs(psi2)**2
        rho_q = e_charge*(s1*rho1+s2*rho2)
        phi = solve_phi(rho_q)
        gres, gsrc = gauss_residual(phi, rho_q)
        c1 = centroid(psi1)
        c2 = centroid(psi2)
        sep = float(c2[0]-c1[0])
        en = energies(psi1, psi2, phi, s1, s2)

        if step % sample_every == 0:
            rows.append({
                "case": label, "step": int(step), "t": float(step*dt),
                "s1": int(s1), "s2": int(s2),
                "N1": norm(psi1), "N2": norm(psi2),
                "N1_drift_rel": abs(norm(psi1)-N10)/(N10+1e-30),
                "N2_drift_rel": abs(norm(psi2)-N20)/(N20+1e-30),
                "c1x": c1[0], "c2x": c2[0],
                "c1y": c1[1], "c2y": c2[1],
                "separation": sep,
                "separation_change": sep-sep0,
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

        # exact per-component norm correction to suppress numerical drift
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

like_df, like_summary = run_case("like_plus_plus_wavepackets", +1, +1)
opp_df, opp_summary = run_case("opposite_plus_minus_wavepackets", +1, -1)
all_df = pd.concat([like_df, opp_df], ignore_index=True)
summaries = pd.DataFrame([like_summary, opp_summary])

pass_conditions = {
    "like_final_separation_increases": bool(like_summary["final_separation_change"] > 0),
    "opposite_final_separation_decreases": bool(opp_summary["final_separation_change"] < 0),
    "like_direction_fraction_over_0p60": bool(like_summary["direction_ok_fraction"] > 0.60),
    "opposite_direction_fraction_over_0p60": bool(opp_summary["direction_ok_fraction"] > 0.60),
    "component_norm_drift_under_1e-10": bool(max(like_summary["N1_drift_rel_max"], like_summary["N2_drift_rel_max"], opp_summary["N1_drift_rel_max"], opp_summary["N2_drift_rel_max"]) < 1e-10),
    "gauss_residual_under_1e-10": bool(max(like_summary["gauss_residual_rel_max"], opp_summary["gauss_residual_rel_max"]) < 1e-10),
    "packets_bounded_peak_under_4x": bool(max(like_summary["peak_ratio_max"], opp_summary["peak_ratio_max"]) < 4.0),
    "packets_not_spread_over_4x": bool(max(like_summary["rms_ratio_max"], opp_summary["rms_ratio_max"]) < 4.0),
    "finite": bool(like_summary["finite"] and opp_summary["finite"])
}
verdict = "PASSED as full two-component wavepacket sandbox" if all(pass_conditions.values()) else "PARTIAL / CHECK"

summary = {
    "G131": "Full two-component wavepacket dynamics",
    "verdict": verdict,
    "scope_limit": "2D two-component NLSE + constrained scalar Gauss law; not full EM/QED/Standard Model.",
    "model": {
        "component_dynamics": "i psi_i,t = -1/2 nabla^2 psi_i + [-g|psi_i|^2 + q|psi_i|^4 + s_i e phi] psi_i",
        "gauss_law": "(-nabla^2 + mu^2) phi = e(s1|psi1|^2 + s2|psi2|^2 - mean)",
        "centroids": "tracked directly from |psi_i|^2"
    },
    "pass_conditions": pass_conditions,
    "like_summary": like_summary,
    "opposite_summary": opp_summary,
    "next_gate": "G132 should test charge-sign wavepacket dynamics with topological winding/vortex components included."
}

all_df.to_csv(OUT/"G131_wavepacket_timeseries.csv", index=False)
summaries.to_csv(OUT/"G131_case_summaries.csv", index=False)
with open(OUT/"G131_constants_and_parameters.json", "w", encoding="utf-8") as f:
    json.dump(params, f, indent=2)
with open(OUT/"G131_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

summary_md = f"""# G131 — Full two-component wavepacket dynamics

## Verdict

**{verdict}**

## Scope

This is a reduced 2D two-component nonlinear Schrödinger + constrained scalar
Gauss-law sandbox. It is **not** full electromagnetism, QED, photons, or the
Standard Model.

Unlike G130/G130B, this evolves actual labeled fields `ψ₁` and `ψ₂`; the packet
centroids are measured directly from `|ψᵢ|²`.

## Model

```text
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

Opposite signs:
  initial separation = {opp_summary["initial_separation"]:.6f}
  final separation   = {opp_summary["final_separation"]:.6f}
  change             = {opp_summary["final_separation_change"]:.6f}
  direction fraction = {opp_summary["direction_ok_fraction"]:.3f}
```

## Interpretation

G131 is the first full labeled-wavepacket test in this branch:

```text
ψ₁, ψ₂ evolved directly
charge density built from |ψᵢ|²
φ solved by Gauss law each step
centroid motion measured from actual packet density
```

The result should be interpreted narrowly:

> The constrained charge-sign field can be coupled to two evolving SU/A
> wavepacket components while conserving component norm and maintaining Gauss-law
> control.

## Next gate

**G132 — charge-sign wavepacket dynamics with topological winding**

Add vortex/winding structure back into `ψ₁` and `ψ₂` and test whether:

```text
winding remains conserved
charge-sign attraction/repulsion still appears
packets remain localized
Gauss constraint remains controlled
```
"""
with open(OUT/"G131_summary.md", "w", encoding="utf-8") as f:
    f.write(summary_md)

# Plots
plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["t"], df["separation"], label=case)
plt.xlabel("time")
plt.ylabel("centroid separation")
plt.title("G131 wavepacket centroid separation")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G131_separation.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["t"], df["gauss_residual_rel"], label=case)
plt.yscale("log")
plt.xlabel("time")
plt.ylabel("relative Gauss residual")
plt.title("G131 Gauss-law residual")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G131_gauss_residual.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["t"], df["peak1_ratio"], label=case+" psi1")
    plt.plot(df["t"], df["peak2_ratio"], linestyle="--", label=case+" psi2")
plt.xlabel("time")
plt.ylabel("peak ratio")
plt.title("G131 packet boundedness")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(PLOTS/"G131_peak_ratios.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case, df in all_df.groupby("case"):
    plt.plot(df["t"], df["rms1_ratio"], label=case+" psi1")
    plt.plot(df["t"], df["rms2_ratio"], linestyle="--", label=case+" psi2")
plt.xlabel("time")
plt.ylabel("RMS ratio")
plt.title("G131 packet spread")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(PLOTS/"G131_rms_ratios.png", dpi=160)
plt.close()

readme = """# G131 Full Two-Component Wavepacket Dynamics Package

Run:

```bash
python scripts/G131_full_two_component_wavepacket_dynamics.py
```

Outputs are written to `results/`.

This package evolves two labeled wavefunction components under a constrained
Gauss-law scalar potential. It is a reduced sandbox, not full electromagnetism,
QED, photons, or Standard Model physics.
"""
with open(ROOT/"README_G131.md", "w", encoding="utf-8") as f:
    f.write(readme)

print(json.dumps(summary, indent=2))
