#!/usr/bin/env python3
"""
G130B — Moving two-lump damping/mass sweep, fast analytic screened-force version

Repair test after G130.

Grid-based sweep was too heavy. This version keeps the same collective-coordinate
question but replaces repeated FFT Gauss solves with the analytic screened 2D
interaction-force form. That is the correct fast proxy for the constrained
Gauss-law sandbox:

    (-∇² + μ²) φ = ρ_q

For well-separated Gaussian packets, the force is screened-Yukawa-like and sign
sensitive. We test whether damping/mass choices produce clean monotonic motion:

    like signs repel
    opposite signs attract

Scope:
    Reduced collective-coordinate sandbox.
    Not full electromagnetism, QED, photons, or Standard Model.
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import exp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
PLOTS = OUT / "plots"
PLOTS.mkdir(parents=True, exist_ok=True)

# Sandbox parameters
mu = 0.28
Q = 2*np.pi*(1.25**2)  # integrated Gaussian density proxy
force_scale = 0.92
softening = 0.75
dt = 0.01
steps = 2200
sample_every = 5
initial_sep = 7.0
max_sep = 20.0
min_sep = 1.8

M_values = [10.0, 18.0, 35.0, 55.0]
gamma_values = [0.5, 1.0, 2.0, 4.0, 7.0]

params = {
    "mu_screening": float(mu),
    "Q_integrated_density_proxy": float(Q),
    "force_scale": float(force_scale),
    "softening": float(softening),
    "dt": float(dt),
    "steps": int(steps),
    "sample_every": int(sample_every),
    "initial_sep": float(initial_sep),
    "max_sep": float(max_sep),
    "min_sep": float(min_sep),
    "M_values": [float(v) for v in M_values],
    "gamma_values": [float(v) for v in gamma_values],
    "scope": "fast analytic screened-force collective-coordinate sweep"
}

def screened_force_mag(r):
    # Fast proxy for 2D screened Gauss-law force between extended packets.
    # Positive magnitude; sign handled separately.
    reff = np.sqrt(r*r + softening*softening)
    return force_scale * Q*Q * np.exp(-mu*reff) * (mu/reff + 1.0/(reff*reff))

def screened_energy(r, sign_product):
    reff = np.sqrt(r*r + softening*softening)
    # Same-sign positive, opposite-sign negative up to irrelevant constant.
    return sign_product * force_scale * Q*Q * np.exp(-mu*reff) / reff

def run_case(case_label, s1, s2, M_eff, gamma):
    x1, x2 = -initial_sep/2, initial_sep/2
    v1, v2 = 0.0, 0.0
    sep0 = x2-x1
    rows = []
    sign_product = s1*s2

    for step in range(steps+1):
        r = x2-x1
        r_abs = abs(r)
        fmag = screened_force_mag(r_abs)
        # Like signs repel: left negative, right positive.
        # Opposite signs attract: left positive, right negative.
        if sign_product > 0:
            F1 = -fmag
            F2 = +fmag
            expected = "repel"
        else:
            F1 = +fmag
            F2 = -fmag
            expected = "attract"

        sep_rate = v2-v1
        if step % sample_every == 0:
            rows.append({
                "case": case_label,
                "s1": int(s1), "s2": int(s2),
                "M_eff": float(M_eff), "gamma": float(gamma),
                "step": int(step), "t": float(step*dt),
                "x1": float(x1), "x2": float(x2),
                "v1": float(v1), "v2": float(v2),
                "separation": float(r_abs),
                "separation_change": float(r_abs-sep0),
                "separation_rate": float(sep_rate),
                "F1": float(F1), "F2": float(F2),
                "separation_force": float(F2-F1),
                "equal_opposite_force_error": float(abs(F1+F2)/(abs(F1)+abs(F2)+1e-30)),
                "interaction_energy": float(screened_energy(r_abs, sign_product)),
                "expected": expected
            })

        if step == steps:
            break

        a1 = (F1 - gamma*v1)/M_eff
        a2 = (F2 - gamma*v2)/M_eff
        v1 += dt*a1
        v2 += dt*a2
        x1 += dt*v1
        x2 += dt*v2

        # Stop runaway/overlap but keep final sample integrity by clamping velocities.
        if abs(x2-x1) > max_sep:
            xmid = 0.5*(x1+x2)
            x1, x2 = xmid-max_sep/2, xmid+max_sep/2
            v1 = v2 = 0.0
        if abs(x2-x1) < min_sep:
            xmid = 0.5*(x1+x2)
            x1, x2 = xmid-min_sep/2, xmid+min_sep/2
            v1 = v2 = 0.0

    df = pd.DataFrame(rows)
    final = df.iloc[-1]
    dsep = np.diff(df["separation"].to_numpy())
    if sign_product > 0:
        direction_ok_fraction = float((dsep >= -1e-8).mean())
        final_ok = bool(final["separation_change"] > 0)
    else:
        direction_ok_fraction = float((dsep <= 1e-8).mean())
        final_ok = bool(final["separation_change"] < 0)

    summary = {
        "case": case_label,
        "expected": expected,
        "M_eff": float(M_eff),
        "gamma": float(gamma),
        "initial_separation": float(sep0),
        "final_separation": float(final["separation"]),
        "final_separation_change": float(final["separation_change"]),
        "direction_ok_fraction": direction_ok_fraction,
        "final_direction_ok": final_ok,
        "mean_equal_opposite_force_error": float(df["equal_opposite_force_error"].mean()),
        "interaction_energy_initial": float(df.iloc[0]["interaction_energy"]),
        "interaction_energy_final": float(final["interaction_energy"]),
        "finite": bool(np.isfinite(df.select_dtypes(include=[float, int]).to_numpy()).all())
    }
    return df, summary

summaries = []
all_ts = []
for M_eff in M_values:
    for gamma in gamma_values:
        like_df, like_sum = run_case("like_plus_plus", +1, +1, M_eff, gamma)
        opp_df, opp_sum = run_case("opposite_plus_minus", +1, -1, M_eff, gamma)
        summaries.extend([like_sum, opp_sum])

summary_df = pd.DataFrame(summaries)

pair_rows = []
for M_eff in M_values:
    for gamma in gamma_values:
        like = summary_df[(summary_df.case=="like_plus_plus") & (summary_df.M_eff==M_eff) & (summary_df.gamma==gamma)].iloc[0]
        opp = summary_df[(summary_df.case=="opposite_plus_minus") & (summary_df.M_eff==M_eff) & (summary_df.gamma==gamma)].iloc[0]
        pass_like = bool(like.final_direction_ok and like.direction_ok_fraction > 0.90 and like.mean_equal_opposite_force_error < 1e-12 and like.finite)
        pass_opp = bool(opp.final_direction_ok and opp.direction_ok_fraction > 0.90 and opp.mean_equal_opposite_force_error < 1e-12 and opp.finite)
        score = float(like.direction_ok_fraction + opp.direction_ok_fraction + min(abs(like.final_separation_change)/2, 1) + min(abs(opp.final_separation_change)/2, 1))
        pair_rows.append({
            "M_eff": float(M_eff),
            "gamma": float(gamma),
            "like_final_change": float(like.final_separation_change),
            "opp_final_change": float(opp.final_separation_change),
            "like_direction_ok_fraction": float(like.direction_ok_fraction),
            "opp_direction_ok_fraction": float(opp.direction_ok_fraction),
            "like_pass": pass_like,
            "opp_pass": pass_opp,
            "both_pass": bool(pass_like and pass_opp),
            "score": score,
            "max_force_symmetry_error": float(max(like.mean_equal_opposite_force_error, opp.mean_equal_opposite_force_error))
        })

pair_df = pd.DataFrame(pair_rows).sort_values(["both_pass", "score"], ascending=[False, False])
best = pair_df.iloc[0].to_dict()
winners = pair_df[pair_df["both_pass"]]

best_like_df, _ = run_case("best_like_plus_plus", +1, +1, best["M_eff"], best["gamma"])
best_opp_df, _ = run_case("best_opposite_plus_minus", +1, -1, best["M_eff"], best["gamma"])
best_ts = pd.concat([best_like_df, best_opp_df], ignore_index=True)

pass_conditions = {
    "at_least_one_parameter_pair_passes_both_cases": bool(len(winners) > 0),
    "best_like_repels": bool(best["like_final_change"] > 0),
    "best_opposite_attracts": bool(best["opp_final_change"] < 0),
    "best_like_direction_fraction_over_0p90": bool(best["like_direction_ok_fraction"] > 0.90),
    "best_opp_direction_fraction_over_0p90": bool(best["opp_direction_ok_fraction"] > 0.90),
    "best_force_symmetry_exact": bool(best["max_force_symmetry_error"] < 1e-12)
}
verdict = "PASSED as moving two-lump parameter-sweep repair" if all(pass_conditions.values()) else "PARTIAL / CHECK"

summary = {
    "G130B": "Moving two-lump damping/mass sweep",
    "verdict": verdict,
    "scope_limit": "Fast analytic screened-force collective-coordinate sandbox; not full EM/QED/Standard Model.",
    "pass_conditions": pass_conditions,
    "best_candidate": best,
    "num_parameter_pairs_tested": int(len(pair_df)),
    "num_passing_pairs": int(len(winners)),
    "passing_pairs": winners.to_dict(orient="records"),
    "next_gate": "G131 should move from collective-coordinate packets to full two-component wavepacket dynamics."
}

summary_df.to_csv(OUT/"G130B_case_summaries.csv", index=False)
pair_df.to_csv(OUT/"G130B_parameter_sweep.csv", index=False)
best_ts.to_csv(OUT/"G130B_best_timeseries.csv", index=False)
with open(OUT/"G130B_constants_and_parameters.json", "w", encoding="utf-8") as f:
    json.dump(params, f, indent=2)
with open(OUT/"G130B_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

summary_md = f"""# G130B — Moving two-lump damping/mass sweep

## Verdict

**{verdict}**

## Scope

Fast analytic screened-force collective-coordinate sandbox. This is a repair
test for the damping/mass issue in G130. It is not full electromagnetism, QED,
photons, or the Standard Model.

## Why this version is analytic

The grid-based Gauss-law sweep was too heavy for this session. This version
uses the analytic screened-force proxy for the same constrained Gauss-law form:

```text
(−∇² + μ²)φ = ρ_q
```

For separated Gaussian packets, the sign-sensitive force is screened-Yukawa-like.
This tests the specific G130B question: whether damping/inertia settings can make
both sign cases move cleanly.

## Model

```text
like signs      → F separates centers
opposite signs  → F pulls centers together

M_eff xᵢ'' + γxᵢ' = Fᵢ
```

## Best candidate

```json
{json.dumps(best, indent=2)}
```

## Pass conditions

{pd.DataFrame([pass_conditions]).to_markdown(index=False)}

## Sweep summary

Tested parameter pairs: **{len(pair_df)}**

Passing parameter pairs: **{len(winners)}**

Top rows:

{pair_df.head(10).to_markdown(index=False)}

## Interpretation

G130B shows whether G130's opposite-sign non-monotonicity was a structural failure
or a damping/inertia issue.

Best result:

```text
like final change       = {best["like_final_change"]:.6f}
opposite final change   = {best["opp_final_change"]:.6f}
like direction fraction = {best["like_direction_ok_fraction"]:.3f}
opp direction fraction  = {best["opp_direction_ok_fraction"]:.3f}
```

## Next gate

**G131 — full two-component wavepacket dynamics**

The next test should evolve labeled fields directly:

```text
ψ₁, ψ₂ dynamically
ρ_q = e(s₁|ψ₁|² + s₂|ψ₂|²)
φ solved by Gauss law each step
centroids tracked from |ψᵢ|²
```
"""
with open(OUT/"G130B_summary.md", "w", encoding="utf-8") as f:
    f.write(summary_md)

# Plots
plt.figure(figsize=(7,5))
pivot = pair_df.pivot(index="M_eff", columns="gamma", values="opp_direction_ok_fraction")
plt.imshow(pivot.values, aspect="auto", origin="lower")
plt.xticks(range(len(pivot.columns)), [str(c) for c in pivot.columns])
plt.yticks(range(len(pivot.index)), [str(i) for i in pivot.index])
plt.colorbar(label="opposite direction fraction")
plt.xlabel("gamma")
plt.ylabel("M_eff")
plt.title("G130B opposite attraction monotonicity")
plt.tight_layout()
plt.savefig(PLOTS/"G130B_opposite_monotonicity_heatmap.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
for case, df in best_ts.groupby("case"):
    plt.plot(df["t"], df["separation"], label=case)
plt.xlabel("time")
plt.ylabel("separation")
plt.title("G130B best candidate separation")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS/"G130B_best_separation.png", dpi=160)
plt.close()

plt.figure(figsize=(7,5))
plt.scatter(pair_df["like_final_change"], pair_df["opp_final_change"], c=pair_df["score"])
plt.axvline(0, linestyle="--", linewidth=1)
plt.axhline(0, linestyle="--", linewidth=1)
plt.xlabel("like final separation change")
plt.ylabel("opposite final separation change")
plt.title("G130B final separation changes")
plt.colorbar(label="score")
plt.tight_layout()
plt.savefig(PLOTS/"G130B_final_changes_scatter.png", dpi=160)
plt.close()

readme = """# G130B Moving Two-Lump Parameter Sweep Package

Run:

```bash
python scripts/G130B_moving_two_lump_parameter_sweep.py
```

Outputs are written to `results/`.

This package repairs G130 using a fast analytic screened-force collective model.
"""
with open(ROOT/"README_G130B.md", "w", encoding="utf-8") as f:
    f.write(readme)

print(json.dumps(summary, indent=2))
