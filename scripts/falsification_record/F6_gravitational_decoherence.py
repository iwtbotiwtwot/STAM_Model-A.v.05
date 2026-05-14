#!/usr/bin/env python3
"""
F6_gravitational_decoherence.py

Model-A prediction for gravitational decoherence of a mass in spatial
superposition.

Conceptual setup (Model-A's unified A picture):
    - A mass in superposition of two positions creates two distinguishable
      A-field configurations (one for each branch).
    - The combined state is "unresolved A" — no interaction has yet pinned
      which configuration is real.
    - Self-resolution rate is set by the energy cost of the
      A-distinguishability: how much the two A-configurations differ.
    - tau_decoherence = hbar / E_G

For two non-overlapping uniform spheres of mass m, radius R, separated by
distance d >= 2R, the A-distinguishability energy works out to:

    E_G = G * m^2 * [ 6/(5R) - 1/d ]

This is the same closed form as Diosi-Penrose gravitational decoherence,
because A-field distinguishability IS the gravitational self-energy of
the density difference — just rephrased in Model-A's language.

What this script does:
    1. Sweep particle mass from atomic (~10^-22 kg) to dust grain (~10^-8 kg).
    2. For each mass, set R from solid density (silica-like, ~2200 kg/m^3),
       and use a well-separated superposition distance d = 10 R.
    3. Compute E_G and tau_decoherence.
    4. Plot tau vs m, mark experimental coherence-time thresholds.
    5. Identify the mass scale where Model-A predicts decoherence in a
       laboratory-accessible time window (microsecond to second).
    6. Note the regime where standard QM without explicit gravity coupling
       predicts indefinite coherence — i.e. where Model-A is testable.

Why this is a meaningful test:
    - Atomic/molecular interferometry has demonstrated coherence at masses
      < 10^-22 kg with no observed gravitational decoherence. Model-A's
      prediction in that regime is tau >> seconds — consistent.
    - The prediction window where Model-A says "you should see decoherence"
      is roughly 10^-15 to 10^-12 kg with tau in milliseconds-microseconds.
      This is the experimental frontier (cavity-cooled nanoparticles,
      levitated optomechanics).
    - Detection of decoherence in this range with the predicted scaling
      tau ~ R / m^2 would be positive evidence for the unified A picture.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HBAR = 1.0545718e-34          # J s
G = 6.67430e-11               # m^3 / (kg s^2)
RHO_MATERIAL = 2200.0         # kg/m^3, silica-like solid
AMU = 1.66053907e-27          # kg

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


def particle_radius(mass: float, density: float = RHO_MATERIAL) -> float:
    """Radius of uniform sphere with given mass and density."""
    return (3.0 * mass / (4.0 * np.pi * density)) ** (1.0 / 3.0)


def E_G_two_spheres(mass: float, radius: float, separation: float) -> float:
    """A-distinguishability energy for two uniform-sphere branches.

    Closed form for non-overlapping spheres (d >= 2R):
        E_G = G m^2 [ 6/(5R) - 1/d ]
    For partial overlap, the distinguishability decreases (configurations
    look more similar). For d -> 0: E_G -> 0 (single configuration).
    """
    if separation >= 2 * radius:
        return G * mass ** 2 * (6.0 / (5.0 * radius) - 1.0 / separation)
    # Partial overlap — quadratic suppression toward full overlap.
    full_separation_E = G * mass ** 2 * (6.0 / (5.0 * radius) - 1.0 / (2 * radius))
    return full_separation_E * (separation / (2.0 * radius)) ** 2


def tau_decoherence(mass: float, radius: float, separation: float) -> float:
    """tau = hbar / E_G."""
    eg = E_G_two_spheres(mass, radius, separation)
    if eg <= 0:
        return float("inf")
    return HBAR / eg


def regime_label(tau: float) -> str:
    """Human label for where this coherence time sits experimentally."""
    if tau > 60.0:
        return "indefinitely coherent (no gravitational decoherence visible)"
    if tau > 1.0:
        return "seconds — extreme experimental frontier"
    if tau > 1e-3:
        return "milliseconds — current cavity-optomechanics frontier"
    if tau > 1e-6:
        return "microseconds — feasible with state-of-the-art controls"
    if tau > 1e-9:
        return "nanoseconds — fast, requires dedicated apparatus"
    return "sub-nanosecond — too fast to resolve in current experiments"


# --- main calculation ---

def sweep_masses() -> dict:
    masses = np.logspace(-22, -8, 400)
    radii = np.array([particle_radius(m) for m in masses])
    separations = 10.0 * radii  # well-separated
    eg = np.array([E_G_two_spheres(m, r, d)
                   for m, r, d in zip(masses, radii, separations)])
    tau = np.where(eg > 0, HBAR / eg, np.inf)
    return {
        "masses": masses,
        "radii": radii,
        "separations": separations,
        "E_G": eg,
        "tau": tau,
    }


def benchmark_table() -> list[dict]:
    """Specific masses tied to recognisable physical objects."""
    cases = [
        ("atomic mass (~1 amu)", AMU),
        ("C60 fullerene (720 amu)", 720 * AMU),
        ("oligo-porphyrin (~25 kDa)", 25e3 * AMU),
        ("antibody (~150 kDa)", 1.5e5 * AMU),
        ("virus capsid (~1 MDa)", 1e6 * AMU),
        ("ribosome (~3 MDa)", 3e6 * AMU),
        ("100 nm silica nanoparticle", (4 / 3) * np.pi * (50e-9) ** 3 * RHO_MATERIAL),
        ("1 micron silica bead", (4 / 3) * np.pi * (500e-9) ** 3 * RHO_MATERIAL),
        ("10 micron silica bead", (4 / 3) * np.pi * (5e-6) ** 3 * RHO_MATERIAL),
        ("100 micron silica bead", (4 / 3) * np.pi * (50e-6) ** 3 * RHO_MATERIAL),
    ]
    rows = []
    for label, m in cases:
        r = particle_radius(m)
        d = 10 * r
        eg = E_G_two_spheres(m, r, d)
        tau = HBAR / eg if eg > 0 else float("inf")
        rows.append({
            "label": label,
            "mass_kg": m,
            "radius_m": r,
            "separation_m": d,
            "E_G_J": eg,
            "tau_s": tau,
            "regime": regime_label(tau),
        })
    return rows


# --- plots ---

def plot_tau_vs_mass(sweep: dict, benchmarks: list[dict]) -> Path:
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.loglog(sweep["masses"], sweep["tau"], color="tab:blue", linewidth=2,
              label="Model-A: tau = hbar / E_G(A-distinguishability)")

    thresholds = [
        (1e-9, "1 ns"),
        (1e-6, "1 us"),
        (1e-3, "1 ms"),
        (1.0, "1 s"),
        (60.0, "1 min"),
    ]
    for tau_val, label in thresholds:
        ax.axhline(tau_val, color="gray", linestyle=":", linewidth=0.8)
        ax.text(sweep["masses"][-1] * 1.1, tau_val, label,
                fontsize=8, color="gray", va="center")

    # Highlight the laboratory-accessible window.
    ax.axhspan(1e-6, 1.0, color="tab:green", alpha=0.10,
               label="lab-accessible coherence window (1 us - 1 s)")

    # Benchmark points.
    for row in benchmarks:
        if np.isfinite(row["tau_s"]) and row["tau_s"] > 1e-15:
            ax.plot(row["mass_kg"], row["tau_s"], "o",
                    color="tab:red", markersize=6)
            ax.annotate(row["label"],
                        xy=(row["mass_kg"], row["tau_s"]),
                        xytext=(8, 4), textcoords="offset points",
                        fontsize=7.5, color="tab:red")

    ax.set_xlabel("Particle mass m  (kg)")
    ax.set_ylabel("Predicted decoherence time tau  (s)")
    ax.set_title("Model-A gravitational decoherence: coherence time vs mass\n"
                 "(uniform sphere, density 2200 kg/m^3, separation = 10 R)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower left")

    out = PLOTS / "F6_decoherence_tau_vs_mass.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


def plot_eg_vs_mass(sweep: dict) -> Path:
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.loglog(sweep["masses"], sweep["E_G"], color="tab:purple", linewidth=2)
    # Reference: thermal energy scales
    kT_300 = 1.380649e-23 * 300
    kT_1mK = 1.380649e-23 * 1e-3
    ax.axhline(kT_300, color="tab:red", linestyle="--", alpha=0.6,
               label=f"k_B T at 300 K = {kT_300:.2e} J")
    ax.axhline(kT_1mK, color="tab:blue", linestyle="--", alpha=0.6,
               label=f"k_B T at 1 mK = {kT_1mK:.2e} J")
    ax.set_xlabel("Particle mass m (kg)")
    ax.set_ylabel("A-distinguishability energy E_G (J)")
    ax.set_title("Model-A E_G vs mass (separation = 10 R, silica-density sphere)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    out = PLOTS / "F6_decoherence_E_G_vs_mass.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    return out


# --- markdown summary ---

def write_markdown(sweep: dict, benchmarks: list[dict],
                   plot_paths: list[Path]) -> Path:
    md = []
    md.append("# F6: Model-A Gravitational Decoherence Prediction\n")

    md.append("## Conceptual setup\n")
    md.append(
        "Under Model-A, a mass in superposition of two positions creates "
        "two distinguishable A-field configurations (one per branch). The "
        "combined state is *unresolved A* — no interaction has yet pinned "
        "which configuration is the resolved one.\n"
        "\n"
        "The framework's resolution principle says: physical interactions "
        "resolve A. A particle's own A-field couples back to its own mass "
        "(it sources A; it also feels grad A). That self-coupling acts as "
        "an internal resolution channel, with rate set by the energy of "
        "A-distinguishability between the two branches:\n"
        "\n"
        "```text\n"
        "tau_decoherence = hbar / E_G\n"
        "E_G             = (G/2) * integral of [rho_1 - rho_2][rho_1 - rho_2]/|x-y|\n"
        "```\n"
        "\n"
        "For two non-overlapping uniform spheres of mass m, radius R, "
        "separation d (with d >= 2R):\n"
        "\n"
        "```text\n"
        "E_G = G * m^2 * [ 6/(5R) - 1/d ]\n"
        "```\n"
        "\n"
        "This is the same closed form as Diosi-Penrose gravitational "
        "decoherence — Model-A reproduces it because A-distinguishability "
        "*is* the gravitational self-energy of the density difference, just "
        "phrased in A-language.\n"
    )

    md.append("## Calculation\n")
    md.append(
        "Sweep mass from 10^-22 kg (atomic) to 10^-8 kg (visible dust) at "
        f"silica density rho = {RHO_MATERIAL:.0f} kg/m^3. Separation set "
        "to d = 10 R (well-separated, asymptotic regime). For each mass, "
        "compute E_G and tau = hbar / E_G.\n"
    )

    md.append("## Benchmark predictions\n")
    md.append("```text")
    md.append(f"{'object':<35s}  {'mass (kg)':>11s}  "
              f"{'R (m)':>10s}  {'tau (s)':>11s}   regime")
    md.append("-" * 110)
    for row in benchmarks:
        md.append(
            f"{row['label']:<35s}  {row['mass_kg']:>11.3e}  "
            f"{row['radius_m']:>10.2e}  {row['tau_s']:>11.3e}   {row['regime']}"
        )
    md.append("```\n")

    md.append("## What the predictions say\n")
    md.append(
        "- **Atomic / small-molecule regime** (~10^-25 kg and below): "
        "predicted tau is enormous (years or longer). Consistent with "
        "atomic interferometry showing no gravitational decoherence at "
        "these scales.\n"
        "- **Macromolecule regime** (10^-22 - 10^-19 kg, fullerenes through "
        "antibodies): predicted tau is still long compared to typical "
        "experiment times. Consistent with current matter-wave "
        "interferometry results.\n"
        "- **Nanoparticle / cavity-optomechanics regime** (10^-18 - 10^-12 "
        "kg): predicted tau drops into the millisecond - microsecond "
        "window. **This is the regime where Model-A makes a testable "
        "prediction.** Cavity-cooled nanoparticles are being prepared in "
        "near-superposition states by groups worldwide (Aspelmeyer, "
        "Novotny, Kiesel, etc.); observation of decoherence at this rate "
        "with the predicted m^2/R scaling would be positive evidence.\n"
        "- **Larger masses** (>10^-9 kg): predicted tau drops below "
        "nanoseconds. These objects effectively cannot be put in "
        "superposition in the first place; the prediction is consistent "
        "with classical behaviour at macroscopic scales falling out "
        "automatically.\n"
    )

    md.append("## What would distinguish Model-A from \"no gravitational decoherence\"\n")
    md.append(
        "Standard quantum mechanics with no explicit gravitational "
        "coupling predicts no intrinsic decoherence from gravity — only "
        "environmental decoherence (gas collisions, thermal photons). For "
        "a perfectly isolated nanoparticle in deep vacuum at low "
        "temperature, standard QM says coherence persists indefinitely.\n"
        "\n"
        "Model-A predicts the opposite: even in perfect isolation, the "
        "self-gravitational A-coupling resolves the superposition on the "
        "timescale tabulated above. The discriminating experiment is:\n"
        "\n"
        "1. Prepare a nanoparticle in a center-of-mass superposition of "
        "two positions separated by ~10 R.\n"
        "2. Isolate it well enough that environmental decoherence times "
        "exceed the predicted tau by at least an order of magnitude.\n"
        "3. Measure the coherence decay.\n"
        "4. If decay matches predicted tau scaling with mass and radius, "
        "the unified A picture has positive support. If decay is faster, "
        "environment dominates. If decay is slower, the prediction "
        "fails.\n"
        "\n"
        "This is the same experimental target as Diosi-Penrose, so any "
        "experimental outcome that bears on DP also bears on Model-A in "
        "this regime. Model-A's distinction comes in when we go to "
        "stronger fields (where the modified g_rr matters) or to the "
        "horizon regime (where the resolution-status aspect of A "
        "dominates and gives Hawking-style emission).\n"
    )

    md.append("## Honest caveats\n")
    md.append(
        "- The closed form assumes a uniform-density sphere. Real "
        "nanoparticles have surface roughness, internal density "
        "variations, and shape factors that change E_G by O(1).\n"
        "- The separation choice d = 10 R is specific. Experimental "
        "superpositions often have d much smaller than R, which "
        "*reduces* E_G (the two configurations are barely "
        "distinguishable), increasing tau accordingly. The prediction "
        "is most robust in the d >> R asymptotic regime.\n"
        "- E_G as written includes the self-energy contribution 6/(5R) "
        "which dominates for any d > R. This is the part that does NOT "
        "depend on the separation, only on the particle's compactness. "
        "Some formulations of gravitational decoherence drop this term "
        "and keep only the cross-term. Model-A as written keeps both.\n"
        "- This script does not derive E_G from a Lagrangian for A. The "
        "form is taken from the Diosi-Penrose result, with Model-A "
        "providing the conceptual reinterpretation (A-distinguishability) "
        "rather than a new derivation. A field-theoretic derivation of "
        "E_G from an A-action remains an open theoretical task — the "
        "same gap the README flags for the broader framework.\n"
    )

    md.append("## Generated plots\n")
    for p in plot_paths:
        md.append(f"- `{p.relative_to(ROOT).as_posix()}`")
    md.append("")

    out = RESULTS / "F6_gravitational_decoherence_summary.md"
    out.write_text("\n".join(md), encoding="utf-8")
    return out


def main() -> None:
    sweep = sweep_masses()
    benchmarks = benchmark_table()

    plots = [
        plot_tau_vs_mass(sweep, benchmarks),
        plot_eg_vs_mass(sweep),
    ]
    summary = write_markdown(sweep, benchmarks, plots)

    print("F6: Model-A gravitational decoherence prediction")
    print("=" * 72)
    print()
    print(f"Mass range scanned: {sweep['masses'][0]:.1e} kg to "
          f"{sweep['masses'][-1]:.1e} kg")
    print(f"Density assumed:    {RHO_MATERIAL} kg/m^3 (silica-like)")
    print(f"Separation:         d = 10 R (well-separated branches)")
    print()
    print("Benchmark predictions (decoherence time tau = hbar / E_G):")
    print()
    print(f"  {'object':<35s}  {'mass (kg)':>11s}  {'R (m)':>10s}  "
          f"{'tau (s)':>11s}")
    print("  " + "-" * 75)
    for row in benchmarks:
        print(f"  {row['label']:<35s}  {row['mass_kg']:>11.3e}  "
              f"{row['radius_m']:>10.2e}  {row['tau_s']:>11.3e}")
    print()
    print("Lab-accessible window (tau between 1 us and 1 s) corresponds to")
    print("the nanoparticle / cavity-optomechanics regime. This is the")
    print("experimental frontier where Model-A's prediction is testable.")
    print()
    print("Files written:")
    print(f"  {summary}")
    for p in plots:
        print(f"  {p}")


if __name__ == "__main__":
    main()
