# F6: Model-A Gravitational Decoherence Prediction

## Conceptual setup

Under Model-A, a mass in superposition of two positions creates two distinguishable A-field configurations (one per branch). The combined state is *unresolved A* — no interaction has yet pinned which configuration is the resolved one.

The framework's resolution principle says: physical interactions resolve A. A particle's own A-field couples back to its own mass (it sources A; it also feels grad A). That self-coupling acts as an internal resolution channel, with rate set by the energy of A-distinguishability between the two branches:

```text
tau_decoherence = hbar / E_G
E_G             = (G/2) * integral of [rho_1 - rho_2][rho_1 - rho_2]/|x-y|
```

For two non-overlapping uniform spheres of mass m, radius R, separation d (with d >= 2R):

```text
E_G = G * m^2 * [ 6/(5R) - 1/d ]
```

This is the same closed form as Diosi-Penrose gravitational decoherence — Model-A reproduces it because A-distinguishability *is* the gravitational self-energy of the density difference, just phrased in A-language.

## Calculation

Sweep mass from 10^-22 kg (atomic) to 10^-8 kg (visible dust) at silica density rho = 2200 kg/m^3. Separation set to d = 10 R (well-separated, asymptotic regime). For each mass, compute E_G and tau = hbar / E_G.

## Benchmark predictions

```text
object                                 mass (kg)       R (m)      tau (s)   regime
--------------------------------------------------------------------------------------------------------------
atomic mass (~1 amu)                   1.661e-27    5.65e-11    2.942e+19   indefinitely coherent (no gravitational decoherence visible)
C60 fullerene (720 amu)                1.196e-24    5.06e-10    5.087e+14   indefinitely coherent (no gravitational decoherence visible)
oligo-porphyrin (~25 kDa)              4.151e-23    1.65e-09    1.377e+12   indefinitely coherent (no gravitational decoherence visible)
antibody (~150 kDa)                    2.491e-22    3.00e-09    6.948e+10   indefinitely coherent (no gravitational decoherence visible)
virus capsid (~1 MDa)                  1.661e-21    5.65e-09    2.942e+09   indefinitely coherent (no gravitational decoherence visible)
ribosome (~3 MDa)                      4.982e-21    8.15e-09    4.715e+08   indefinitely coherent (no gravitational decoherence visible)
100 nm silica nanoparticle             1.152e-18    5.00e-08    5.413e+04   indefinitely coherent (no gravitational decoherence visible)
1 micron silica bead                   1.152e-15    5.00e-07    5.413e-01   milliseconds — current cavity-optomechanics frontier
10 micron silica bead                  1.152e-12    5.00e-06    5.413e-06   microseconds — feasible with state-of-the-art controls
100 micron silica bead                 1.152e-09    5.00e-05    5.413e-11   sub-nanosecond — too fast to resolve in current experiments
```

## What the predictions say

- **Atomic / small-molecule regime** (~10^-25 kg and below): predicted tau is enormous (years or longer). Consistent with atomic interferometry showing no gravitational decoherence at these scales.
- **Macromolecule regime** (10^-22 - 10^-19 kg, fullerenes through antibodies): predicted tau is still long compared to typical experiment times. Consistent with current matter-wave interferometry results.
- **Nanoparticle / cavity-optomechanics regime** (10^-18 - 10^-12 kg): predicted tau drops into the millisecond - microsecond window. **This is the regime where Model-A makes a testable prediction.** Cavity-cooled nanoparticles are being prepared in near-superposition states by groups worldwide (Aspelmeyer, Novotny, Kiesel, etc.); observation of decoherence at this rate with the predicted m^2/R scaling would be positive evidence.
- **Larger masses** (>10^-9 kg): predicted tau drops below nanoseconds. These objects effectively cannot be put in superposition in the first place; the prediction is consistent with classical behaviour at macroscopic scales falling out automatically.

## What would distinguish Model-A from "no gravitational decoherence"

Standard quantum mechanics with no explicit gravitational coupling predicts no intrinsic decoherence from gravity — only environmental decoherence (gas collisions, thermal photons). For a perfectly isolated nanoparticle in deep vacuum at low temperature, standard QM says coherence persists indefinitely.

Model-A predicts the opposite: even in perfect isolation, the self-gravitational A-coupling resolves the superposition on the timescale tabulated above. The discriminating experiment is:

1. Prepare a nanoparticle in a center-of-mass superposition of two positions separated by ~10 R.
2. Isolate it well enough that environmental decoherence times exceed the predicted tau by at least an order of magnitude.
3. Measure the coherence decay.
4. If decay matches predicted tau scaling with mass and radius, the unified A picture has positive support. If decay is faster, environment dominates. If decay is slower, the prediction fails.

This is the same experimental target as Diosi-Penrose, so any experimental outcome that bears on DP also bears on Model-A in this regime. Model-A's distinction comes in when we go to stronger fields (where the modified g_rr matters) or to the horizon regime (where the resolution-status aspect of A dominates and gives Hawking-style emission).

## Honest caveats

- The closed form assumes a uniform-density sphere. Real nanoparticles have surface roughness, internal density variations, and shape factors that change E_G by O(1).
- The separation choice d = 10 R is specific. Experimental superpositions often have d much smaller than R, which *reduces* E_G (the two configurations are barely distinguishable), increasing tau accordingly. The prediction is most robust in the d >> R asymptotic regime.
- E_G as written includes the self-energy contribution 6/(5R) which dominates for any d > R. This is the part that does NOT depend on the separation, only on the particle's compactness. Some formulations of gravitational decoherence drop this term and keep only the cross-term. Model-A as written keeps both.
- This script does not derive E_G from a Lagrangian for A. The form is taken from the Diosi-Penrose result, with Model-A providing the conceptual reinterpretation (A-distinguishability) rather than a new derivation. A field-theoretic derivation of E_G from an A-action remains an open theoretical task — the same gap the README flags for the broader framework.

## Generated plots

- `plots/F6_decoherence_tau_vs_mass.png`
- `plots/F6_decoherence_E_G_vs_mass.png`
