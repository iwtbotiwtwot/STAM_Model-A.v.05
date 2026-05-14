# F5c: GW170817 Multi-Messenger Test of Model-B's 'GW = c, light slowed by A' Claim

## The question

F5b's summary asserts that GW170817's ~1.74 s arrival offset between the GW signal and the gamma-ray counterpart is consistent with Model-B **by construction** — GWs propagate at c intrinsically, and light is slowed only by `(1/c) * integral A ds` along the path. This script tests that claim numerically: take the strong form of the claim at face value, compute the predicted delay over the 40 Mpc path to NGC 4993, and compare to 1.74 s.

## Setup

Model-B propagation:
```text
Delta_t_GW    = D / c                          (no A coupling)
Delta_t_light = D / c + (1/c) * int A(s) ds    (Shapiro-slowed)
Differential  = (1/c) * int A(s) ds
```

A is matter-sourced (Framework C), so each mass along the line of sight contributes A(r) = 2GM/(c^2 r). Two dominant contributions:

- **Milky Way**: M = 1.0e+12 M_sun at the galactic center; photon path with impact parameter b = R_sun * sin(70 deg) = 7.52 kpc.
  Integral closed form: (2GM/c^2) * asinh(D/b).

- **NGC 4993**: M = 4.0e+10 M_sun; kilonova at 2.0 kpc offset, photon traveling radially outward.
  Integral closed form: (2GM/c^2) * ln(D/r_off).

Both integrals are logarithmically divergent for an isolated point mass; we cut off at the source distance D = 40 Mpc.

## Results

```text
MW   integral A ds = 2.739e+16 m   -> Delta_t = 9.137e+07 s  (2.90 yr)
NGC  integral A ds = 1.170e+15 m   -> Delta_t = 3.903e+06 s  (45.2 days)
TOTAL              = 2.856e+16 m   -> Delta_t = 9.527e+07 s  (3.02 yr)
```

Observed GW170817 delay: **1.74 s**.

Ratio predicted / observed: **5.48e+07**.

## What value of A would match the observation?

Inverting `Delta_t = (1/c) * <A> * D`:
```text
<A>_required = c * Delta_t / D = 4.226e-16
A_MW at the Sun (reference)   = 1.197e-05
```

The path-average A required to reproduce 1.74 s is about 2.8e+10 times **smaller** than the local Milky Way potential A at the Sun. For a photon to traverse the Galaxy and have an effective path-average A this small, light must NOT couple to the local matter-sourced A as written.

## Verdict

**The strong-form Model-B claim — 'GW always c; light slowed by the full matter-sourced A integral' — predicts a GW-vs-light arrival offset of ~years for GW170817. Observed: 1.74 s. The model is INCONSISTENT with the data by roughly eight orders of magnitude.**

The reason standard GR matches GW170817 cleanly is that BOTH messengers experience identical Shapiro delay through galactic potentials, so the differential cancels and only the astrophysical jet-launch ~1.7 s remains. Any framework that decouples GW from the gravitational potential while keeping light coupled will, generically, predict large differential delays from local potentials.

F5b's 'passes by construction' claim does not survive a numerical evaluation under the construction it specifies.

## What Model-B has to do to actually pass GW170817

Three structural options, each with consequences:

**Option A — A is *not* the full Newtonian potential.** Restrict A to a cosmological background (or a mode that decouples from local virialized matter), so light traveling through galactic potentials sees ordinary GR Shapiro delay (felt equally by GW), and only a tiny extra A-induced delay from cosmological background remains. **Cost:** breaks the README's local A results — GPS clock-rate offset, solar-system Shapiro delay, the gravity bridge `g = (c^2/2) grad A` — all of which assume A = 2GM/(c^2 r) for local masses.

**Option B — GW also feels A, identically to light.** Restore messenger equivalence; the differential delay collapses to the astrophysical jet-launch ~1.7 s. **Cost:** Model-B's narrative claim that 'GW = space, light = matter on space' becomes purely interpretive — there is no observable predictive distinction between the two messengers' propagation.

**Option C — A as a perturbation, not the full potential.** Define A as the deviation from a fiducial GR background, so local matter still sources the GR metric (and both messengers feel it equally), with A capturing only an additional small-amplitude scalar mode. **Cost:** real theoretical work — specifying the coupling structure is essentially the open problem flagged in the README ('Deriving the full field equations').

Until one of these is committed and the local A-based calculations re-derived consistently, GW170817's 1.74 s remains an unresolved tension for Model-B as currently stated.

## Caveats and assumptions

- Milky Way modeled as a point mass at the galactic center. A realistic NFW halo would change the prefactor by O(1) but not the order of magnitude.
- NGC 4993 modeled as a point mass at the galaxy center. The kilonova offset of 2 kpc is consistent with observed offsets of short-GRB hosts; varying it by an order of magnitude moves the NGC delay by a factor of ~ln(10) = 2.3, not enough to change the verdict.
- Intergalactic structure (Local Group, Virgo Cluster, intervening haloes) is neglected. Including them would only increase the predicted delay.
- The line-of-sight geometry uses a single representative angle (70 deg from Sun-GC). Realistic geometry gives a similar impact parameter to within a factor of 2.
- Cosmological expansion (proper vs comoving distance over 40 Mpc) is a few-percent correction at this redshift and doesn't affect the conclusion.

## Generated plots

- `plots/F5c_gw170817_arrival_delay.png`
- `plots/F5c_A_along_path.png`
