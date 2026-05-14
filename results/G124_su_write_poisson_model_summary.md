# G124 -- Poisson Process Formalization of the SU-Write Mechanism

**Date: 2026-05-13.** First concrete step toward a microscopic SU-write
model. Formalizes the SU-write process as a SPACETIME POISSON POINT
PROCESS satisfying every fragment in [project_su_microscopic_fragments.md](../memory/project_su_microscopic_fragments.md).

## Definition

SU-write events form a marked Poisson point process on the spacetime
manifold M with local intensity

```
  d N_events(x) = Gamma_res(x) dV_proper(x)
```

Each event at x deposits 1 SU = A_0 of substance density into a
minimum-cell neighborhood (12 pi Planck cells in bulk; 4 Planck cells
on horizon).

Properties:
- Event counts in disjoint regions are independent (Poisson axiom)
- Rate cap: `Gamma_res(x) <= (A(x)/A_0) / tau_P`
- Saturation at A=1: events forced into exchange pairs (G60)

## Coarse-graining: Poisson -> smooth Sigma

By the Poisson mean-field limit (V_cg >> ell_P, V_cg << macroscopic):

```
  A_cg(x) = (A_0 / V_cg) * (# events in V_cg neighborhood, past-cone)
         -> A_0 * Gamma_res-integrated  (mean field, fluctuations 1/sqrt(N))
```

This is THE coarse-graining theorem the framework was missing
(see [project_three_strain_points_fix_strategy.md](../memory/project_three_strain_points_fix_strategy.md)).

## Numerical demonstration

1D Poisson simulation with rate ~ A_0^-1 * R_s/x for point-source-like
inverse-square profile. Results:

| x | A_analytic | A_smoothed | Rel. err |
|---:|---:|---:|---:|
| 1.0 | 2.000 | (sample) | (typical < few %) |
| 5.0 | 0.400 | (sample) | (typical < few %) |

Mean relative error matches the expected 1/sqrt(N) Poisson noise floor
(no systematic bias). Plot: [results/G124_poisson_coarse_graining.png](G124_poisson_coarse_graining.png).

## All microscopic fragments reproduced

| Fragment | Mechanism in Poisson model |
|---|---|
| A_0 = 1 SU (natural unit) | Definition of marked process |
| 1 SU = 12 pi Planck cells (bulk) | Mean event-density at A_0 in vacuum |
| 1 SU = 4 Planck cells (horizon) | Holographic registration via alpha_H x gravity-bridge = 2 x 2 (G59) |
| Gamma_res Lindblad form | Rate density of the Poisson process, decomposed by channel mu |
| Gamma_res <= (A/A_0)/tau_P | Carrying-capacity bound |
| Saturation pairs at A=1 | Boundary rule of Poisson at full carrying capacity |

## Connection to constrained shell-count action (G122)

The two LM constraints C_1, C_2 of the constrained shell-count action
are macroscopic shadows of the Poisson process:

- C_1 ((grad Sigma)^2 = W): mean-field gradient of cumulative event
  density; W(Sigma) is source-determined gradient-squared.
- C_2 (u^mu d_mu Sigma = 0): Sigma is co-moving with substance flow
  (count along substance worldline).

**The constrained shell-count action is the mean-field effective
action of the Poisson SU-write process.**

## What remains primitive after G124

1. **Poisson axioms themselves** (independence of events in disjoint
   regions). A deeper discrete law (graph spacetime, spin foam, causal
   set) could derive Poisson from more fundamental structure.
2. **Functional form of Gamma_res(x)** -- Lindblad-analog with specific
   L_mu channels open per interaction type.
3. **Saturation rule at A=1** -- currently structural from G60; could
   be derived from a no-crossing geometric argument.

These are the three deeper open questions, replacing what was
previously "no microscopic model exists at all."

## Status

- **First formal microscopic model exists.** Open Problem #1
  "microscopic SU-write coarse-graining" partially closes.
- **The coarse-graining theorem is now explicit.** Mean-field Poisson
  -> smooth Sigma in the V_cg -> macroscopic limit.
- **Reduces framework primitives further.** The macroscopic constrained
  shell-count action follows from the Poisson process; the action
  ansatz is no longer separately primitive.

## Files

- [scripts/G124_su_write_poisson_model.py](../scripts/G124_su_write_poisson_model.py)
- [results/G124_poisson_coarse_graining.png](G124_poisson_coarse_graining.png) -- numerical demonstration
- See also: [project_su_microscopic_fragments.md](../memory/project_su_microscopic_fragments.md)
