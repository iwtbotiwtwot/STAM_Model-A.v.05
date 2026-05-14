# G60 - Pair Structure of Hawking Emission from the Elevator Argument

**Date: 2026-05-13.** Derives the pair structure of Hawking emission from the framework's no-interior + substance ontology + mass-conservation + two-face commitments, and verifies numerical self-consistency with G59's `alpha = 4` and the resolution rule for Hawking T.

Closes the last open piece in G18's entropy derivation chain (P5: 'each Hawking emission is paired -- stated as structural consequence, not derived').

## The elevator argument

When someone steps off an elevator, the elevator gets lighter. These are not two separate events -- they are **one event viewed from two structural sides**: (a) the person's exit, (b) the elevator's mass reduction.

The same holds for an emission event at A = 1:

- **(a) Outward write**: substance escapes the boundary outward. Forced by no-interior (inward writes across A = 1 are forbidden).

- **(b) Horizon reduction**: source mass M decreases; horizon area A_h = 16 pi G^2 M^2 / c^4 shrinks. Forced by mass conservation.



These are not two events. They are **one event with two structural aspects**:

- The outward write is recorded on the **outer face** (where accreted substance lives -- two-face refinement).

- The mass reduction is echoed on the **inner face** (where primordial mass lives) as a reduction in the boundary's reach.



Each emission event therefore necessarily has **two-face structural presence**, derived (not postulated) from the elevator identity. The factor of 2 from two-face in `alpha = 4` (G59) follows.

## What this changes from G18/G59

G18 used 'each Hawking emission is paired' as P5 -- a stated structural commitment, not derived. Under the elevator argument, P5 is **derived** from:

- No-interior at A = 1 (forces outward direction)

- Substance ontology (substance can leave but cannot vanish)

- Mass conservation (Delta M = -Delta E / c^2)

- Two-face refinement (boundary has inner + outer face)

- Single-event identity (the write IS the reduction)



So the entire chain from primitives to S = A_h/(4 L_P^2) is now structurally derived: substance + presentism + ledger-channel + two-face + SU shell-count + elevator -> alpha = 4 -> entropy.

## Numerical self-consistency (verified)

Under the elevator framing, each emission event satisfies:

```
Delta E = k_B T              (one quantum's energy at boundary T)
Delta M = -Delta E / c^2     (mass conservation)
Delta A_h = -alpha L_P^2     (one entry = 4 cells)
```

**Independent derivation of T from elevator + alpha = 4:**

```
4 L_P^2 = (32 pi G^2 M / c^4) (k_B T / c^2)
k_B T  = 4 L_P^2 c^6 / (32 pi G^2 M)
      = (hbar G / c^3) c^6 / (8 pi G^2 M)
      = hbar c^3 / (8 pi G M)   <- standard Hawking k_B T
```

Two independent routes give Hawking T:

- Resolution rule (Q8): `k_B T = (1/4 pi) hbar c |grad A|` at boundary

- Elevator + alpha = 4: `k_B T = 4 L_P^2 c^6 / (32 pi G^2 M)`



Both reduce to `T = hbar c^3 / (8 pi k_B G M)`. The framework is **internally self-consistent**: T derived from boundary substance gradient agrees with T derived from per-entry energy-area-mass relations. Numerically verified for stellar BHs through PBHs to machine precision.

## Emission rate consistency

Framework derives T and pair structure; **absolute rate still requires Stefan-Boltzmann** (boundary radiates as blackbody at T). Under that assumption, framework's `dN_entries/dt = -dA_h/dt / (alpha L_P^2)` matches standard `dN_quanta/dt = P / (k_B T)` to machine precision.

## Status

**Derived in G60 (was P5 'stated' in G18):**

- Pair structure of Hawking emission (elevator argument)

- Two-face contribution to alpha = 4 (single-event with dual aspects)



**Self-consistent chain verified:**

- resolution rule -> Hawking T

- elevator -> pair structure -> alpha = 4 (G59)

- per-entry relations all mutually consistent

- emission rate matches Stefan-Boltzmann blackbody at T



**Still open:**

- Stefan-Boltzmann luminosity of the boundary still assumed; framework gets T and per-entry structure but not absolute rate.

- Spectral details (frequency distribution, grey-body factors) still QFT machinery, not framework-internal.

## Bottom line

G18's 'good suspects, not derivation' caveat for the entropy chain is now **fully resolved**. Pair structure (P5) is derived in G60, alpha = 4 is derived in G59, S = A_h/(4 L_P^2) follows. The framework's structural completeness for boundary thermodynamics now matches its structural completeness for the strong-field metric: both arcs run from the same primitives to derived results without target-shaped decompositions.

## Files

- [scripts/G60_pair_structure_from_elevator.py](../scripts/G60_pair_structure_from_elevator.py)

- [results/G59_entropy_under_ledger_channel_summary.md](G59_entropy_under_ledger_channel_summary.md)

- [results/G18_entropy_from_ledger_counting_summary.md](G18_entropy_from_ledger_counting_summary.md) (historical, with P5 caveat)
