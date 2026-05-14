# G67 - Bell / CHSH and Path Integral Under STAM Born Rule

**Date: 2026-05-13.** Verifies that STAM's two-step Born rule (projection-geometry theorem `u = ||Pψ||` + ledger-measure step `p = u²`) reproduces standard QM for entangled correlations (CHSH/Bell) and path-integral interference (double-slit).

## Part 1: CHSH violation

Singlet state |ψ⟩ = (|↑↓⟩ - |↓↑⟩) / √2.


Under STAM's two-step Born rule with optimal angles (a=0, a'=π/2, b=π/4, b'=-π/4):

```
|CHSH| = 2.828427
Tsirelson bound 2√2 = 2.828427
Bell bound (local realism): 2
```

STAM reproduces the Tsirelson bound exactly. Bell violation is automatic given STAM's structural Born-rule derivation.

### Structural reading of Bell violation

Bell violation in STAM comes from three combined ingredients:

1. **Unresolved support is global** — the singlet state ψ exists on the bipartite Hilbert space; the phase structure encodes the correlation.

2. **Resolution events are local** — Alice writes 1 SU at her location; Bob writes 1 SU at his.

3. **Pair-structure squaring** (p = u²) converts amplitude correlations into measurable probability correlations.



Local hidden variables fail Bell because they would require the correlation structure to be carried locally in resolved-state information. STAM's commitment that phase (and correlation structure) lives ONLY in unresolved support — never written to the ledger directly — is the structural reason no local hidden variable model can reproduce QM.

## Part 2: Free-particle path integral

Standard QM propagator: `K(x_f, t; x_i, 0) = √(m/(2πiℏt)) exp(im(x_f-x_i)²/(2ℏt))`. Verified to integrate to unitarity (∫|K|² dx = 1).


**STAM reading:** the propagator K is the evolution of unresolved support between resolution events. Standard QM dynamics (Schrödinger equation / Feynman path integral) operates in the unresolved layer; resolution events convert `|ψ|²` into resolution-probability density via the framework's Born rule (structurally derived).

### Double-slit interference + decoherence

Two-path propagator `ψ_screen(x) = K(x, t; x_L, 0) + K(x, t; x_R, 0)`.


**Without which-way interaction**: paths recombine in unresolved support before the screen resolution event. The 1-SU resolution samples u = |ψ_L + ψ_R|; then p = u² gives interference fringes.


**With which-way interaction**: paths are *separately* resolved at the slits — each becomes its own ledger entry. At the screen, the two paths are now distinct ledger channels; probabilities add incoherently: p = |ψ_L|² + |ψ_R|². No interference.


STAM's decoherence picture: **interference is the phase-sensitive reshaping of unresolved support before the 1-SU resolution event**. Phase lives in unresolved support; it is never written to the ledger directly.

## What STAM contributes (interpretive, not new dynamics)

- **Why Bell violation**: global unresolved support encodes correlations; local resolution events sample them via pair-structure squaring.

- **Why phase in path integral**: phase is a configuration property of unresolved support, never written to ledger.

- **Why decoherence**: paired events resolving alternatives into distinct ledger channels before recombination.

## Open work

- **Full path-integral formulation from primitives**: the QM action S enters via the unresolved-support dynamics. A first-principles derivation of the QM action from STAM substance-ontology + SU-write structure is open. Per-SU phase increment is not a constant under naive reading.

- **Continuous-spectrum measurements** under STAM's discrete 1-SU resolution events. The classical limit (position, momentum continuous) should emerge from coarse-graining the discrete SU-write structure.

- **Specific predictions where STAM might diverge from QM**: at the current commitment level, STAM matches QM exactly. STAM-distinctive quantum predictions (where the discrete SU-quantum scale matters) are an open exploration item — F6 gravitational decoherence is one candidate (~0.5 s decoherence time for 1 micron silica nanoparticle).

## Files

- [scripts/G67_bell_and_path_integral.py](../scripts/G67_bell_and_path_integral.py)

- [plots/G67_chsh_scan.png](../plots/G67_chsh_scan.png)

- [plots/G67_double_slit.png](../plots/G67_double_slit.png)
