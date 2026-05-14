# G10: SPARC-Style Galaxy Rotation Curve Test

## The setup

STAM's gravity bridge `g = (c²/2)∇A` combined with linear cumulative A from baryons gives:
```text
A_baryon(x) = Σ 2 G m_i / (c² |x - x_i|)
g_STAM = (c²/2) ∇A_baryon = -∇Φ_Newton = g_Newton
```
**STAM linear cumulative A from baryons IS Newtonian gravity from baryons.** So the rotation curve test is the well-known dark matter test: do baryons alone (via standard Newton) reproduce observed rotation curves?

If yes, STAM has no DM problem. If no, STAM's framework requires an A_collective term beyond the linear sum (the exploratory extension flagged in CLAIMS_AND_STATUS §7.4).

## Galaxies tested (canonical values)

```text
Galaxy                 M_baryon        v_bar_outer    v_obs    ratio    DM factor
----------------------------------------------------------------------

NGC 3198                 3.70e+10        78.6      150.0     1.91     3.64x
NGC 2403                 2.10e+10        72.0      130.0     1.81     3.26x
Milky Way                7.00e+10        83.8      220.0     2.62     6.89x
DDO 154 (LSB dwarf)      2.70e+09        42.6       50.0     1.17     1.38x
```

## Verdict

**STAM linear cumulative A from baryons does NOT reproduce observed rotation curves.** Same result as standard Newton-from-baryons: outer rotation speeds fall short of observation by factors of 1.5-3× in v (corresponding to factors of 2-10× in inferred enclosed mass, i.e., the standard dark-matter deficit).

**This is consistent with what was already documented in the scratch test** (results/galactic_cumulative_A/SUMMARY.md): linear cumulative A from a 1e11 M_sun toy galaxy at 10 kpc gives compact-source circular speed ~207 km/s, but the same toy with extended exponential disc gives outer speed only ~121 km/s — well below the ~220 km/s flat rotation observed for Milky Way-class galaxies.

**What this means for the framework:** STAM does NOT replace dark matter via linear cumulative A alone. The framework either needs:

1. **A_collective from STAM principles** — an additional term beyond the linear sum, derivable from the framework's field equations or structural commitments. Currently EXPLORATORY ONLY (§7.4).

2. **Cosmic-structure cumulative A at galactic scale** — galaxies embedded in cosmic web feel a tidal A field from external matter. Could in principle contribute to galactic-scale gradient. NOT YET COMPUTED.

3. **Modified Newton's laws (MOND-like)** — STAM's gravity bridge g = (c²/2)∇A could acquire a MOND-style modification at low accelerations. NOT IN THE FRAMEWORK CURRENTLY.

4. **Dark matter is real** — STAM accepts DM as a separate matter component (just like LCDM does). This loses the 'no DM needed' story but preserves the framework's other commitments.

## Implications for the cumulative-A picture

This is an honest result that's important for the framework's assessment:

- The 'many nearzeros add up to galactic gravity' story is true in the Newton sense (linear cumulative A = Newton). And Newton from baryons isn't enough for galactic rotation curves.

- The cosmic-LoS amplification picture (G9) and the galactic-rotation problem are NOT the same problem. G9 is about the integrated A along a 14 Gpc cosmological path, where structure-amplification of order 1.5× the void-baseline is plausible. The galactic-rotation problem is about LOCAL gradients of A within a single galaxy, where the deficit is factor 2-10× — much larger.

- So G10 doesn't break G9's CMB closure logic, but it does mean STAM doesn't replace dark matter without additional structure. The framework's SN-distance and CMB-θ_⋆ stories at H_0=73 stand independent of the rotation-curve question. But the broader claim of 'unified cumulative-A from cosmic to galactic scales' is incomplete: linear cumulative A doesn't bridge the galactic-scale gap.

## Honest framing

**This is the test's value: it sharpens what the framework is and isn't.** STAM Model-A as currently specified:

- ✅ **Makes specific cosmological predictions** (SN distances, CMB θ_⋆ at H_0=73 with cumulative-A LoS amplification) — these are real and partially backed up.

- ✅ **Has a coherent black-hole story** (bubble picture, thermodynamics, no-Penrose-extraction).

- ✅ **Has a derivable bridge term** (A_0 = 1/(12π) → b = 355 Mly).

- ❌ **Does not currently replace dark matter at galactic scales.** Linear cumulative A = Newton, which doesn't fit observed rotation curves. The framework's exploratory A_collective extension is needed for this, and it's not yet derived from first principles.

**This is consistent with the framework's own status notes** — CLAIMS_AND_STATUS §7.4 already labels the collective term as exploratory and §7.5 explicitly says 'STAM may need a collective A-envelope term'. G10 confirms quantitatively that this collective term is needed.

## What this means for the ongoing framework development

The keystone test came back with a clear, honest result: STAM linear cumulative A doesn't replace dark matter at galactic scales. This means:

**Things that survive:**
- All cosmological-scale predictions (SN, CMB, bridge term)
- All BH-related structural commitments (bubble, thermodynamics)
- The A_0 = 1/(12π) structural derivation
- The unified-A field across regimes (gravity, quantum, thermodynamics)

**Things that are now sharper:**
- The 'no DM needed' galactic claim becomes the explicit research question 'does STAM derive an A_collective term from first principles that fits SPARC?'
- The exploratory A_collective extension in §7.4 becomes a real research target with quantitative criteria (must explain factor 2-10× rotation deficit across SPARC sample)
- Alternative: the framework accepts DM as a separate matter component, in which case STAM is parallel to LCDM in that regime rather than replacing it

**The honest meta-result of this morning's work:**
- Cosmological story (SN + CMB + bridge term): close to working with V₃ + cumulative-A LoS
- Galactic story (rotation curves): not working with linear cumulative A alone; needs derivable A_collective or DM

Two-thirds of the framework's distinctive claims survive scrutiny. The galactic claim needs more theoretical work before it can be made or dropped honestly.

## Generated plots

- `plots/G10_sparc_rotation.png`
