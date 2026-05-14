"""
G42: Doughnut Threshold with STAM-Corrected Orbital Speed
==========================================================

Strip the calculation down to ONE question:

  If the orbiting binary moves slower than v_Newton because the substance
  resists motion through elevated A, where does the doughnut threshold sit?

The threshold is the orbital separation r at which the binary's kinematic
A first exceeds the cosmic baseline A_0.

Naive (no STAM correction):
  v² = GM/r (Keplerian, no substance effect)
  A_kinematic = v²/c² = R_s/(2r)
  At threshold: A_kinematic = A_0  →  r = R_s/(2·A_0) = R_s × 6π

STAM-corrected (substance slows motion):
  v_actual² = v_Newton² × (1-A_local)
  A_kinematic_corrected = (R_s/(2r)) × (1-A)   (self-consistent)
  Solving:  A = R_s / (2r + R_s)
  At threshold: A = A_0  →  r = R_s × (1-A_0) / (2·A_0)

Then convert to doughnut spatial location at merger time:
  r_doughnut_shell = c · τ_critical
  where τ_critical is Peters-Mathews chirp time from threshold to merger.
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

# Constants
c = 299792458.0
G = 6.67430e-11
M_sun = 1.989e30

# STAM
A_0 = 1.0 / (12 * np.pi)

# GW170817
M_total = 2.7 * M_sun
R_s = 2 * G * M_total / c**2

print("="*78)
print("G42: Doughnut threshold with STAM-corrected orbital speed")
print("="*78)
print(f"\n  M_total = 2.7 M_sun, R_s = {R_s/1000:.4f} km")
print(f"  A_0 = 1/(12π) = {A_0:.6f}")

# Step 1: Naive threshold (no correction)
print("\n" + "-"*78)
print("Step 1: NAIVE threshold (Keplerian, no STAM correction)")
print("-"*78)
r_naive = R_s / (2 * A_0)  # = R_s × 6π
print(f"  v² = GM/r = R_s c²/(2r)")
print(f"  A_kinematic = v²/c² = R_s/(2r)")
print(f"  At threshold A_kinematic = A_0:")
print(f"    r_threshold_naive = R_s/(2A_0) = R_s × 6π = {r_naive/1000:.4f} km")
print(f"    r_naive/R_s = {r_naive/R_s:.4f}")

# Step 2: STAM-corrected threshold
print("\n" + "-"*78)
print("Step 2: STAM-CORRECTED threshold (substance-velocity coupling)")
print("-"*78)
print(f"  v_actual² = v_Newton² × (1-A)")
print(f"  A_kinematic_corrected = (R_s/2r) × (1-A)")
print(f"  Self-consistent: A = R_s/(2r + R_s)")
print(f"  At threshold A = A_0:")
print(f"    r_threshold_corrected = R_s × (1-A_0) / (2·A_0) = R_s × 6π × (1-A_0)")

r_corrected = R_s * (1 - A_0) / (2 * A_0)
print(f"    = {r_corrected/1000:.4f} km")
print(f"    r_corrected/R_s = {r_corrected/R_s:.4f}")

print(f"\n  Shift relative to naive:")
print(f"    Δr/r_naive = -{(r_naive - r_corrected)/r_naive * 100:.2f}%")
print(f"    (threshold moves INWARD)")

# Step 3: Peters-Mathews chirp time from each threshold
print("\n" + "-"*78)
print("Step 3: Time from threshold to merger (Peters-Mathews)")
print("-"*78)
print(f"  t_chirp(r) = (5/8) × (r/R_s)⁴ × R_s/c    (equal-mass binary)")

R_s_over_c = R_s / c
print(f"  R_s/c = {R_s_over_c:.4e} s for this mass")

tau_naive = (5/8) * (r_naive/R_s)**4 * R_s_over_c
tau_corrected = (5/8) * (r_corrected/R_s)**4 * R_s_over_c

print(f"\n  τ_critical_NAIVE       = (5/8) × {(r_naive/R_s):.4f}⁴ × {R_s_over_c:.3e}")
print(f"                          = {tau_naive:.4f} s")
print(f"\n  τ_critical_CORRECTED   = (5/8) × {(r_corrected/R_s):.4f}⁴ × {R_s_over_c:.3e}")
print(f"                          = {tau_corrected:.4f} s")

print(f"\n  Ratio corrected/naive = {tau_corrected/tau_naive:.4f}")
print(f"  Reduction: {(tau_naive - tau_corrected)/tau_naive * 100:.2f}%")

# Step 4: Doughnut shell at merger time
print("\n" + "-"*78)
print("Step 4: Doughnut shell location at merger time")
print("-"*78)
print(f"  Shell radius at merger = c · τ_critical")

r_shell_naive = c * tau_naive
r_shell_corrected = c * tau_corrected

print(f"\n  r_shell_NAIVE     = c × {tau_naive:.4f} = {r_shell_naive:.3e} m")
print(f"                    = {tau_naive:.4f} light-seconds")
print(f"\n  r_shell_CORRECTED = c × {tau_corrected:.4f} = {r_shell_corrected:.3e} m")
print(f"                    = {tau_corrected:.4f} light-seconds")

# Step 5: Compare to observation
print("\n" + "="*78)
print("Comparison to GW170817 observation")
print("="*78)

t_observed = 1.74

print(f"\n  Observed GW-EM gap (or bubble retraction time): 1.74 s")
print(f"\n  NAIVE prediction:    {tau_naive:.4f} s     ({(tau_naive - t_observed)/t_observed * 100:+.1f}%)")
print(f"  CORRECTED prediction: {tau_corrected:.4f} s     ({(tau_corrected - t_observed)/t_observed * 100:+.1f}%)")

print(f"\n  The STAM correction moves the prediction CLOSER to observed by")
print(f"  shifting the threshold inward (substance can't be moved as fast as Newton).")

# Show the doughnut shell location side-by-side with observed
print(f"\n  Doughnut shell location at merger (substance-corrected):")
print(f"    r_shell = {r_shell_corrected/1e6:.2f} × 10⁶ km")
print(f"           = {tau_corrected:.4f} light-seconds")
print(f"")
print(f"  Observed GW-light gap = {t_observed} light-seconds")
print(f"  Discrepancy: {(tau_corrected - t_observed)/t_observed * 100:+.1f}%")
