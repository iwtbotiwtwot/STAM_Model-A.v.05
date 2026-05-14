"""
G43: Doughnut Position and Contraction — Full STAM Correction
==============================================================

Same single question as G42, but with the FULL self-consistent A_local
at the orbital position (static + kinematic, both substance-corrected).

Picture:
  1. Inspiral builds doughnut to position d at orbit-stop time
  2. After orbit-stop, no more orbital driving → doughnut contracts inward
  3. Contraction at STAM substance speed (the wave speed in the medium)
  4. Eventually meets the explosion shell from the merger event

All in STAM time and speed (no naive Newtonian or c assumptions).
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
print("G43: Doughnut position + contraction — full STAM correction")
print("="*78)
print(f"\n  M_total = 2.7 M_sun, R_s = {R_s/1000:.4f} km")
print(f"  A_0 = 1/(12π) = {A_0:.6f}")

# --- Step 1: Self-consistent threshold (with static + kinematic A) ---
# At orbital separation r, each NS sees:
#   A_static = R_s_other/r = R_s/(2r)
#   A_kinematic_naive = R_s/(2r)
# With STAM correction:
#   v_actual² = v_Newton² × (1-A_local)
#   A_kinematic_corrected = (R_s/2r) × (1 - A_static - A_kinematic)
#   Let x = R_s/(2r). Then:
#     A_kin = x(1 - x - A_kin)
#     A_kin = x(1-x) / (1+x)
# At threshold A_kin = A_0:
#   A_0 = x(1-x)/(1+x)
#   x² + x(A_0 - 1) + A_0 = 0
#   x = [1-A_0 - √((1-A_0)² - 4A_0)] / 2

print("\n" + "-"*78)
print("Step 1: Self-consistent threshold (full STAM A)")
print("-"*78)

discriminant = (1 - A_0)**2 - 4*A_0
x_threshold = (1 - A_0 - np.sqrt(discriminant)) / 2
r_threshold = R_s / (2 * x_threshold)

print(f"  Solving x² + x(A_0 - 1) + A_0 = 0 for the threshold radius:")
print(f"    x = R_s/(2r) = {x_threshold:.6f}")
print(f"    r_threshold = R_s/(2x) = {r_threshold/1000:.4f} km")
print(f"    r_threshold/R_s = {r_threshold/R_s:.4f}")

print(f"\n  Compare:")
print(f"    Naive (no STAM):              r/R_s = 6π = {6*np.pi:.4f}")
print(f"    Partial STAM (kinematic only): r/R_s = 6π·(1-A_0) = {6*np.pi*(1-A_0):.4f}")
print(f"    Full STAM (this calc):         r/R_s = {r_threshold/R_s:.4f}")

# --- Step 2: τ_critical (Peters-Mathews chirp from threshold) ---
print("\n" + "-"*78)
print("Step 2: τ_critical = chirp time from threshold to merger")
print("-"*78)

R_s_over_c = R_s / c
tau_critical = (5/8) * (r_threshold/R_s)**4 * R_s_over_c

tau_naive = (5/8) * (6*np.pi)**4 * R_s_over_c
tau_partial = (5/8) * (6*np.pi*(1-A_0))**4 * R_s_over_c

print(f"  τ_critical (full STAM):     {tau_critical:.4f} s")
print(f"  Compared to:")
print(f"    Naive:        {tau_naive:.4f} s")
print(f"    Partial STAM: {tau_partial:.4f} s")

# --- Step 3: Doughnut location at orbit-stop ---
print("\n" + "-"*78)
print("Step 3: Doughnut location at orbit-stop (wave propagation through substance)")
print("-"*78)
print(f"  Wave coordinate speed through baseline A_0: v_substance = c × (1-A_0)² × (1+A_0)")

v_sub_factor = (1 - A_0)**2 * (1 + A_0)
v_sub = c * v_sub_factor
print(f"    Speed factor: {v_sub_factor:.6f}")
print(f"    v_substance ≈ {v_sub_factor:.3f} c")

d_doughnut_c = c * tau_critical  # if waves propagated at c
d_doughnut_sub = v_sub * tau_critical  # at substance speed

print(f"\n  Doughnut outer edge at orbit-stop:")
print(f"    With c-propagation:        d = c × τ = {tau_critical:.4f} light-seconds")
print(f"    With substance-propagation: d = {v_sub_factor:.4f}c × τ = {tau_critical*v_sub_factor:.4f} light-seconds")
print(f"                                  = {d_doughnut_sub:.3e} m")

# --- Step 4: Contraction back toward the event ---
print("\n" + "-"*78)
print("Step 4: Contraction time (all in STAM time and speed)")
print("-"*78)
print(f"  After orbit-stop, substance is no longer driven.")
print(f"  Doughnut contracts inward at v_substance toward the merger.")
print(f"  Time to reach r=0: τ_contract = d_doughnut / v_substance")

tau_contract_c = d_doughnut_c / c  # If we use c for both location and speed
tau_contract_sub = d_doughnut_sub / v_sub  # Both use substance speed

print(f"\n  If c-based throughout: τ_contract = {tau_critical*1:.4f} s = τ_critical")
print(f"  If substance-based throughout: τ_contract = {tau_contract_sub:.4f} s")
print(f"    (The substance factors cancel — both location and speed scale together.)")

# --- Step 5: Compare to observation ---
print("\n" + "="*78)
print("Comparison to GW170817")
print("="*78)

t_observed = 1.74

print(f"\n  Observed engine time (orbit-stop → light arrival): 1.74 s")
print(f"\n  Prediction chain:")
print(f"    Naive (Newton, no STAM):        {tau_naive:.4f} s  ({(tau_naive-t_observed)/t_observed*100:+.1f}%)")
print(f"    Partial STAM (kinematic only):  {tau_partial:.4f} s  ({(tau_partial-t_observed)/t_observed*100:+.1f}%)")
print(f"    Full STAM (static + kinematic): {tau_critical:.4f} s  ({(tau_critical-t_observed)/t_observed*100:+.1f}%)")
print(f"")
print(f"  ⇒ Doughnut at orbit-stop is at {tau_critical*v_sub_factor:.3f} light-seconds out (substance speed).")
print(f"  ⇒ Contraction takes the same τ_critical seconds (factors cancel).")
print(f"  ⇒ Total engine time (if light emitted at contraction completion) = {tau_critical:.3f} s")
print(f"")
print(f"  Gap to observed 1.74 s: {(tau_critical - t_observed)/t_observed*100:+.1f}%")
print(f"  (~ {(tau_critical - t_observed)*1000:.0f} ms)")
