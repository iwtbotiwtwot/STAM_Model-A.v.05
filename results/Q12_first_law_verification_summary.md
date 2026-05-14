# STAM Model-A Q12: First Law, Smarr Formula, Evaporation, GSL

## Purpose

Close the bold-STAM thermodynamic sector. Q11 made the differential first law hold by construction; Q12 confirms it numerically and adds three real checks that are *not* by construction:
1. **Smarr formula** `M c^2 = 2 T S` — integrated form of the first law.
2. **Hawking evaporation dynamics** — M(t), T(t), S(t), and total lifetime.
3. **Generalized second law** — total entropy grows by a 1/3 surplus during evaporation.

## Numerical results

|   M_over_Msun |    T_STAM_K |   S_STAM_J_per_K |   first_law_residual |   smarr_residual |   2TS_over_Mc2 |   GSL_surplus |         L_W |   tau_evap_yr |
|--------------:|------------:|-----------------:|---------------------:|-----------------:|---------------:|--------------:|------------:|--------------:|
|         1e-06 | 0.0617007   |      1.44824e+42 |          0           |      2.16466e-16 |              1 |      0.333333 | 9.00761e-17 |   2.09568e+49 |
|         0.001 | 6.17007e-05 |      1.44824e+48 |          1.78024e-16 |      2.21661e-16 |              1 |      0.333333 | 9.00761e-23 |   2.09568e+58 |
|         1     | 6.17007e-08 |      1.44824e+54 |          0           |     -1.1349e-16  |              1 |      0.333333 | 9.00761e-29 |   2.09568e+67 |
|         1e+06 | 6.17007e-14 |      1.44824e+66 |          1.78024e-16 |      1.19003e-16 |              1 |      0.333333 | 9.00761e-41 |   2.09568e+85 |
|         1e+10 | 6.17007e-18 |      1.44824e+74 |         -1.78024e-16 |      0           |              1 |      0.333333 | 9.00761e-49 |   2.09568e+97 |



**Reading the table:**
- `first_law_residual ~ 0`: differential first law `c^2 = T dS/dM` holds to machine precision.
- `smarr_residual ~ 0` and `2TS / (Mc^2) ~ 1`: Smarr formula confirmed.
- `GSL_surplus = 1/3`: total entropy grows by 1/3 of |dS_BH| during evaporation.
- `tau_evap_yr`: total Hawking evaporation lifetime in years (Stefan-Boltzmann photon-only).

## What's tested vs what's verified

- **By construction (passes trivially):** differential first law `c^2 = T dS/dM`. Q11 derived S by integrating exactly this equation, so it must hold. Confirmed.
- **Real check, not by construction (Smarr):** `M c^2 = 2 T S`. This is a non-trivial integrated relation, not the same as the differential first law. It holds for bold STAM exactly. Confirmed.
- **Real check, not by construction (GSL):** `dS_total = dS_BH + dS_rad >= 0`. The BH's entropy decreases during evaporation; the radiation carries entropy `(4/3)(c^2 dM)/T`. Net entropy grows by 1/3 of |dS_BH|. The factor 4/3 is standard thermal-radiation thermodynamics; bold STAM inherits it because it inherits the Planckian spectrum from Q10. Confirmed.
- **Real check, not by construction (negative heat capacity):** as M decreases, T increases. Visible in the evaporation trajectory plot. Confirmed.

## What the evaporation lifetime says

For a solar-mass BH using Stefan-Boltzmann photon-only emission:
- Lifetime: ~2.1e+67 years.
- Standard textbook Hawking-evaporation lifetime for solar-mass BHs is roughly `10^67 years`, vastly longer than the age of the universe. STAM and standard agree on the scaling; the absolute value depends on which species are radiated (photons only vs all massless fields).
- For primordial BHs of ~10^11 to 10^12 kg, the evaporation completes within the age of the universe — these are the targets for any observational tests of Hawking radiation.

## Bold-STAM thermodynamic sector status

The thermodynamic weapons rack is now full. Bold STAM has, with no fitted parameters and no GR-imported machinery:
- Hawking temperature `T = hbar c |grad A| / (4 pi k_B)` — Q8/Q10
- Planckian spectrum at T — Q10
- Bekenstein-Hawking entropy `S = k_B A / (4 ell_P^2)` — Q11
- Differential first law `dE = T dS` — Q11/Q12
- Smarr formula `M c^2 = 2 T S` — Q12
- Hawking evaporation dynamics — Q12
- Generalized second law (1/3 surplus) — Q12
- Negative heat capacity (BH gets hotter as it shrinks) — Q12

Every result matches the standard prediction. Mechanism is bold-STAM-native: phase boundary, bubble, no interior. There is no observational test that would distinguish bold STAM from standard BH thermodynamics; the distinction is interpretive.

## Generated plots

- `plots/Q12_smarr_check.png`
- `plots/Q12_evaporation_trajectory.png`
- `plots/Q12_gsl_check.png`
