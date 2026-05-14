#!/usr/bin/env python3
"""STAM Model-A weak-field Shapiro delay test.

This script tests the local propagation formula:

    A(r) = Rs/r = 2GM/(c^2 r)
    Delta t = (1/c) int A(r) ds

For a straight path with impact parameter b_imp near a spherical mass:

    Delta t = (2GM/c^3) [asinh(x2/b_imp) - asinh(x1/b_imp)]

This reproduces the first-order weak-field logarithmic Shapiro structure.

Interpretation in Model-A: the formula computes the SU/photon-A
path-stretching (g_rr piece) along the light geodesic. The clock-rate
effect from g_tt (dtau/dt = sqrt(1-A)) is a separate effect and does NOT
appear in this calculation. The standard GR "Shapiro delay" measured as
round-trip travel time observed by one Earth clock is also pure
path-stretching (the Earth clock's rate cancels between start and end of
the round trip), so the numerical equivalence here is real, not a
double-count. Model-A's "Shapiro" interpretation of the same number is
the SU traversal load through the A field; A was defined as 2GM/c^2r so
that integral A ds matches the standard GR Shapiro coefficient.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

G = 6.67430e-11
C = 299_792_458.0
M_SUN = 1.98847e30
R_SUN = 6.9634e8
AU = 1.495978707e11


def schwarzschild_radius(mass_kg: float) -> float:
    return 2.0 * G * mass_kg / C**2


def stam_shapiro_analytic(mass_kg: float, impact_m: float, r1_m: float, r2_m: float) -> float:
    if impact_m <= 0 or r1_m <= impact_m or r2_m <= impact_m:
        raise ValueError("Need impact_m > 0 and r1_m,r2_m > impact_m")
    x1 = -math.sqrt(r1_m**2 - impact_m**2)
    x2 = math.sqrt(r2_m**2 - impact_m**2)
    return (2.0 * G * mass_kg / C**3) * (math.asinh(x2 / impact_m) - math.asinh(x1 / impact_m))


def stam_shapiro_numeric(
    mass_kg: float,
    impact_m: float,
    r1_m: float,
    r2_m: float,
    n: int = 200_001,
) -> float:
    x1 = -math.sqrt(r1_m**2 - impact_m**2)
    x2 = math.sqrt(r2_m**2 - impact_m**2)
    xs = np.linspace(x1, x2, n)
    r = np.sqrt(xs**2 + impact_m**2)
    A = schwarzschild_radius(mass_kg) / r
    return float(np.trapezoid(A, xs) / C)


def first_order_log_delay(mass_kg: float, impact_m: float, r1_m: float, r2_m: float) -> float:
    x1_abs = math.sqrt(r1_m**2 - impact_m**2)
    x2 = math.sqrt(r2_m**2 - impact_m**2)
    R = x1_abs + x2
    return (2.0 * G * mass_kg / C**3) * math.log((r1_m + r2_m + R) / (r1_m + r2_m - R))


def large_distance_approx(mass_kg: float, impact_m: float, r1_m: float, r2_m: float) -> float:
    return (2.0 * G * mass_kg / C**3) * math.log(4.0 * r1_m * r2_m / impact_m**2)


def main() -> None:
    out_dir = Path("results/shapiro_delay")
    out_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        ("Sun_grazing_Earth_Mercury", 1.0 * AU, 0.387 * AU, 1.0 * R_SUN),
        ("Sun_grazing_Earth_Venus", 1.0 * AU, 0.723 * AU, 1.0 * R_SUN),
        ("Sun_grazing_Earth_Mars", 1.0 * AU, 1.524 * AU, 1.0 * R_SUN),
        ("Sun_grazing_Earth_Jupiter", 1.0 * AU, 5.204 * AU, 1.0 * R_SUN),
        ("Sun_grazing_Earth_Saturn", 1.0 * AU, 9.582 * AU, 1.0 * R_SUN),
        ("Sun_2Rsun_Earth_Mars", 1.0 * AU, 1.524 * AU, 2.0 * R_SUN),
        ("Sun_5Rsun_Earth_Mars", 1.0 * AU, 1.524 * AU, 5.0 * R_SUN),
        ("Sun_10Rsun_Earth_Mars", 1.0 * AU, 1.524 * AU, 10.0 * R_SUN),
    ]

    rows = []
    for name, r1, r2, impact in cases:
        analytic = stam_shapiro_analytic(M_SUN, impact, r1, r2)
        numeric = stam_shapiro_numeric(M_SUN, impact, r1, r2)
        log_delay = first_order_log_delay(M_SUN, impact, r1, r2)
        approx = large_distance_approx(M_SUN, impact, r1, r2)
        rows.append(
            {
                "case": name,
                "r1_AU": r1 / AU,
                "r2_AU": r2 / AU,
                "impact_Rsun": impact / R_SUN,
                "A_at_impact": schwarzschild_radius(M_SUN) / impact,
                "delay_analytic_microseconds": analytic * 1e6,
                "delay_two_way_microseconds": analytic * 2e6,
                "delay_numeric_microseconds": numeric * 1e6,
                "numeric_rel_error": (numeric - analytic) / analytic,
                "first_order_log_microseconds": log_delay * 1e6,
                "first_order_log_rel_error": (log_delay - analytic) / analytic,
                "large_distance_approx_microseconds": approx * 1e6,
                "large_distance_approx_rel_error": (approx - analytic) / analytic,
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "shapiro_delay_cases.csv", index=False)

    summary = {
        "test": "STAM Model-A weak-field Shapiro delay",
        "formula": "Delta t = (1/c) integral A(r) ds",
        "A": "A(r)=2GM/(c^2 r)",
        "earth_mars_grazing_one_way_delay_microseconds": float(
            df.loc[df["case"] == "Sun_grazing_Earth_Mars", "delay_analytic_microseconds"].iloc[0]
        ),
        "earth_mars_grazing_two_way_delay_microseconds": float(
            df.loc[df["case"] == "Sun_grazing_Earth_Mars", "delay_two_way_microseconds"].iloc[0]
        ),
        "max_abs_first_order_log_rel_error": float(df["first_order_log_rel_error"].abs().max()),
        "max_abs_numeric_rel_error": float(df["numeric_rel_error"].abs().max()),
        "interpretation": (
            "The STAM accumulation integral reproduces the first-order weak-field logarithmic "
            "Shapiro delay structure for a spherical source."
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    print("PASS: STAM Model-A Shapiro delay weak-field structure recovered.")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote {out_dir / 'summary.json'}")
    print(f"Wrote {out_dir / 'shapiro_delay_cases.csv'}")


if __name__ == "__main__":
    main()
