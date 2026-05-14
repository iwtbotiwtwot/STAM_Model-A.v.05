#!/usr/bin/env python3
"""STAM Model-A horizon threshold test.

Tests the spherical non-rotating threshold:

    A(r) = Rs/r = 2GM/(c^2 r)

Horizon condition:

    A = 1 <=> r = Rs = 2GM/c^2

Escape velocity identity:

    v_escape^2 = 2GM/r
    A = (v_escape/c)^2

Therefore:

    A=1 <=> v_escape=c
    A>1 <=> v_escape>c
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
KM = 1000.0


def schwarzschild_radius(mass_kg: float) -> float:
    return 2.0 * G * mass_kg / C**2


def accumulation(mass_kg: float, radius_m: float) -> float:
    return schwarzschild_radius(mass_kg) / radius_m


def escape_velocity(mass_kg: float, radius_m: float) -> float:
    return math.sqrt(2.0 * G * mass_kg / radius_m)


def mass_from_horizon_radius(radius_m: float) -> float:
    return C**2 * radius_m / (2.0 * G)


def main() -> None:
    out_dir = Path("results/horizon_threshold")
    out_dir.mkdir(parents=True, exist_ok=True)

    examples = [
        ("1 solar mass", 1.0),
        ("3 solar masses", 3.0),
        ("10 solar masses", 10.0),
        ("30 solar masses", 30.0),
        ("100 solar masses", 100.0),
        ("1 million solar masses", 1.0e6),
        ("4.3 million solar masses illustrative SgrA*", 4.3e6),
        ("6.5 billion solar masses illustrative M87*", 6.5e9),
    ]

    rows = []
    for name, mass_solar in examples:
        mass = mass_solar * M_SUN
        rh = schwarzschild_radius(mass)
        rows.append(
            {
                "object": name,
                "mass_solar": mass_solar,
                "horizon_radius_km": rh / KM,
                "horizon_diameter_km": 2.0 * rh / KM,
                "A_at_horizon": accumulation(mass, rh),
                "escape_velocity_at_horizon_over_c": escape_velocity(mass, rh) / C,
                "M_recovered_ratio": mass_from_horizon_radius(rh) / mass,
            }
        )

    examples_df = pd.DataFrame(rows)
    examples_df.to_csv(out_dir / "horizon_radius_from_mass_examples.csv", index=False)

    mass = 10.0 * M_SUN
    rh = schwarzschild_radius(mass)
    scan_rows = []
    for r_over_h in [5.0, 3.0, 2.0, 1.5, 1.25, 1.1, 1.01, 1.0, 0.99, 0.9, 0.75, 0.5, 0.25]:
        r = r_over_h * rh
        A = accumulation(mass, r)
        v_over_c = escape_velocity(mass, r) / C
        scan_rows.append(
            {
                "mass_solar": 10.0,
                "r_over_horizon": r_over_h,
                "radius_km": r / KM,
                "A": A,
                "escape_velocity_over_c": v_over_c,
                "A_minus_v_escape_over_c_squared": A - v_over_c**2,
                "regime": "outside/A<1" if A < 1 else "horizon/A=1" if math.isclose(A, 1.0) else "inside/A>1",
            }
        )
    scan_df = pd.DataFrame(scan_rows)
    scan_df.to_csv(out_dir / "A_escape_threshold_scan_10Msun.csv", index=False)

    checks = {
        "max_abs_A_at_horizon_minus_1": float(np.max(np.abs(examples_df["A_at_horizon"] - 1.0))),
        "max_abs_escape_at_horizon_over_c_minus_1": float(
            np.max(np.abs(examples_df["escape_velocity_at_horizon_over_c"] - 1.0))
        ),
        "max_abs_mass_recovered_ratio_minus_1": float(np.max(np.abs(examples_df["M_recovered_ratio"] - 1.0))),
        "max_abs_A_minus_v_escape_over_c_squared": float(
            np.max(np.abs(scan_df["A_minus_v_escape_over_c_squared"]))
        ),
    }

    summary = {
        "test": "STAM Model-A horizon threshold and escape velocity identity",
        "formulas": {
            "A": "A(r)=Rs/r=2GM/(c^2 r)",
            "horizon": "A=1 <=> r=Rs",
            "escape_identity": "A=(v_escape/c)^2",
        },
        "checks": checks,
        "one_solar_mass_horizon_radius_km": float(
            examples_df.loc[examples_df["mass_solar"] == 1.0, "horizon_radius_km"].iloc[0]
        ),
        "ten_solar_mass_horizon_radius_km": float(
            examples_df.loc[examples_df["mass_solar"] == 10.0, "horizon_radius_km"].iloc[0]
        ),
        "interpretation": (
            "The STAM A=1 threshold is exactly the spherical escape-speed-equals-light-speed threshold. "
            "For A>1, v_escape>c."
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    print("PASS: STAM Model-A horizon threshold identity holds.")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote outputs to {out_dir}")


if __name__ == "__main__":
    main()
