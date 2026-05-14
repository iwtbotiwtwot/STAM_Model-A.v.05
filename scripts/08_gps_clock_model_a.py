#!/usr/bin/env python3
"""STAM Model-A GPS satellite clock adjustment test.

Non-distance Model-A test.

STAM local accumulation:
    A(r) = 2GM/(c^2 r)

Weak-field clock mapping tested:
    dτ/dt ≈ 1 - A/2

Then the gravitational clock-rate difference between Earth geoid/surface and
GPS orbit is approximately:
    Δf/f_grav = (A_geoid - A_orbit)/2

A circular-orbit kinematic term is then added:
    Δf/f_kin = -v^2/(2c^2)
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd

C = 299_792_458.0
MU_EARTH = 3.986004418e14
R_EARTH = 6_378_137.0
GPS_RADIUS = 26_560_000.0
DAY = 86_400.0
F_NOMINAL = 10.23e6

# Effective geoid potential magnitude divided by c^2 used in GPS relativity treatments.
# This captures the real geoid reference rather than a perfectly spherical Earth surface.
U0_OVER_C2 = 6.969290134e-10


def accumulation_spherical(r_m: float) -> float:
    return 2.0 * MU_EARTH / (C**2 * r_m)


def gps_clock_geoid_reference(radius_m: float = GPS_RADIUS) -> dict[str, float]:
    A_geoid_eff = 2.0 * U0_OVER_C2
    A_orbit = accumulation_spherical(radius_m)
    grav_fraction = (A_geoid_eff - A_orbit) / 2.0

    v = math.sqrt(MU_EARTH / radius_m)
    kin_fraction = -v**2 / (2.0 * C**2)
    net_fraction = grav_fraction + kin_fraction

    return {
        "orbit_radius_m": radius_m,
        "orbit_altitude_above_wgs84_equator_m": radius_m - R_EARTH,
        "orbit_speed_m_s": v,
        "A_geoid_effective": A_geoid_eff,
        "A_orbit": A_orbit,
        "delta_A_geoid_minus_orbit": A_geoid_eff - A_orbit,
        "stam_grav_fractional_shift": grav_fraction,
        "stam_grav_us_per_day": grav_fraction * DAY * 1e6,
        "kinematic_fractional_shift": kin_fraction,
        "kinematic_us_per_day": kin_fraction * DAY * 1e6,
        "net_fractional_shift_satellite_vs_geoid": net_fraction,
        "net_us_per_day_satellite_gains": net_fraction * DAY * 1e6,
        "factory_clock_fractional_offset_needed": -net_fraction,
        "factory_adjusted_frequency_mhz": F_NOMINAL * (1.0 - net_fraction) / 1e6,
    }


def main() -> None:
    out_dir = Path("results/gps_clock")
    out_dir.mkdir(parents=True, exist_ok=True)

    result = gps_clock_geoid_reference()
    (out_dir / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True))

    rows = []
    for altitude_km in [400, 1000, 5000, 10000, 15000, (GPS_RADIUS - R_EARTH) / 1000, 20200, 25000, 30000, 35786]:
        radius = R_EARTH + altitude_km * 1000.0
        r = gps_clock_geoid_reference(radius)
        rows.append(
            {
                "altitude_km": altitude_km,
                "A_orbit": r["A_orbit"],
                "stam_grav_us_day": r["stam_grav_us_per_day"],
                "kinematic_us_day": r["kinematic_us_per_day"],
                "net_us_day": r["net_us_per_day_satellite_gains"],
                "factory_offset_fraction_needed": r["factory_clock_fractional_offset_needed"],
            }
        )
    pd.DataFrame(rows).to_csv(out_dir / "altitude_scan_clock_effects.csv", index=False)

    print("PASS: STAM Model-A GPS weak-field clock adjustment computed.")
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"Wrote {out_dir / 'summary.json'}")
    print(f"Wrote {out_dir / 'altitude_scan_clock_effects.csv'}")


if __name__ == "__main__":
    main()
