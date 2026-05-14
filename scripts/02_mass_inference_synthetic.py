#!/usr/bin/env python3
"""Synthetic mass-inference consistency check for STAM_model-A."""

from __future__ import annotations

from pathlib import Path
import json
import math

from stam_model_a.constants import C_SI, G_SI
from stam_model_a.local import schwarzschild_radius
from stam_model_a.mass_estimators import (
    mass_from_acceleration,
    mass_from_gravitational_redshift,
    mass_from_horizon_radius,
    mass_from_lensing_deflection,
    mass_from_orbital_velocity,
    mass_from_shapiro_coefficient,
)


def main() -> None:
    out_dir = Path("results/identity_checks")
    out_dir.mkdir(parents=True, exist_ok=True)

    mass = 7.5e30
    r = 2.0e11
    rs = schwarzschild_radius(mass)

    observables = {
        "horizon_radius_m": rs,
        "acceleration_m_s2": G_SI * mass / r**2,
        "orbital_velocity_m_s": math.sqrt(G_SI * mass / r),
        "shapiro_coefficient_s": 2.0 * G_SI * mass / C_SI**3,
        "gravitational_redshift_weak": G_SI * mass / (C_SI**2 * r),
        "lensing_deflection_rad": 4.0 * G_SI * mass / (C_SI**2 * r),
    }

    estimates = {
        "horizon": mass_from_horizon_radius(observables["horizon_radius_m"]),
        "acceleration": mass_from_acceleration(observables["acceleration_m_s2"], r),
        "orbital_velocity": mass_from_orbital_velocity(observables["orbital_velocity_m_s"], r),
        "shapiro_coefficient": mass_from_shapiro_coefficient(observables["shapiro_coefficient_s"]),
        "gravitational_redshift": mass_from_gravitational_redshift(observables["gravitational_redshift_weak"], r),
        "lensing_deflection": mass_from_lensing_deflection(observables["lensing_deflection_rad"], r),
    }

    ratios = {name: value / mass for name, value in estimates.items()}
    result = {
        "input_mass_kg": mass,
        "test_radius_m": r,
        "observables": observables,
        "mass_estimates_kg": estimates,
        "ratios_to_input_mass": ratios,
        "max_abs_ratio_error": max(abs(x - 1.0) for x in ratios.values()),
    }

    out_path = out_dir / "synthetic_mass_inference.json"
    out_path.write_text(json.dumps(result, indent=2, sort_keys=True))
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
