"""Bounded time-dependent probes with the two #38 roots' own boundary sources.

This is a convergence probe, not certification. The existing oracle uses warm
Maxwellian ions, so the cold-ion limit must also be assessed explicitly.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import replace
from pathlib import Path

from beach.outer_kinetic import (
    EPS0,
    ELEMENTARY_CHARGE,
    ELECTRON_MASS,
    KineticQuery,
    VelocityGridConfig,
    load_outer_kinetic_config,
    write_kinetic_atlas,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    base = load_outer_kinetic_config(Path(__file__).resolve().parents[1] / "outer_kinetic_reference.toml")
    density, te = 8.7e6, 12.0
    velocity = math.sqrt(ELEMENTARY_CHARGE * te / ELECTRON_MASS)
    debye = math.sqrt(EPS0 * te / (density * ELEMENTARY_CHARGE))
    ion_drift = 10 * math.sqrt(ELEMENTARY_CHARGE * te / base.ion.mass_kg)
    query = KineticQuery(0.1 * EPS0 * te / debye, 0.3 * density * velocity, 0.2 * te)
    rows = []
    for branch, amplitude, phi, q in (
        ("B", 0.5003210848357, 0.0228867312916, 0.2675613319),
        ("N", 0.8974844341329, -0.1329182234611, 0.2781919119),
    ):
        # Change one axis per variant; no failed case is converted to a table.
        for name, nz, nv, length, ti, duration, vmax in (
            ("base", 32, 128, 3, 0.12, 6, 6),
            ("space", 64, 128, 3, 0.12, 6, 6),
            ("velocity", 32, 256, 3, 0.12, 6, 6),
            ("length", 64, 128, 6, 0.12, 6, 6),
            ("cold", 32, 256, 3, 0.03, 6, 6),
            ("time", 32, 128, 3, 0.12, 12, 6),
            ("range", 32, 192, 3, 0.12, 6, 9),
        ):
            transit = length * debye / ion_drift
            window = duration * transit / 3
            grid = VelocityGridConfig(nv, -vmax * velocity, vmax * velocity)
            config = replace(
                base, nz=nz, z_length_m=length * debye, max_time_s=duration * transit,
                electron=replace(base.electron, number_density_m3=amplitude * density, temperature_ev=te,
                                 drift_velocity_mps=0.0, grid=grid),
                ion=replace(base.ion, number_density_m3=density, temperature_ev=ti,
                            drift_velocity_mps=-ion_drift,
                            grid=VelocityGridConfig(4 * nv, -0.5 * velocity, 0.05 * velocity)),
                photoelectron_grid=grid,
                certification=replace(base.certification, warmup_time_s=window, averaging_window_s=window,
                                      sample_interval_s=window / 40, far_field_abs_v_m=0.01 * te / debye),
            )
            result = write_kinetic_atlas(config, [query], args.output / f"{branch}-{name}")[0]
            row = dict(branch=branch, variant=name, source_amplitude=amplitude,
                       static_phi_v=phi * te, static_escape_flux=q * density * velocity, **result.raw_row())
            rows.append(row)
            (args.output / "comparison.json").write_text(json.dumps(rows, indent=2))
            print(branch, name, result.classification, result.response[0], result.failure_reason, flush=True)


if __name__ == "__main__":
    main()
