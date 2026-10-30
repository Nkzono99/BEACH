#!/usr/bin/env python3
"""Compare the floating-plate tutorial run with the current-balance estimate.

Usage::

    python tutorial_floating_analysis.py outputs/tutorial --config beach.toml --save-plot outputs/tutorial/floating.png

Electrons and protons enter through the top face of the box as one-sided Maxwellian
fluxes ``J0 = e n sqrt(T / (2 pi m))``. The plate collects every proton and the electrons
whose normal energy exceeds ``-e dphi``, i.e. ``J_e0 exp(e dphi / Te)``, where ``dphi`` is
the plate potential relative to the top face. It charges until the two currents balance:

    dphi = (Te / e) ln(J_i0 / J_e0)

The script prints ``dphi`` and the plate charge density over time next to the analytic
relaxation ``exp(-e dphi / Te) = R - (R - 1) exp(-t / tau)`` (``R = J_e0 / J_i0``), their
steady-state averages, and the absorbed counts of each species between the periodic
checkpoint and the end of the run.
"""

from __future__ import annotations

import argparse
import math
import tomllib
from pathlib import Path

import numpy as np

from beach import load_run

ELEMENTARY_CHARGE = 1.602176634e-19
EPSILON_0 = 8.8541878188e-12


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output", type=Path, help="output.dir of the run")
    parser.add_argument("--config", type=Path, default=Path("beach.toml"), help="the beach.toml that was run")
    parser.add_argument("--every", type=int, default=1000, help="print every N batches")
    parser.add_argument("--save-plot", type=Path, help="save a two-panel plot (needs matplotlib)")
    return parser.parse_args()


def inflow_current_density(species: dict) -> tuple[float, float]:
    """Return (one-sided Maxwellian current density [A/m^2], temperature [eV]) of one species."""

    distribution = species["distribution"]
    density = distribution.get("number_density_m3", distribution.get("number_density_cm3", 0.0) * 1.0e6)
    temperature_ev = distribution["temperature_ev"]
    speed = math.sqrt(temperature_ev * ELEMENTARY_CHARGE / (2.0 * math.pi * species["mass_kg"]))
    return abs(species["charge_c"]) * density * speed, temperature_ev


def window_ledgers(output: Path, final_batch: int):
    """Return the ledger of the earliest periodic checkpoint before the final batch, if any."""

    candidates = []
    for slot in sorted((output / "checkpoints").glob("slot*")):
        if not (slot / "summary.txt").is_file():
            continue
        ledger = load_run(slot).ledger
        if ledger and ledger[0].batch < final_batch:
            candidates.append(ledger)
    return min(candidates, key=lambda rows: rows[0].batch) if candidates else None


def main() -> int:
    args = parse_args()
    config = tomllib.loads(args.config.read_text(encoding="utf-8"))
    species = config["particles"]["species"]
    electron = next(s for s in species if s["charge_c"] < 0.0)
    proton = next(s for s in species if s["charge_c"] > 0.0)
    j_e0, te = inflow_current_density(electron)
    j_i0, _ = inflow_current_density(proton)
    batch_s = config["run"]["batch"]["duration_s"]
    box_min = np.asarray(config["domain"]["box_min"], dtype=float)
    box_max = np.asarray(config["domain"]["box_max"], dtype=float)
    area = float(np.prod((box_max - box_min)[:2]))
    if config["fields"]["periodic"]["lower_boundary_model"] != "e_bottom_zero":
        # With no field below the plate, all of sigma / epsilon_0 points into the gap.
        raise ValueError('the estimate assumes fields.periodic.lower_boundary_model = "e_bottom_zero"')

    run = load_run(args.output)
    vertices = run.mesh.triangles
    charge_history = run.history.charge
    potential_history = run.history.potential
    top = run.top_reference_history
    if vertices is None or charge_history is None or potential_history is None or not top:
        raise ValueError("the tutorial needs saved mesh, charge, potential and top reference histories")
    batches = charge_history.batch_indices
    if not np.array_equal(batches, potential_history.batch_indices) or not np.array_equal(
        batches, [row.batch for row in top]
    ):
        raise ValueError("charge, potential and top reference histories must share their batch index")

    gap = float(box_max[2] - vertices[:, :, 2].mean())
    capacitance_per_area = EPSILON_0 / gap
    ratio = j_e0 / j_i0
    dphi_eq = te * math.log(j_i0 / j_e0)
    tau = capacitance_per_area * te / j_i0

    sigma = charge_history.as_array().sum(axis=0) / area
    dphi = potential_history.as_array().mean(axis=0) - np.array([row.potential_mean_v for row in top])
    time_s = batches * batch_s
    dphi_model = -te * np.log(ratio - (ratio - 1.0) * np.exp(-time_s / tau))

    print(f"electron_temperature_eV={te:.3f} gap_m={gap:.4f} area_m2={area:.4e}")
    print(f"inflow_current_density_A_m2 electron={j_e0:.4e} proton={j_i0:.4e} ratio={ratio:.2f}")
    print(f"expected_dphi_V={dphi_eq:.3f} expected_charge_density_C_m2={capacitance_per_area * dphi_eq:.4e}")
    print(f"relaxation_time_s={tau:.4f} relaxation_batches={tau / batch_s:.0f}")
    print("batch   time_s   dphi_V  model_V  charge_density_C_m2")
    rows = {n for n, batch in enumerate(batches) if (batch - batches[0]) % args.every == 0} | {len(batches) - 1}
    for n in sorted(rows):
        print(f"{batches[n]:6d} {time_s[n]:7.3f} {dphi[n]:8.3f} {dphi_model[n]:8.3f}  {sigma[n]: .4e}")

    final_batch = run.summary.batches
    start = window_ledgers(args.output, final_batch)
    start_batch = start[0].batch if start else int(batches[len(batches) // 2])
    steady = batches > start_batch
    print(f"steady_window_batches={start_batch + 1}-{final_batch}")
    print(
        f"steady_dphi_V mean={dphi[steady].mean():.4f} std={dphi[steady].std():.4f} expected={dphi_eq:.4f}"
    )
    print(f"steady_charge_density_C_m2 mean={sigma[steady].mean():.4e}")

    if start and run.ledger:
        window_s = (final_batch - start_batch) * batch_s
        expected_electron = float(np.exp(dphi[steady] / te).mean())
        for name, spec, j0, expected in (
            ("electron", electron, j_e0, expected_electron),
            ("proton", proton, j_i0, 1.0),
        ):
            index = species.index(spec) + 1
            end_row = next(row for row in run.ledger if row.species_idx == index)
            start_row = next(row for row in start if row.species_idx == index)
            injected = end_row.injected_count - start_row.injected_count
            absorbed = end_row.absorbed_count - start_row.absorbed_count
            current = abs(end_row.absorbed_on_surface_c - start_row.absorbed_on_surface_c) / (area * window_s)
            print(
                f"window_{name} injected={injected} absorbed={absorbed} fraction={absorbed / injected:.5f} "
                f"expected_fraction={expected:.5f} absorbed_current_A_m2={current:.4e} "
                f"expected_current_A_m2={j0 * expected:.4e}"
            )
    print(f"survived_max_step={run.summary.survived_max_step}")

    if args.save_plot:
        import matplotlib.pyplot as plt

        net = np.diff(sigma) / np.diff(time_s) / j_i0
        mid = 0.5 * (time_s[1:] + time_s[:-1])
        model_mid = 1.0 - ratio / (ratio - (ratio - 1.0) * np.exp(-mid / tau))
        fig, (upper, lower) = plt.subplots(2, 1, sharex=True, figsize=(6.0, 5.0))
        upper.plot(time_s, dphi, label="BEACH")
        upper.plot(time_s, dphi_model, "k--", lw=1.0, label="current-balance model")
        upper.axhline(dphi_eq, color="tab:red", lw=0.8, label=rf"$(T_e/e)\ln(J_{{i0}}/J_{{e0}})$ = {dphi_eq:.2f} V")
        upper.set_ylabel(r"$\phi_{\rm plate}-\phi_{\rm top}$ [V]")
        upper.legend(loc="upper right")
        lower.plot(mid, net, lw=0.8, label="BEACH")
        lower.plot(mid, model_mid, "k--", lw=1.0, label="model")
        lower.axhline(0.0, color="tab:red", lw=0.8)
        lower.set_ylabel(r"net current / $J_{i0}$")
        lower.set_xlabel("time [s]")
        lower.set_xscale("log")
        lower.legend(loc="lower right")
        fig.tight_layout()
        fig.savefig(args.save_plot, dpi=150)
        print(f"saved_plot={args.save_plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
