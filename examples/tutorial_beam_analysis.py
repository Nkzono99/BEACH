#!/usr/bin/env python3
"""Compare the electron-beam tutorial run with the expected charging physics.

Usage::

    python tutorial_beam_analysis.py outputs/tutorial --config beach.toml --save-plot outputs/tutorial/charging.png

The script reads the per-batch charge and potential histories and prints, for every
reported batch, the absorbed fraction of the beam, the plate potential, the potential
at the centre of the beam source plane, and their difference. An electron of normal
energy E_b reaches the plate only while ``phi_plate - phi_source > -E_b/e``.

The source-plane potential is a sum over element charges treated as point charges.
This is accurate when the source plane is far from the plate compared with the element
size (50 mm against 2 mm in the tutorial); use the field kernel of the Python API for
points close to the surface.
"""

from __future__ import annotations

import argparse
import math
import tomllib
from pathlib import Path

import numpy as np

ELEMENTARY_CHARGE = 1.602176634e-19
EPSILON_0 = 8.8541878188e-12
COULOMB_K = 1.0 / (4.0 * math.pi * EPSILON_0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output", type=Path, help="output.dir of the run")
    parser.add_argument("--config", type=Path, default=Path("beach.toml"), help="the beach.toml that was run")
    parser.add_argument("--every", type=int, default=25, help="print every N batches")
    parser.add_argument("--save-plot", type=Path, help="save a two-panel charging plot (needs matplotlib)")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = tomllib.loads(args.config.read_text(encoding="utf-8"))
    species = config["particles"]["species"][0]
    source = species["source"]
    sampling = species["sampling"]
    mass = species["mass_kg"]
    charge = species["charge_c"]
    drift = np.asarray(species["distribution"]["drift_velocity_m_s"], dtype=float)
    batch_s = config["run"]["batch"]["duration_s"]
    macro = sampling["volume_macro_particles_per_batch"]
    weight = sampling["weight"]

    beam_energy_ev = 0.5 * mass * float(drift @ drift) / ELEMENTARY_CHARGE
    beam_current_a = abs(charge) * weight * macro / batch_s
    source_low = np.asarray(source["pos_low"], dtype=float)
    source_high = np.asarray(source["pos_high"], dtype=float)
    source_centre = 0.5 * (source_low + source_high)
    beam_area = float(np.prod((source_high - source_low)[:2]))

    triangles = np.loadtxt(args.output / "mesh_triangles.csv", delimiter=",", skiprows=1)
    vertices = triangles[:, 1:10].reshape(-1, 3, 3)
    centroids = vertices.mean(axis=1)
    areas = 0.5 * np.linalg.norm(
        np.cross(vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]), axis=1
    )
    nelem = len(centroids)

    charges = np.loadtxt(args.output / "charge_history.csv", delimiter=",", skiprows=1)
    potentials = np.loadtxt(args.output / "potential_history.csv", delimiter=",", skiprows=1)
    nbatch = int(charges[:, 0].max())
    q = charges[:, 4].reshape(nbatch, nelem)
    phi = potentials[:, 2].reshape(nbatch, nelem)

    total = q.sum(axis=1)
    absorbed = -np.diff(np.concatenate([[0.0], total])) / (abs(charge) * weight * macro)
    distance = np.linalg.norm(centroids - source_centre, axis=1)
    phi_source = COULOMB_K * (q / distance).sum(axis=1)
    relative = phi - phi_source[:, None]
    time_us = np.arange(1, nbatch + 1) * batch_s * 1.0e6

    plate_area = float(areas.sum())
    disk_capacitance = 8.0 * EPSILON_0 * math.sqrt(plate_area / math.pi)
    print(f"beam_energy_eV={beam_energy_ev:.3f}")
    print(f"beam_current_A={beam_current_a:.4e} beam_current_density_uA_cm2={beam_current_a / beam_area * 1.0e2:.3f}")
    print(f"geometric_hit_fraction={min(1.0, plate_area / beam_area):.3f}")
    print("batch time_us absorbed_fraction plate_charge_C phi_plate_V phi_source_V phi_plate_minus_source_V(min/mean/max)")
    for n in sorted(set(range(0, nbatch, args.every)) | {nbatch - 1}):
        print(
            f"{n + 1:5d} {time_us[n]:7.3f} {absorbed[n]:6.3f} {total[n]: .4e} {phi[n].mean():7.3f} "
            f"{phi_source[n]:7.3f} {relative[n].min():7.3f}/{relative[n].mean():7.3f}/{relative[n].max():7.3f}"
        )
    print(f"capacitance_from_run_pF={abs(total[-1] / phi[-1].mean()) * 1.0e12:.3f}")
    print(f"capacitance_equal_area_disk_pF={disk_capacitance * 1.0e12:.3f}")

    if args.save_plot:
        import matplotlib.pyplot as plt

        fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(6.0, 5.0))
        top.plot(time_us, relative.mean(axis=1), label=r"$\phi_{\rm plate}-\phi_{\rm source}$ (mean)")
        top.fill_between(time_us, relative.min(axis=1), relative.max(axis=1), alpha=0.25, label="min to max")
        top.axhline(-beam_energy_ev, color="k", linestyle="--", label=r"$-E_b/e$")
        top.set_ylabel("potential [V]")
        top.legend(loc="upper right")
        bottom.plot(time_us, absorbed)
        bottom.set_ylabel("absorbed fraction")
        bottom.set_xlabel(r"time [$\mu$s]")
        bottom.set_ylim(bottom=0.0)
        fig.tight_layout()
        fig.savefig(args.save_plot, dpi=150)
        print(f"saved_plot={args.save_plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
