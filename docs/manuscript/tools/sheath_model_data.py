"""Compute the numerical examples of the outer-sheath-model chapter.

Run on a compute node with the sheath-model Python package importable, e.g.

    PYTHONPATH=/path/to/sheath-model python tools/sheath_model_data.py figures/data/sheath_model

The CSV files written here are read by tools/sheath_model_figures.py, which emits the
TikZ sources included by the manuscript. Only public sheath-model APIs are used.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.integrate import quad

import sheath_model
from sheath_model import (ProfileOptions, SearchFailure, SearchOptions, SheathSolver, ZhaoParams)

EPS0 = 8.8541878128e-12
QE = 1.602176634e-19
ME = 9.1093837015e-31

# Common illumination and solar-wind conditions of the examples (sheath-model README).
BASE = ZhaoParams(sun_elevation_deg=60.0, ion_density_m3=8.7e6, photoelectron_reference_density_m3=64e6,
                  electron_temperature_ev=12.0, ion_temperature_ev=0.0, photoelectron_temperature_ev=2.2,
                  solar_wind_speed_mps=468e3, electron_drift_mode="zero", ion_drift_mode="normal")
EXAMPLES = (("A", 60.0), ("B", 20.0), ("C", 10.0))
FIELD_FLOOR_V_M = 2e-3
SOLVER = SheathSolver(profile=ProfileOptions(zmax_hat=200.0, n_profile_grid=1200, n_type_a_grid=8000))
DRIFT_SOLVER = SheathSolver(search=SearchOptions(upstream_field_floor_v_m=FIELD_FLOOR_V_M),
                            profile=ProfileOptions(zmax_hat=200.0))


def write_csv(path: Path, header, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        for row in rows:
            writer.writerow([f"{v:.9g}" if isinstance(v, float) else v for v in row])


def fluxes(root):
    """Number fluxes [m^-2 s^-1] at the wall from the closed-form moments of the chapter."""
    p = root.inputs
    vte = math.sqrt(2 * QE * p.electron_temperature_ev / ME)
    vtp = math.sqrt(2 * QE * p.photoelectron_temperature_ev / ME)
    u = p.electron_drift_mps / vte
    b = math.sqrt(max(0.0, -root.minimum_potential_v / p.electron_temperature_ev))
    gamma_e = root.ambient_electron_density_m3 * vte / (2 * math.sqrt(math.pi)) * (
        math.exp(-(b - u) ** 2) + math.sqrt(math.pi) * u * math.erfc(b - u))
    gamma_i = p.ion_density_m3 * p.ion_entry_speed_mps
    gamma_emit = p.emission_density_m3 * vtp / (2 * math.sqrt(math.pi))
    barrier = root.surface_potential_v - root.minimum_potential_v
    gamma_esc = gamma_emit * math.exp(-barrier / p.photoelectron_temperature_ev)
    return gamma_e, gamma_i, gamma_emit, gamma_esc


def field_squared_from_upstream(root, phi, side="upper"):
    """E^2(phi) [V^2/m^2] = (2/eps0) int_phi^0 rho dpsi on a segment connected to upstream."""
    def rho(psi):
        return float(root.density(psi, side=side).charge_c_m3) if root.branch == "A" else \
            float(root.density(psi).charge_c_m3)
    value, _ = quad(rho, phi, 0.0, limit=200, epsabs=0.0, epsrel=1e-11)
    return 2.0 / EPS0 * value


def field_squared_from_minimum(root, phi, side):
    def rho(psi):
        return float(root.density(psi, side=side).charge_c_m3)
    value, _ = quad(rho, root.minimum_potential_v, phi, limit=200, epsabs=0.0, epsrel=1e-11)
    return -2.0 / EPS0 * value


def profiles(out: Path):
    meta = []
    type_a = None
    for branch, alpha in EXAMPLES:
        profile = SOLVER.solve_profile(replace(BASE, sun_elevation_deg=alpha), branch=branch)
        root = profile.equilibrium
        d = profile.density
        columns = [profile.z_m, profile.potential_v, profile.electric_field_v_m, d.ion_m3, d.electron_free_m3,
                   d.electron_reflected_m3, d.emitted_free_m3, d.emitted_captured_m3,
                   d.charge_c_m3 / QE]
        # Keep the reconstructed nodes but thin them to about 0.1 m (plot resolution); keep both ends and
        # the potential minimum so that no feature of the solution is interpolated.
        z = np.asarray(profile.z_m)
        keep = np.zeros(len(z), dtype=bool)
        keep[[0, -1, int(np.argmin(profile.potential_v))]] = True
        last = -math.inf
        for i, value in enumerate(z):
            if value - last >= 0.1 or i < 40:
                keep[i], last = True, value
        rows = zip(*(np.asarray(c)[keep] for c in columns))
        write_csv(out / f"profile_{branch}.csv",
                  ["z_m", "phi_v", "e_v_m", "n_i_m3", "n_e_pass_m3", "n_e_refl_m3", "n_pe_free_m3",
                   "n_pe_capt_m3", "rho_over_e_m3"], ([float(x) for x in r] for r in rows))
        ge, gi, gemit, gesc = fluxes(root)
        meta.append(dict(branch=branch, alpha_deg=alpha, phi0_v=root.surface_potential_v,
                         phim_v=root.minimum_potential_v, n_e_norm_m3=root.ambient_electron_density_m3,
                         z_min_m=profile.turning_height_m, wall_field_v_m=float(profile.electric_field_v_m[0]),
                         length_scale_m=profile.equilibrium.inputs.length_scale_m,
                         gamma_e=ge, gamma_i=gi, gamma_emit=gemit, gamma_esc=gesc,
                         net_current_a_m2=QE * (gi + gesc - ge), residual=root.residual_norm,
                         z_end_m=float(profile.z_m[-1])))
        if branch == "A":
            type_a = profile
    (out / "examples.json").write_text(json.dumps(meta, indent=1))
    return type_a


def first_integral(out: Path, profile):
    """rho(phi) and E^2(phi) on both Type A segments, by quadrature independent of the profile grid."""
    root = profile.equilibrium
    phi0, phim = root.surface_potential_v, root.minimum_potential_v
    rows = []
    for side, hi in (("lower", phi0), ("upper", 0.0)):
        t = np.linspace(0.0, 1.0, 241)
        for phi in phim + (hi - phim) * t ** 2:
            rho = float(root.density(phi, side=side).charge_c_m3) / QE
            e2 = 0.0 if phi == phim else field_squared_from_minimum(root, phi, side)
            rows.append((side, float(phi), rho, e2))
    write_csv(out / "first_integral_A.csv", ["side", "phi_v", "rho_over_e_m3", "e2_v2_m2"], rows)


def vdfs(out: Path, profile):
    zm = profile.turning_height_m
    rows = []
    for label, z in (("wall", 0.0), ("lower", 0.5 * zm), ("upper", 40.0)):
        for species, keys in (("swe", ("g_free_incoming", "g_reflected_incoming", "g_reflected_outgoing")),
                              ("phe", ("g_free_outgoing", "g_captured_outgoing", "g_captured_returning"))):
            data = profile.vdf(z, species=species, n_v=1601)[species]
            for i, v in enumerate(data["vz_mps"]):
                rows.append((label, float(z), species, float(v), *[float(data[k][i]) for k in keys]))
    write_csv(out / "vdf_A.csv", ["position", "z_m", "species", "vz_mps", "g1_s_m4", "g2_s_m4", "g3_s_m4"], rows)


def _sweep_branch(branch):
    rows, previous = [], None
    for alpha in np.arange(1.0, 89.01, 1.0):
        inputs = replace(BASE, sun_elevation_deg=float(alpha))
        root = None
        for guess in ((previous,) if previous is not None else ()) + (None,):
            try:
                root = SOLVER.solve_equilibrium(inputs, branch=branch, initial_guess=guess)
                break
            except SearchFailure:
                continue
        if root is None:
            previous = None
            continue
        previous = root
        ge, gi, gemit, gesc = fluxes(root)
        rows.append((branch, float(alpha), root.surface_potential_v, root.minimum_potential_v,
                     root.ambient_electron_density_m3, ge, gi, gemit, gesc))
    return rows


def sweep(out: Path, pool):
    rows = [r for part in pool.map(_sweep_branch, "ABC") for r in part]
    write_csv(out / "alpha_sweep.csv", ["branch", "alpha_deg", "phi0_v", "phim_v", "n_e_norm_m3", "gamma_e",
                                        "gamma_i", "gamma_emit", "gamma_esc"], rows)


def _map_point(args):
    alpha, ratio = args
    inputs = replace(BASE, sun_elevation_deg=alpha, photoelectron_reference_density_m3=ratio * BASE.ion_density_m3)
    found, rejected = 0, 0
    for bit, branch in ((1, "A"), (2, "B"), (4, "C")):
        try:
            SOLVER.solve_equilibrium(inputs, branch=branch)
            found |= bit
        except SearchFailure as exc:
            if exc.diagnostics.rejected["ABC".index(branch)] > 0:
                rejected |= bit
    return alpha, ratio, found, rejected


def type_map(out: Path, pool, n_alpha=89, n_ratio=49):
    alphas = np.linspace(1.0, 89.0, n_alpha)
    ratios = np.exp(np.linspace(math.log(0.5), math.log(16.0), n_ratio))
    tasks = [(float(a), float(r)) for r in ratios for a in alphas]
    rows = list(pool.map(_map_point, tasks, chunksize=8))
    write_csv(out / "type_map.csv", ["alpha_deg", "source_ratio", "types_found", "types_rejected"], rows)


def upstream_band(out: Path):
    """E^2 next to upstream for the Type A root with and without inward electron drift."""
    rows, meta = [], []
    for mode, solver in (("zero", SOLVER), ("normal", DRIFT_SOLVER)):
        root = solver.solve_equilibrium(replace(BASE, electron_drift_mode=mode), branch="A")
        meta.append(dict(drift=mode, phi0_v=root.surface_potential_v, phim_v=root.minimum_potential_v,
                         band_v=root.upstream_negative_band_v, deficit_v_m=root.upstream_field_deficit_v_m,
                         n_e_norm_m3=root.ambient_electron_density_m3))
        depth_max = -root.minimum_potential_v
        for depth in np.exp(np.linspace(math.log(1e-5), math.log(0.999 * depth_max), 161)):
            rows.append((mode, float(depth), field_squared_from_upstream(root, -float(depth))))
    write_csv(out / "upstream_band.csv", ["drift", "depth_v", "e2_v2_m2"], rows)
    (out / "upstream_band.json").write_text(json.dumps(meta, indent=1))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    parser.add_argument("--only", nargs="*", default=None)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    jobs = set(args.only or ("profiles", "sweep", "map", "band"))
    package = Path(sheath_model.__file__).resolve().parents[1]
    try:
        commit = subprocess.run(["git", "-C", str(package), "rev-parse", "HEAD"], capture_output=True,
                                text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    (args.output / "provenance.json").write_text(json.dumps(dict(
        sheath_model_commit=commit, base_inputs={k: getattr(BASE, k) for k in BASE.__dataclass_fields__},
        field_floor_v_m=FIELD_FLOOR_V_M), indent=1))
    profile = None
    if "profiles" in jobs:
        profile = profiles(args.output)
        first_integral(args.output, profile)
        vdfs(args.output, profile)
    if "band" in jobs:
        upstream_band(args.output)
    if jobs & {"sweep", "map"}:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            if "sweep" in jobs:
                sweep(args.output, pool)
            if "map" in jobs:
                type_map(args.output, pool)


if __name__ == "__main__":
    main()
