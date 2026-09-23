#!/usr/bin/env python3
"""Independently diagnose a failed Type-B matching-plane endpoint (NumPy/SciPy required).

Run on a compute node:
  python tools/diagnose_matching_plane_failure.py run.err --config beach.toml --output diagnosis
  python tools/diagnose_matching_plane_failure.py spectrum.csv --config beach.toml --query captured-query.txt

CSV columns are energy_low_ev, energy_high_ev, flux_m2_s.
Bin fluxes are integrated number fluxes; no production Python/Fortran is imported.
Ne follows neutral infinity exactly. Electron density integrals use upstream
velocity, PE integrals use exact constant-flux-bin energy primitives. The BE
condition is D - D_before = dt*qe*(Gamma_i + Gamma_PE_escape - Gamma_e).
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
import re
import time
from functools import lru_cache
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

try:
    import numpy as np
    from scipy.integrate import quad
    from scipy.optimize import brentq, minimize_scalar
except ImportError as error:
    raise SystemExit('This optional analysis tool requires NumPy and SciPy.') from error

QE = 1.602176634e-19
EPS0 = 8.8541878128e-12
# Species and normalization constants are populated from the required config.
ME = MI = NI = TE = VE = VI = 0.0
VPE = VTE = U = F0 = ELECTRON_FLUX_PER_DENSITY = ION_FLUX = ION_LIMIT = 0.0


def read_failure_capture(path):
    text = path.read_text()
    marker = 'matching-plane failed spectrum begin'
    if marker not in text:
        return text, {}
    before, block = text.rsplit(marker, 1)
    end_marker = 'matching-plane failed spectrum end'
    if end_marker not in block:
        raise ValueError('the last spectrum receipt is incomplete')
    spectrum = block.split(end_marker, 1)[0].strip()
    query = {}
    def last_numbers(prefix):
        matches = re.findall(re.escape(prefix)+r'([^\n]*)', before)
        if not matches:
            raise ValueError(f'missing failure receipt: {prefix}')
        return [float(value.replace('D', 'E')) for value in matches[-1].split()]
    endpoint = last_numbers('matching-plane failed endpoint: D_before, D_seed [C/m2], duration [s]=')
    feedback = last_numbers('matching-plane failed feedback: PE flux, PE energy, electron flux, ion flux=')
    if len(endpoint) != 3 or len(feedback) != 4:
        raise ValueError('malformed endpoint/feedback receipt')
    query.update(displacement_before_c_m2=str(endpoint[0]), displacement_seed_c_m2=str(endpoint[1]),
                 trial_batch_duration_s=str(endpoint[2]), guess=' '.join(map(str, feedback)))
    seed_prefix = 'matching-plane failed seed: phi_H, phi_m [V], electron density [m-3]='
    if seed_prefix in before:
        seed = last_numbers(seed_prefix)
        if len(seed) == 3:
            query.update(root_trial_phi0_v=str(seed[0]), root_trial_phi_m_v=str(seed[1]),
                         root_trial_electron_density_m3=str(seed[2]))
    return spectrum, query


def configure_from_toml(path):
    global ME, MI, NI, TE, VE, VI, VPE, VTE, U, F0, ELECTRON_FLUX_PER_DENSITY, ION_FLUX, ION_LIMIT
    config = tomllib.loads(path.read_text())
    model = config['surface_current_model']
    if model.get('response_backend') != 'zhao_online' or model.get('photoelectron_closure') != 'energy_spectrum':
        raise ValueError('diagnosis requires zhao_online with the recorded energy_spectrum closure')
    species = {item['species_key']: item for item in config['particles']['species']}
    electron = species[model['electron_species']]
    ion = species[model['ion_species']]
    photo = species[model['photoelectron_species']]
    for item, charge in [(electron, -QE), (ion, QE), (photo, -QE)]:
        if not math.isclose(float(item['q_particle']), charge, rel_tol=1e-12, abs_tol=0.0):
            raise ValueError('this Type-B oracle currently supports singly charged electrons and ions')
    ME, MI = float(electron['m_particle']), float(ion['m_particle'])
    NI = float(ion['number_density_cm3'])*1e6 if 'number_density_cm3' in ion else float(ion['number_density_m3'])
    # Match BEACH's documented eV-to-K input conversion before SI normalization.
    temperature_k = float(electron['temperature_ev'])*1.160451812e4 if 'temperature_ev' in electron else float(electron['temperature_k'])
    TE = temperature_k*1.380649e-23/QE
    VE, VI = -float(electron['drift_velocity'][2]), -float(ion['drift_velocity'][2])
    if min(ME, MI, NI, TE, VI) <= 0.0:
        raise ValueError('masses, density, temperature, and incoming ion drift must be positive')
    VPE, VTE = math.sqrt(2.0*QE/ME), math.sqrt(2.0*QE*TE/ME)
    U = VE/VTE
    F0 = 0.5*math.erfc(-U)
    ELECTRON_FLUX_PER_DENSITY = VTE/(2.0*math.sqrt(math.pi))*(math.exp(-U*U)+2.0*math.sqrt(math.pi)*U*F0)
    ION_FLUX, ION_LIMIT = NI*VI, 0.5*MI*VI*VI/QE
    branch = model.get('zhao_branch', 'auto').removeprefix('zhao_')
    return dict(config_path=str(path.resolve()), config_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                branch=branch, all_configured_branches_covered=branch == 'b' or (branch == 'auto' and U > 0.0),
                branch_scope='Only Type B is numerically diagnosed. For u>0, A/C cannot connect to strict neutral field-free infinity.')


def load_spectrum(text):
    rows = list(csv.DictReader(io.StringIO(text)))
    if not rows:
        raise ValueError('empty spectrum CSV')
    def column(name):
        if name not in rows[0]:
            raise ValueError(f'missing {name}; columns={list(rows[0])}')
        return np.asarray([float(row[name].replace('D', 'E')) for row in rows])
    lo, hi, flux = [column(name) for name in ('energy_low_ev', 'energy_high_ev', 'flux_m2_s')]
    order = np.argsort(lo)
    lo, hi, flux = lo[order], hi[order], flux[order]
    assert np.all(np.isfinite(lo)) and np.all(np.isfinite(hi)) and np.all(np.isfinite(flux))
    assert np.all(lo >= 0.0) and np.all(hi > lo) and np.all(flux >= 0.0)
    assert np.all(lo[1:] >= hi[:-1] - 1e-12*np.maximum(1.0, hi[:-1]))
    return lo, hi, flux


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('spectrum', type=Path, help='stderr log with exact spectrum markers, or a spectrum CSV')
    parser.add_argument('--config', type=Path, required=True, help='the beach.toml used by the failed run')
    parser.add_argument('--output', type=Path, default=Path('scalar-b'))
    parser.add_argument('--query', type=Path, help='captured-query.txt supplies duration and D_before')
    parser.add_argument('--duration', type=float)
    parser.add_argument('--displacement-before', type=float)
    args = parser.parse_args()
    configuration = configure_from_toml(args.config)
    spectrum_text, query = read_failure_capture(args.spectrum)
    if args.query:
        for line in args.query.read_text().splitlines():
            if '=' in line:
                key, value = line.split('=', 1)
                query[key.strip()] = value.strip()
    if args.duration is None:
        if 'trial_batch_duration_s' not in query:
            parser.error('CSV input requires --query or explicit --duration and --displacement-before.')
        args.duration = float(query['trial_batch_duration_s'].replace('D', 'E'))
    if args.displacement_before is None:
        if 'displacement_before_c_m2' not in query:
            parser.error('CSV input requires --query or explicit --duration and --displacement-before.')
        args.displacement_before = float(query['displacement_before_c_m2'].replace('D', 'E'))
    assert args.duration > 0.0 and math.isfinite(args.displacement_before)
    lo, hi, flux = load_spectrum(spectrum_text)
    spectral_flux = flux/(hi-lo)
    total_flux = float(flux.sum())
    mean_energy = float(np.sum(flux*(lo+hi)/2.0)/total_flux) if total_flux > 0.0 else 0.0
    support_max = float(hi[flux > 0.0].max()) if np.any(flux > 0.0) else 0.0
    if 'guess' in query:
        feedback = [float(value.replace('D', 'E')) for value in query['guess'].split()]
        if not math.isclose(total_flux, feedback[0], rel_tol=1e-10, abs_tol=1e-20):
            raise ValueError('captured spectrum total flux disagrees with the failure query')
        if total_flux > 0.0 and not math.isclose(mean_energy, feedback[1], rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError('captured spectrum mean energy disagrees with the failure query')

    def tail(h):
        return float(np.sum(spectral_flux*np.maximum(hi-np.maximum(lo, h), 0.0)))

    def pe_density_infinity(h):
        return float(2.0/VPE*np.sum(spectral_flux*(np.sqrt(np.maximum(hi-h, 0.0)) - np.sqrt(np.maximum(lo-h, 0.0)))))

    def neutrality(h):
        return (NI-pe_density_infinity(h))/F0

    def upstream_sqrt_rho_coefficient(h, ne):
        # At phi=0+, Delta ne=-Ne exp(-u^2)/sqrt(pi Te)*sqrt(phi).
        # PE passing and returning legs contribute (-2 j_right+4 j_left)/VPE.
        # A positive rho coefficient makes E^2 negative arbitrarily near infinity.
        left_bin = (lo < h) & (hi >= h)
        right_bin = (lo <= h) & (hi > h)
        j_left = float(spectral_flux[left_bin].sum())
        j_right = float(spectral_flux[right_bin].sum())
        electron = ne*math.exp(-U*U)/math.sqrt(math.pi*TE)
        photo = (4.0*j_left-2.0*j_right)/VPE
        return electron-photo, abs(electron)+abs(photo)

    def current_number(h):
        ne = neutrality(h)
        return ION_FLUX+tail(h)-ne*ELECTRON_FLUX_PER_DENSITY

    @lru_cache(maxsize=131072)
    def electron_integral_unit(phi):
        """Integral_0^phi n_e/Ne dV, integrating upstream velocity first."""
        if phi == 0.0:
            return 0.0
        psi = phi/TE
        edge = math.sqrt(psi)
        top = max(12.0, U+12.0)
        value = quad(lambda a: a/(math.sqrt(a*a+psi)+a)*math.exp(-(a-U)**2),
                     0.0, top, points=sorted(set([min(edge, top/2), U])),
                     epsabs=2e-13, epsrel=2e-12, limit=250)[0]
        return 2.0*phi*value/math.sqrt(math.pi)

    def power_difference(x, phi):
        """(x+phi)_+^(3/2)-x_+^(3/2), stable for phi << positive x."""
        positive = np.maximum(x, 0.0)
        shifted = np.maximum(x+phi, 0.0)
        a, b = np.sqrt(shifted), np.sqrt(positive)
        difference = np.zeros_like(x)
        nonzero = a+b > 0.0
        difference[nonzero] = (shifted[nonzero]-positive[nonzero])*(
            shifted[nonzero]+a[nonzero]*b[nonzero]+positive[nonzero])/(a[nonzero]+b[nonzero])
        mask = x >= 0.0
        denominator = a[mask]+b[mask]
        difference[mask] = phi*np.divide(shifted[mask]+a[mask]*b[mask]+positive[mask], denominator,
                                          out=np.zeros_like(denominator), where=denominator > 0.0)
        return difference

    def pe_integral(phi, h):
        """Integral_0^phi PE density dV, outward plus reflected leg."""
        outward = power_difference(hi-h, phi)-power_difference(lo-h, phi)
        reflected_hi = np.minimum(hi, h)
        reflected_lo = np.minimum(lo, h)
        returning = power_difference(reflected_hi-h, phi)-power_difference(reflected_lo-h, phi)
        return float(4.0/(3.0*VPE)*np.sum(spectral_flux*(outward+returning)))

    def field_squared(phi, h, ne):
        if phi == 0.0:
            return 0.0
        ion_integral = 2.0*NI*phi/(1.0+math.sqrt(1.0-phi/ION_LIMIT))
        electron_integral = ne*electron_integral_unit(float(phi))
        photo_integral = pe_integral(phi, h)
        return 2.0*QE/EPS0*(electron_integral+photo_integral-ion_integral)

    @lru_cache(maxsize=131072)
    def endpoint(h):
        ne = neutrality(h)
        e2 = field_squared(h, h, ne)
        d = EPS0*math.sqrt(max(0.0, e2))
        residual = d-args.displacement_before-args.duration*QE*current_number(h)
        return ne, e2, d, residual

    def be_residual(h):
        ne, e2, d, residual = endpoint(float(h))
        if ne <= 0.0 or e2 < 0.0:
            return math.nan
        return residual

    # Cover the entire cold-ion-accessible open interval, resolving both the
    # compact source support and the cold-ion endpoint independently of guesses.
    max_h = float(np.nextafter(ION_LIMIT, 0.0))
    def make_grid(subdivisions, thermal_count):
        theta = np.linspace(0.0, 0.5*math.pi, thermal_count)
        thermal = np.minimum(ION_LIMIT*np.sin(theta)**2, max_h)
        bins = [lo+(hi-lo)*fraction for fraction in np.linspace(0.0, 1.0, subdivisions+1)]
        near_ion_end = ION_LIMIT*(1.0-10.0**(-np.arange(2.0, 15.0)))
        grid = np.unique(np.r_[thermal, *bins, near_ion_end, 0.0, max_h])
        return grid[(grid >= 0.0) & (grid <= max_h)]

    coarse_grid = make_grid(4, 513)
    grid = make_grid(16, 1025)
    current_scale = max(ION_FLUX, total_flux, NI*ELECTRON_FLUX_PER_DENSITY/F0)
    field_scale_squared = NI*QE*TE/EPS0

    # Add points where current, neutrality, or endpoint E^2 changes sign or bends
    # substantially. This is a finite adaptive search, not an interval proof.
    refinements = []
    def refine_interval(left, right, depth):
        middle = 0.5*left+0.5*right
        if middle == left or middle == right:
            return
        values = []
        for point in (left, middle, right):
            ne, e2, _, _ = endpoint(float(point))
            values.append(np.asarray([current_number(point)/current_scale, ne/NI, e2/field_scale_squared]))
        a, m, b = values
        scale = np.maximum.reduce([np.abs(a), np.abs(m), np.abs(b), np.full(3, 1e-8)])
        crossing = np.any(a*b < 0.0) or np.any(a*m < 0.0) or np.any(m*b < 0.0)
        curvature = np.any(np.abs(m-0.5*(a+b)) > 0.02*scale)
        near_zero = np.any(np.minimum.reduce([np.abs(a), np.abs(m), np.abs(b)]) < 0.05*scale)
        if not (crossing or curvature or near_zero):
            return
        refinements.append(middle)
        if depth < 4:
            refine_interval(left, middle, depth+1)
            refine_interval(middle, right, depth+1)
    for left, right in zip(grid[:-1], grid[1:]):
        refine_interval(float(left), float(right), 0)
    grid = np.unique(np.r_[grid, refinements])
    failed_brackets = []

    def roots_on_grid(function, points, label):
        roots = []
        prev_x = float(points[0])
        prev_y = function(prev_x)
        for x in points[1:]:
            x = float(x)
            y = function(x)
            if math.isfinite(prev_y) and prev_y == 0.0:
                roots.append(prev_x)
            if math.isfinite(prev_y) and math.isfinite(y) and prev_y*y < 0.0:
                try:
                    roots.append(brentq(function, prev_x, x, xtol=3e-13, rtol=1e-14))
                except ValueError as error:
                    failed_brackets.append(dict(function=label, left=prev_x, right=x, reason=str(error)))
            prev_x, prev_y = x, y
        if math.isfinite(prev_y) and prev_y == 0.0:
            roots.append(prev_x)
        unique = []
        for root in sorted(roots):
            if not unique or abs(root-unique[-1]) > 1e-9*max(1.0, abs(root)):
                unique.append(root)
        return unique

    coarse_current_roots = roots_on_grid(current_number, coarse_grid, 'coarse_current')
    current_roots = roots_on_grid(current_number, grid, 'current')
    neutrality_boundaries = roots_on_grid(neutrality, grid, 'neutrality')
    field_boundaries = roots_on_grid(lambda h: endpoint(float(h))[1], grid, 'endpoint_field_squared')
    # Explicitly insert domain boundaries and points on both sides before BE
    # root search, so an invalid midpoint is not silently bridged by a bracket.
    boundary_samples = []
    for value in neutrality_boundaries+field_boundaries:
        for offset in (-1e-8, -1e-10, -1e-12, 0.0, 1e-12, 1e-10, 1e-8):
            point = value+offset*max(1.0, value)
            if 0.0 <= point <= max_h:
                boundary_samples.append(point)
    grid = np.unique(np.r_[grid, boundary_samples])
    coarse_be_roots = roots_on_grid(be_residual, coarse_grid, 'coarse_BE')
    be_roots = roots_on_grid(be_residual, grid, 'BE')
    # At dt=2 s, D corrections can be much smaller than one bin: refine around
    # every current root even if the coarse BE scan never resolved the interval.
    current_neighborhoods = []
    for h in current_roots:
        points = [h]
        for radius in np.geomspace(1e-12, max(0.1, 0.05*h), 25):
            points.extend([max(0.0, h-radius), min(max_h, h+radius)])
        points = np.unique(points)
        current_neighborhoods.extend(points)
        be_roots.extend(roots_on_grid(be_residual, points, 'BE_near_current'))
    be_roots = sorted(be_roots)
    be_roots = [root for index, root in enumerate(be_roots)
                if index == 0 or abs(root-be_roots[index-1]) > 1e-9*max(1.0, abs(root))]
    grid = np.unique(np.r_[grid, current_neighborhoods])

    def density_rho(phi, h, ne):
        # Independent velocity integral for stationary points of E^2(phi).
        if phi == 0.0:
            return 0.0
        psi = phi/TE
        edge = math.sqrt(psi)
        ambient_unit = quad(lambda a: a/math.sqrt(a*a+psi)*math.exp(-(a-U)**2),
                            0.0, max(12.0, U+12.0), points=[min(edge, 6.0), U],
                            epsabs=2e-12, epsrel=2e-12, limit=250)[0]/math.sqrt(math.pi)
        shift = h-phi
        outward = np.sqrt(np.maximum(hi-shift, 0.0))-np.sqrt(np.maximum(lo-shift, 0.0))
        returning = np.sqrt(np.maximum(np.minimum(hi, h)-shift, 0.0))-np.sqrt(np.maximum(np.minimum(lo, h)-shift, 0.0))
        photo = 2.0/VPE*float(np.sum(spectral_flux*(outward+returning)))
        return NI/math.sqrt(1.0-phi/ION_LIMIT)-ne*ambient_unit-photo

    def inspect_profile(h, ne):
        turns = h-np.unique(np.r_[lo, hi])
        turns = turns[(turns >= 0.0) & (turns <= h)]
        knots = np.unique(np.r_[0.0, h, turns, h*np.geomspace(1e-10, 1e-2, 65), np.linspace(0.0, h, 257)])
        # Every PE turning-energy edge is represented; add intermediate values
        # and independently locate rho=0 extrema of the Sagdeev field.
        points = np.unique(np.r_[knots, 0.5*(knots[:-1]+knots[1:])])
        stationary = roots_on_grid(lambda phi: density_rho(phi, h, ne), points, 'profile_rho')
        points = np.unique(np.r_[points, stationary])
        values = np.asarray([field_squared(float(phi), h, ne) for phi in points])
        minima = []
        for index in range(1, len(points)-1):
            if values[index] < values[index-1] and values[index] < values[index+1]:
                result = minimize_scalar(lambda phi: field_squared(float(phi), h, ne),
                                         bounds=(points[index-1], points[index+1]), method='bounded',
                                         options={'xatol': 1e-12, 'maxiter': 120})
                minima.append(float(result.x))
        points = np.unique(np.r_[points, minima])
        values = np.asarray([field_squared(float(phi), h, ne) for phi in points])
        return points, values, stationary

    profile_rows = []
    candidates = []
    for kind, roots in [('current_balance', current_roots), ('backward_euler', be_roots)]:
        for index, h in enumerate(roots):
            ne, eh2, displacement, residual = endpoint(float(h))
            points, e2, stationary = inspect_profile(h, ne)
            tol = 1e-9*max(1.0, float(np.max(np.abs(e2))))
            upstream_coefficient, upstream_scale = upstream_sqrt_rho_coefficient(h, ne)
            upstream_obstructed = h > 0.0 and upstream_coefficient > 1e-10*max(1.0, upstream_scale)
            reasons = []
            if ne <= 0.0:
                reasons.append('nonpositive_neutral_electron_amplitude')
            if eh2 < 0.0:
                reasons.append('negative_endpoint_field_squared')
            if float(e2.min()) < -tol:
                reasons.append('negative_profile_field_squared')
            if upstream_obstructed:
                reasons.append('upstream_analytic_obstruction')
            if not all(math.isfinite(value) for value in (ne, eh2, float(e2.min()), upstream_coefficient)):
                reasons.append('nonfinite_evaluation')
            candidates.append(dict(kind=kind, index=index, phi_H=h, Ne=ne, E_H_squared=eh2,
                                   displacement=displacement if eh2 >= 0.0 else None,
                                   Gamma_e=ne*ELECTRON_FLUX_PER_DENSITY, Gamma_i=ION_FLUX,
                                   Gamma_PE_escape=tail(h), current_A_m2=QE*current_number(h),
                                   BE_residual_C_m2=residual if eh2 >= 0.0 else None,
                                   minimum_profile_E_squared=float(e2.min()),
                                   minimum_phi=float(points[int(e2.argmin())]),
                                   profile_tolerance_E_squared=tol, profile_points=len(points),
                                   profile_stationary_potentials=stationary,
                                   nearest_bin_edge_distance=float(np.min(np.abs(np.r_[lo, hi]-h))),
                                   upstream_sqrt_rho_coefficient_m3_Vmhalf=upstream_coefficient,
                                   upstream_analytic_obstruction=upstream_obstructed,
                                   rejection_reasons=reasons, physical=not reasons))
            for phi, value in zip(points, e2):
                profile_rows.append([kind, index, h, ne, float(phi), float(value)])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    Path(str(args.output)+'-spectrum.csv').write_text(spectrum_text.strip()+'\n')
    states = []
    with Path(str(args.output)+'-scan.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['phi_H', 'Ne', 'E_H_squared', 'current_A_m2', 'BE_residual_C_m2',
                         'upstream_sqrt_rho_coefficient', 'neutrality_valid', 'endpoint_field_valid'])
        for h in grid:
            ne, e2, d, residual = endpoint(float(h))
            coefficient, coefficient_scale = upstream_sqrt_rho_coefficient(h, ne)
            states.append((float(h), ne, e2, residual, coefficient))
            writer.writerow([h, ne, e2, QE*current_number(h), residual if ne > 0.0 and e2 >= 0.0 else '',
                             coefficient, ne > 0.0, e2 >= 0.0])
    with Path(str(args.output)+'-profiles.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['kind', 'index', 'phi_H', 'Ne', 'phi', 'E_squared'])
        writer.writerows(profile_rows)
    physical_roots = [c for c in candidates if c['kind']=='backward_euler' and c['physical']]
    valid_endpoint_states = [state for state in states if state[1] > 0.0 and state[2] >= 0.0]
    best = min(valid_endpoint_states, key=lambda state: abs(state[3])) if valid_endpoint_states else None
    root_stable = (len(coarse_current_roots) == len(current_roots) and
                   all(abs(a-b) <= 1e-7*max(1.0, abs(b)) for a, b in zip(coarse_current_roots, current_roots)))
    # Optional absence certificate using conservative binwise interval bounds.
    # This uses no root count or sampled profile to exclude an interval.
    def certify_absence():
        started = time.monotonic()
        epsilon_pad = 2048.0*np.finfo(float).eps
        exclusions = dict(nonpositive_Ne=0, upstream_analytic=0, positive_BE_residual=0, negative_BE_residual=0)
        unresolved = []
        known_root_intervals = []
        leaves = []
        checked = 0
        ceiling = min(ION_LIMIT, support_max)
        boundaries = np.unique(np.r_[0.0, lo[lo <= ceiling], hi[hi <= ceiling], ceiling])
        stack = [(float(left), float(right), 0) for left, right in zip(boundaries[:-1], boundaries[1:])]

        def bin_densities(h):
            a = np.maximum(lo-h, 0.0)
            b = np.maximum(hi-h, 0.0)
            denominator = np.sqrt(a)+np.sqrt(b)
            difference = np.divide(b-a, denominator, out=np.zeros_like(a), where=denominator > 0.0)
            return 2.0/VPE*spectral_flux*difference

        def photo_sqrt_coefficient(h):
            j_left = float(spectral_flux[(lo < h) & (hi >= h)].sum())
            j_right = float(spectral_flux[(lo <= h) & (hi > h)].sum())
            return (4.0*j_left-2.0*j_right)/VPE

        while stack:
            left, right, depth = stack.pop()
            checked += 1
            a, b = bin_densities(left), bin_densities(right)
            lower_pe = float(np.minimum(a, b).sum())
            upper_bins = np.maximum(a, b)
            peaks = (lo >= left) & (lo <= right)
            upper_bins[peaks] = np.maximum(upper_bins[peaks], 2.0/VPE*spectral_flux[peaks]*np.sqrt(hi[peaks]-lo[peaks]))
            upper_pe = float(upper_bins.sum())
            padding = epsilon_pad*max(NI, upper_pe, 1.0)
            lower_ne = (NI-upper_pe-padding)/F0
            upper_ne = (NI-lower_pe+padding)/F0
            reason = None
            if upper_ne <= 0.0:
                reason = 'nonpositive_Ne'
            else:
                lower_ne = max(0.0, lower_ne)
                # Include both interval endpoints: j(H-) and j(H+) differ at
                # exact bin boundaries, unlike the interior coefficient 2j/VPE.
                upper_photo_coefficient = max(photo_sqrt_coefficient(left), photo_sqrt_coefficient(right),
                                              photo_sqrt_coefficient(0.5*left+0.5*right))
                electron_coefficient = lower_ne*math.exp(-U*U)/math.sqrt(math.pi*TE)
                c_padding = epsilon_pad*max(1.0, abs(electron_coefficient), abs(upper_photo_coefficient))
                if left > 0.0 and electron_coefficient-upper_photo_coefficient > c_padding:
                    reason = 'upstream_analytic'
                else:
                    current_min = ION_FLUX+tail(right)-upper_ne*ELECTRON_FLUX_PER_DENSITY
                    current_max = ION_FLUX+tail(left)-lower_ne*ELECTRON_FLUX_PER_DENSITY
                    current_padding = epsilon_pad*max(current_scale, abs(current_min), abs(current_max))
                    current_min -= current_padding
                    current_max += current_padding
                    # n_e(phi)<=Ne*F0<=ni; omit the negative ion integral.
                    # The total PE endpoint integral increases monotonically in H.
                    photo_integral = pe_integral(right, right)
                    integral_upper = NI*right+photo_integral
                    integral_upper += epsilon_pad*max(1.0, abs(integral_upper), abs(photo_integral))
                    displacement_upper = math.sqrt(max(0.0, 2.0*EPS0*QE*integral_upper))
                    lower_f = -args.displacement_before-args.duration*QE*current_max
                    upper_f = displacement_upper-args.displacement_before-args.duration*QE*current_min
                    f_padding = epsilon_pad*max(abs(args.displacement_before), displacement_upper,
                                                args.duration*QE*current_scale, 1e-30)
                    if lower_f > f_padding:
                        reason = 'positive_BE_residual'
                    elif upper_f < -f_padding:
                        reason = 'negative_BE_residual'
            known_root_inside = any(left <= root['phi_H'] <= right for root in physical_roots)
            if reason:
                assert not known_root_inside, 'an interval bound excluded a known physical endpoint'
                exclusions[reason] += 1
                leaves.append([left, right, reason])
            elif known_root_inside:
                # Exercise the exclusion bounds first, then stop splitting an
                # interval already known to contain a physical counterexample.
                unresolved.append([left, right])
                known_root_intervals.append([left, right])
            elif depth >= 40 or right-left <= 1e-10*max(1.0, right) or checked >= 200000 or time.monotonic()-started > 60.0:
                unresolved.append([left, right])
                if checked >= 200000 or time.monotonic()-started > 60.0:
                    unresolved.extend([[a, b] for a, b, _ in stack])
                    stack.clear()
            else:
                middle = 0.5*left+0.5*right
                stack.extend([(left, middle, depth+1), (middle, right, depth+1)])
        # Zero support has no nonflat B state; the flat state remains a separate
        # endpoint and is checked from its exact neutrality and BE residual.
        ne_zero, _, _, f_zero = endpoint(0.0)
        zero_margin = epsilon_pad*max(abs(args.displacement_before), args.duration*QE*current_scale, 1e-30)
        flat_excluded = bool(ne_zero <= 0.0 or abs(f_zero) > zero_margin)
        unresolved_leaf_count = len(unresolved)
        merged_unresolved = []
        for left, right in sorted(unresolved):
            if merged_unresolved and left <= merged_unresolved[-1][1]:
                merged_unresolved[-1][1] = max(right, merged_unresolved[-1][1])
            else:
                merged_unresolved.append([left, right])
        unresolved = merged_unresolved
        complete = not unresolved and flat_excluded
        # A previously located physical root is an independent counterexample
        # to any false absence certificate; exercise this on the 3-root fixture.
        if physical_roots:
            assert not complete, 'interval exclusions contradicted a physical endpoint'
            for root in physical_roots:
                assert not any(left <= root['phi_H'] <= right for left, right, _ in leaves), \
                    'an interval bound excluded a known physical endpoint'
        with Path(str(args.output)+'-intervals.csv').open('w') as stream:
            writer = csv.writer(stream)
            writer.writerow(['phi_left', 'phi_right', 'exclusion_reason'])
            writer.writerows(sorted(leaves))
            writer.writerows([left, right, 'unresolved'] for left, right in unresolved)
            writer.writerows([left, right, 'contains_validated_endpoint'] for left, right in known_root_intervals)
        return dict(attempted=True, complete=complete, checked_intervals=checked, exclusions=exclusions,
                    unresolved=unresolved, unresolved_leaf_count=unresolved_leaf_count,
                    contains_validated_endpoint=known_root_intervals, flat_state_excluded=flat_excluded,
                    bounded_range=[0.0, ceiling], source_free_open_range_excluded=[ceiling, ION_LIMIT] if ceiling < ION_LIMIT else [],
                    elapsed_seconds=time.monotonic()-started, relative_roundoff_margin=epsilon_pad,
                    meaning='Binwise neutrality/current/field bounds plus the upstream analytic obstruction cover the interval; no root enumeration is used in these exclusions.',
                    numerical_assumption='IEEE float64 elementary arithmetic with the stated conservative roundoff margins; this is not directed-rounding interval arithmetic.')

    certificate = certify_absence()
    summary = dict(spectrum=str(args.spectrum.resolve()), sha256=hashlib.sha256(args.spectrum.read_bytes()).hexdigest(),
                   query=query, configuration=configuration,
                   ambient_outward_feedback_ignored_by_model=True,
                   constants=dict(ni=NI, Te=TE, me=ME, mi=MI, ve=VE, vi=VI, u=U, qe=QE, eps0=EPS0),
                   duration=args.duration, displacement_before=args.displacement_before,
                   total_flux=total_flux, mean_energy=mean_energy, bins=len(flux),
                   scan_phi_max=max_h, cold_ion_limit=ION_LIMIT, nonzero_spectrum_support_max=support_max,
                   analytic_support_exclusion='For phi_H > support_max and Ne > 0, the upstream sqrt(phi) rho coefficient is positive; a nonflat B profile is impossible.',
                   search=dict(coarse_points=len(coarse_grid), refined_points=len(grid),
                               adaptive_added_points=len(refinements), coarse_current_roots=coarse_current_roots,
                               refined_current_roots=current_roots, coarse_BE_roots=coarse_be_roots,
                               refined_BE_roots=be_roots, current_roots_stable_under_refinement=root_stable,
                               neutrality_boundaries=neutrality_boundaries, endpoint_field_boundaries=field_boundaries,
                               failed_brackets=failed_brackets,
                               sampled_nonpositive_Ne=sum(state[1] <= 0.0 for state in states),
                               sampled_negative_endpoint_E2=sum(state[2] < 0.0 for state in states),
                               closest_valid_sample_to_BE_zero=best),
                   candidates=candidates, physical_BE_candidates=len(physical_roots),
                   absence_certificate=certificate,
                   diagnosis=('physical_BE_endpoint_found' if physical_roots else
                              'type_B_absence_certified_in_model_with_roundoff_margin' if certificate['complete'] else
                              'located_BE_endpoints_rejected_by_physical_conditions' if be_roots else
                              'no_BE_endpoint_located_by_finite_search'),
                   limitation='The finite adaptive candidate/profile search is not a proof that all roots or all sign changes were located.')
    Path(str(args.output)+'.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
