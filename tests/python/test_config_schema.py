from __future__ import annotations

import copy
from pathlib import Path

import pytest

from beach.config import default_config, to_grouped_layout
from beach.config._toml import load_toml_file
from beach.config.schema import load_schema, schema_errors

ROOT = Path(__file__).resolve().parents[2]


def _example(name):
    return load_toml_file(ROOT / 'examples' / name)


def test_public_schema_accepts_only_grouped_layout():
    schema, _ = load_schema()
    assert schema['$ref'] == '#/$defs/groupedConfig'
    assert 'legacyConfig' not in schema['$defs']
    assert schema_errors(default_config(), schema)
    assert schema_errors(to_grouped_layout(default_config()), schema) == []


@pytest.mark.parametrize('path', sorted(ROOT.joinpath('examples').glob('*.toml')), ids=lambda p: p.name)
def test_schema_accepts_bundled_grouped_examples(path):
    schema, _ = load_schema()
    assert schema_errors(load_toml_file(path), schema) == []


def test_schema_distinguishes_fortran_integer_storage_from_real_values():
    schema, _ = load_schema()
    config = to_grouped_layout(default_config())
    for value in (-(2**31), 2**31 - 1):
        config['run']['rng_seed'] = value
        assert schema_errors(config, schema) == []
    for value in (-(2**31) - 1, 2**31, 1.0, True):
        config['run']['rng_seed'] = value
        assert schema_errors(config, schema)
    config = to_grouped_layout(default_config())
    config['particles']['species'][0]['distribution']['number_density_m3'] = 10**12
    assert schema_errors(config, schema) == []


def test_tree_parameter_schema_uses_runtime_omission_semantics():
    schema, _ = load_schema()
    solver = schema['$defs']['groupedConfig']['properties']['fields']['properties']['solver']['properties']
    for key in ('theta', 'leaf_max'):
        assert 'default' not in solver['tree']['properties'][key]
        assert 'element-count heuristic' in solver['tree']['properties'][key]['description']
    assert 'direct evaluation to FMM' in solver['auto_fmm_min_elements']['description']


@pytest.mark.parametrize('key', ['field', 'outer_plasma', 'coupling', 'external_boundary', 'schema_version'])
def test_schema_rejects_unknown_top_level_keys(key):
    schema, _ = load_schema()
    config = to_grouped_layout(default_config())
    config[key] = {}
    assert schema_errors(config, schema)


@pytest.mark.parametrize('legacy_mesh', [
    {'mode': 'template'}, {'obj': {'path': 'unused.obj'}}, {'templates': [{'kind': 'plane'}]},
])
def test_schema_rejects_removed_mesh_keys(legacy_mesh):
    # The Fortran loader's rejection is test_app_config_parser mesh_parts_public_contract.
    schema, _ = load_schema()
    config = to_grouped_layout(default_config())
    config['mesh'] = legacy_mesh
    assert schema_errors(config, schema)


def test_schema_enforces_soft_discard_bounds_only_when_active():
    schema, _ = load_schema()
    config = to_grouped_layout(default_config())
    events = config['particles']['tracking'].setdefault('events', {})
    events.update(discard_count_grace=-1, discard_fraction_limit=-1.0, discard_charge_warning_c=-1.0)
    assert schema_errors(config, schema) == []
    events['policy'] = 'abort'
    assert schema_errors(config, schema) == []
    events.update(policy='soft_discard', discard_count_grace=0, discard_fraction_limit=1.0, discard_charge_warning_c=1e-30)
    assert schema_errors(config, schema) == []
    for key, value in [('discard_count_grace', -1), ('discard_count_grace', 0.5),
                       ('discard_fraction_limit', 0.0), ('discard_fraction_limit', 1.000001),
                       ('discard_fraction_limit', '1e-6'), ('discard_charge_warning_c', 0.0)]:
        invalid = copy.deepcopy(config)
        invalid['particles']['tracking']['events'][key] = value
        assert schema_errors(invalid, schema)


@pytest.mark.parametrize('action,valid', [('reflect', True), ('redistributed_reflect', True), ('unknown', False)])
def test_schema_species_boundary_actions(action, valid):
    schema, _ = load_schema()
    config = _example('beach.toml')
    config['particles']['species'][0]['boundary'] = {'z_high': action}
    assert (schema_errors(config, schema) == []) is valid


def test_schema_checks_source_and_current_channels():
    schema, _ = load_schema()
    config = _example('closed_photoelectron.toml')
    for path, value in [('charging.closure', 'unknown'), ('source.mode', 'volume_seed'),
                        ('source.deposit_opposite_charge_on_emit', False), ('charge_c', 1e-19)]:
        invalid = copy.deepcopy(config)
        node = invalid['particles']['species'][-1]
        parts = path.split('.')
        for part in parts[:-1]:
            node = node[part]
        node[parts[-1]] = value
        assert schema_errors(invalid, schema)
    species = config['particles']['species'][-1]
    species['charging'].update(closure='fixed_current', target_absorbed_current_a=-2e-6, target_emission_current_a=3e-6)
    assert schema_errors(config, schema) == []
    del species['charging']['closure']
    assert schema_errors(config, schema)


@pytest.mark.parametrize('key,values', [
    ('upstream_band_tolerance', [(0.2, True), (-0.1, False), (1.0, False)]),
    ('upstream_field_floor_v_m', [(2e-3, True), (0.0, True), (-1e-3, False)]),
    ('ion_pressure_factor', [(3.0, True), (0.0, False), (-1.0, False)]),
])
def test_schema_zhao_ranges(key, values):
    schema, _ = load_schema()
    for value, valid in values:
        config = _example('zero_current_refresh.toml')
        config['sheath']['zhao'][key] = value
        assert (schema_errors(config, schema) == []) is valid


def test_schema_checks_reflected_electron_switch():
    schema, _ = load_schema()
    config = _example('zero_current_refresh.toml')
    config['sheath']['coupling']['include_reflected_electrons'] = True
    assert schema_errors(config, schema) == []
    assert load_schema()[0]['$defs']['surfaceCurrentModel']['properties']['include_reflected_electrons']['default'] is True
    config['sheath']['coupling']['include_reflected_electrons'] = 'yes'
    assert schema_errors(config, schema)
    disabled = _example('zero_current.toml')
    disabled['sheath'] = {'closure': 'none', 'coupling': {'include_reflected_electrons': True}}
    assert schema_errors(disabled, schema)


def test_schema_checks_zhao_closure_and_refresh_keys():
    schema, _ = load_schema()
    config = _example('zero_current_refresh.toml')
    config['sheath']['coupling']['outflow_refresh_batches'] = -1
    assert schema_errors(config, schema)
    for key in ('response_backend', 'response_table_path', 'implicit_zero_mode', 'coupling_rtol'):
        invalid = _example('zero_current.toml')
        invalid['sheath'][key] = 'removed'
        assert schema_errors(invalid, schema)
    disabled = _example('zero_current.toml')
    disabled['sheath'] = {'closure': 'none'}
    assert schema_errors(disabled, schema) == []
    disabled['sheath']['zhao'] = {'branch': 'a'}
    assert schema_errors(disabled, schema)
    no_photo = _example('zero_current_no_photo.toml')
    no_photo['sheath']['zhao']['branch'] = 'c'
    assert schema_errors(no_photo, schema) == []
    no_photo['sheath']['zhao']['branch'] = 'a'
    assert schema_errors(no_photo, schema)


def test_schema_requires_kind_and_describes_all_part_key_types():
    schema, _ = load_schema()
    config = _example('beach.toml')
    config['mesh']['parts'].append({'enabled': False})
    assert schema_errors(config, schema)
    config['mesh']['parts'][-1]['kind'] = 'obj'
    assert schema_errors(config, schema) == []
    config['mesh']['parts'][-1]['rotation'] = [0, 1]
    assert schema_errors(config, schema)


@pytest.mark.parametrize('height', ['off', 'auto', 5e-4])
def test_upper_fourier_schema_types(height):
    schema, _ = load_schema()
    config = _example('zero_current_refresh.toml')
    config['fields']['periodic']['upper_fourier'] = {'switch_height': height, 'mode_layers': 7}
    assert schema_errors(config, schema) == []
    for invalid in ('on', True, [], {}):
        config['fields']['periodic']['upper_fourier']['switch_height'] = invalid
        assert schema_errors(config, schema)
    config['fields']['periodic']['upper_fourier']['switch_height'] = height
    for invalid in (-1, 33, 2.5, True):
        config['fields']['periodic']['upper_fourier']['mode_layers'] = invalid
        assert schema_errors(config, schema)
