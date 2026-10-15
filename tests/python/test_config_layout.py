from __future__ import annotations

import copy
from pathlib import Path
import subprocess
import tomllib

import pytest

from beach.config import ConfigError, default_config, dump_beach_toml, load_config_file, normalize_config_document
from beach.config._layout import to_grouped_layout, to_runtime_layout

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = sorted(ROOT.joinpath('examples').glob('*.toml'))


@pytest.mark.parametrize('path', EXAMPLES, ids=lambda path: path.name)
def test_grouped_example_roundtrip_keeps_runtime_dictionary(path, tmp_path):
    before = load_config_file(path)
    destination = tmp_path / path.name
    destination.write_text(dump_beach_toml(before), encoding='utf-8')
    assert load_config_file(destination) == before
    assert set(tomllib.loads(destination.read_text())) <= {
        'run', 'domain', 'mesh', 'particles', 'fields', 'sheath', 'output'
    }


def test_init_emits_grouped_input_without_changing_runtime_defaults():
    runtime = default_config()
    grouped = tomllib.loads(dump_beach_toml(runtime))
    assert 'sim' not in grouped
    assert grouped['particles']['tracking']['dt_s'] == runtime['sim']['dt']
    assert normalize_config_document(grouped) == runtime


@pytest.mark.parametrize('runtime', [default_config(), {'sim': {'dt': 1e-9}, 'particles': {'species': [{}]}}])
def test_flat_input_is_rejected(runtime, tmp_path):
    import tomli_w

    with pytest.raises(ConfigError, match='flat beach.toml layout is not supported'):
        normalize_config_document(runtime)
    with pytest.raises(ConfigError, match='flat beach.toml layout is not supported'):
        to_runtime_layout(runtime)
    path = tmp_path / 'flat.toml'
    path.write_text(tomli_w.dumps(runtime))
    with pytest.raises(ConfigError, match='flat beach.toml layout is not supported'):
        load_config_file(path)


def test_boundary_only_species_has_no_standalone_source():
    runtime = default_config()
    runtime['sim']['batch_duration'] = 1e-6
    item = runtime['particles']['species'][0]
    for key in ('pos_low', 'pos_high', 'drift_velocity', 'temperature_k'):
        item.pop(key, None)
    item.update(source_mode='volume_seed', npcls_per_step=0, number_density_m3=1e6,
                temperature_ev=1.0, boundary_inflow={'z_high': 'reservoir'})
    grouped = to_grouped_layout(runtime)
    assert 'source' not in grouped['particles']['species'][0]
    assert 'volume_macro_particles_per_batch' not in grouped['particles']['species'][0]['sampling']
    assert normalize_config_document(grouped)['particles']['species'][0]['boundary_inflow'] == {'z_high': 'reservoir'}


def test_outflow_refresh_uses_the_grouped_coupling_table():
    grouped = tomllib.loads((ROOT / 'examples/zero_current_refresh.toml').read_text())
    assert grouped['sheath']['coupling'] == {'outflow_refresh_batches': 2}
    assert grouped['sheath']['photoelectrons'] == {'solar_elevation_deg': 60.0, 'ref_density_m3': 6.4e7}
    assert to_runtime_layout(grouped)['surface_current_model']['outflow_refresh_batches'] == 2


@pytest.mark.parametrize('mutation', [
    lambda c: c.update(sim={'dt': 1e-9}),
    lambda c: c['particles']['species'][0].update(q_particle=-1e-19),
    lambda c: c['run'].update(unknown=1),
    lambda c: c['particles']['species'][0].update(source={}),
    lambda c: c.setdefault('fields', {}).update(periodic={'backend': 'finite_images', 'lower_boundary_model': 'e_bottom_zero'}),
    lambda c: c['particles']['tracking'].update(dt_s=float('nan')),
])
def test_grouped_rejects_unknown_mixed_or_invalid_values(mutation):
    config = to_grouped_layout(default_config())
    mutation(config)
    original = copy.deepcopy(config)
    with pytest.raises(ConfigError):
        normalize_config_document(config)
    assert repr(config) == repr(original)


def test_restart_presence_controls_internal_resume_flag():
    grouped = to_grouped_layout(default_config())
    grouped['run']['restart'] = {}
    assert to_runtime_layout(grouped)['output']['resume'] is True
    runtime = default_config()
    runtime['output']['resume'] = False
    assert 'restart' not in to_grouped_layout(runtime)['run']


@pytest.mark.parametrize('path', EXAMPLES, ids=lambda path: path.name)
def test_grouped_examples_pass_native_config_check(path, native_checker):
    result = subprocess.run([str(native_checker), '--check-config', str(path)],
                            cwd=ROOT, text=True, capture_output=True, timeout=20)
    assert result.returncode == 0, result.stdout + result.stderr


def test_native_rejects_flat_input(native_checker, tmp_path):
    import tomli_w

    path = tmp_path / 'flat.toml'
    path.write_text(tomli_w.dumps(default_config()))
    result = subprocess.run([str(native_checker), '--check-config', str(path)],
                            cwd=ROOT, text=True, capture_output=True, timeout=20)
    assert result.returncode != 0
    assert 'flat beach.toml layout is not supported' in result.stdout + result.stderr


def test_secondary_emission_layout_declares_public_keys_and_modes():
    runtime = load_config_file(ROOT / "examples/zero_current_secondary.toml")
    incident = runtime["particles"]["species"][0]
    target = runtime["particles"]["species"][-1]
    assert incident["secondary_emission_enabled"] is True
    assert incident["secondary_emission_target_species"] == target["species_key"]
    assert target["source_mode"] == "secondary_emission"
    assert runtime["surface_current_model"]["secondary_electron_species"] == target["species_key"]
    grouped = to_grouped_layout(runtime)
    assert grouped["particles"]["species"][0]["surface_interaction"]["secondary_emission"]["model"] == "whipple"
    grouped["particles"]["secondary_emission"] = {"max_generations": 2}
    assert to_runtime_layout(grouped)["particles"]["secondary_emission_max_generations"] == 2
