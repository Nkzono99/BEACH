"""Grouped input preserves the physical runtime contract and has one key owner."""
from pathlib import Path
import copy
import tomllib

import pytest

from beach.config import default_config, dump_beach_toml, normalize_config_document
from beach.config._layout import to_grouped_layout, to_runtime_layout
from beach.config._shared import ConfigError

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize('name', ['beach.toml', 'periodic2_zhao_fixed_current.toml',
                                 'periodic2_zhao_outflow_refresh.toml',
                                 'periodic2_matching_plane_zhao_online.toml'])
def test_released_input_roundtrip_keeps_physics(name):
    path=ROOT/'examples'/name
    if name=='beach.toml':
        path=ROOT/'tests/fixtures/config_layout/beach_legacy.toml'
    if not path.exists():
        path=ROOT/'examples/periodic2_matching_plane_quasistatic.toml'
    old=tomllib.loads(path.read_text())
    grouped=tomllib.loads(dump_beach_toml(old))
    assert set(grouped) <= {'run','domain','mesh','particles','fields','sheath','output'}
    before=normalize_config_document(old)
    after=normalize_config_document(grouped)
    # Grouped split backends explicitly materialize the old effective k=0 policy.
    if 'periodic2' in after and 'periodic2' not in before:
        before['periodic2']={k:v for k,v in after['periodic2'].items()}
    before.get('periodic2',{}).setdefault('zero_mode_policy','exclude_k0') if 'periodic2' in before else None
    for item, new_item in zip(before['particles']['species'], after['particles']['species']):
        if item.get('source_mode', 'volume_seed') == 'volume_seed' and item.get('npcls_per_step', 0) == 0:
            for key in ('source_mode', 'npcls_per_step', 'pos_low', 'pos_high'):
                if key not in new_item:
                    item.pop(key, None)
    assert before == after


def test_init_emits_grouped_input_without_changing_runtime_defaults():
    old=default_config()
    grouped=tomllib.loads(dump_beach_toml(old))
    assert 'sim' not in grouped
    assert grouped['particles']['tracking']['dt_s']==old['sim']['dt']
    assert normalize_config_document(grouped)==normalize_config_document(old)


def test_boundary_only_species_has_no_standalone_source():
    old=default_config()
    old['sim']['batch_duration']=1e-6
    item=old['particles']['species'][0]
    for key in ('pos_low','pos_high','drift_velocity','temperature_k'):
        item.pop(key,None)
    item.update(source_mode='volume_seed',npcls_per_step=0, number_density_m3=1e6,
                temperature_ev=1.0,boundary_inflow={'z_high':'reservoir'})
    grouped=to_grouped_layout(old)
    assert 'source' not in grouped['particles']['species'][0]
    assert 'volume_macro_particles_per_batch' not in grouped['particles']['species'][0]['sampling']
    resolved=normalize_config_document(grouped)
    assert resolved['particles']['species'][0]['boundary_inflow']=={'z_high':'reservoir'}


def test_named_atol_is_independent_of_order_and_defaults_missing_components_to_zero():
    grouped=to_grouped_layout(default_config())
    grouped['sheath']={'closure':'none','coupling':{'atol':{
        'ion_outward_flux_m2_s':4.0,'photoelectron_outward_flux_m2_s':1.0}}}
    flat=to_runtime_layout(grouped)
    assert flat['surface_current_model']['coupling_atol']==[1.0,0.0,0.0,4.0]


@pytest.mark.parametrize('mutation', [
    lambda c:c.update(sim={'dt':1e-9}),
    lambda c:c['particles']['species'][0].update(q_particle=-1e-19),
    lambda c:c['run'].update(unknown=1),
    lambda c:c['particles']['species'][0].update(source={}),
    lambda c:c.setdefault('fields',{}).update(periodic={'backend':'finite_images','lower_boundary_model':'e_bottom_zero'}),
    lambda c:c['particles']['tracking'].update(dt_s=float('nan')),
])
def test_grouped_rejects_unknown_mixed_or_invalid_values(mutation):
    config=to_grouped_layout(default_config())
    mutation(config)
    original=copy.deepcopy(config)
    with pytest.raises(ConfigError):
        normalize_config_document(config)
    assert repr(config)==repr(original)


def test_restart_presence_and_explicit_false_legacy_are_distinct():
    grouped=to_grouped_layout(default_config())
    grouped['run']['restart']={}
    assert to_runtime_layout(grouped)['output']['resume'] is True
    old=default_config()
    old['output']['resume']=False
    assert 'restart' not in to_grouped_layout(old)['run']


def test_migrate_keeps_response_target_when_output_directory_changes(tmp_path):
    from beach.cli.main import main
    old=tomllib.loads((ROOT/'examples/periodic2_matching_plane_quasistatic.toml').read_text())
    old['surface_current_model']['response_table_path']='response.csv'
    source=tmp_path/'source/input.toml'
    source.parent.mkdir()
    import tomli_w
    source.write_text(tomli_w.dumps(old))
    destination=tmp_path/'dest/grouped.toml'
    main(['config','migrate',str(source),str(destination)])
    grouped=tomllib.loads(destination.read_text())
    assert grouped['sheath']['table']['path']==str(source.parent/'response.csv')
    with pytest.raises(SystemExit, match='already exists'):
        main(['config','migrate',str(source),str(destination)])


def test_grouped_native_loading_uses_the_same_config_contract(tmp_path):
    import os
    import subprocess
    executable=os.environ.get('BEACH_CONFIG_CHECK_EXE')
    if not executable:
        pytest.skip('native config check executable required')
    for old_path in ROOT.joinpath('tests/fixtures/config_layout').glob('*_legacy.toml'):
        grouped=tmp_path/(old_path.stem+'.toml')
        grouped.write_text(dump_beach_toml(tomllib.loads(old_path.read_text())))
        result=subprocess.run([executable,'--check-config',str(grouped)],
            cwd=ROOT, text=True, capture_output=True, timeout=20)
        assert result.returncode==0, result.stdout+result.stderr
    for grouped in ROOT.joinpath('examples/grouped').glob('*.toml'):
        result=subprocess.run([executable,'--check-config',str(grouped)],
            cwd=ROOT, text=True, capture_output=True, timeout=20)
        assert result.returncode==0, result.stdout+result.stderr
