"""Generate grouped schema and Fortran paths from the Python layout declaration.

Run after editing beach/config/_layout.py. The scalar reader and its preflight
checks remain the owners of units, defaults, and physical constraints.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from beach.config._layout import CHOICES, PATHS, PRESERVED, SPECIES_PATHS


def object_rule():
    return {"type": "object", "additionalProperties": False, "properties": {}}


def insert(rule, path, value):
    parts = path.split('.')
    for part in parts[:-1]:
        rule = rule['properties'].setdefault(part, object_rule())
    rule['properties'][parts[-1]] = copy.deepcopy(value)


def dereference(rule, schema):
    while '$ref' in rule:
        node = schema
        for segment in rule['$ref'][2:].split('/'):
            node = node[segment]
        rule = node
    return rule


def lookup(rule, path, schema):
    for segment in path.split('.'):
        rule = dereference(rule, schema)['properties'][segment]
    return rule


def generate_schema(schema):
    legacy = schema.get('$defs', {}).get('legacyConfig')
    if legacy is None:
        legacy = {k: copy.deepcopy(v) for k,v in schema.items()
                  if k in {'type', 'properties', 'required', 'additionalProperties', 'allOf'}}
    grouped = object_rule()
    grouped['required'] = ['particles']
    for old,new in PATHS.items():
        try:
            rule = lookup(legacy, old, schema)
        except KeyError:  # Online-only APIs are not advertised in released main.
            continue
        if new in CHOICES:
            rule = {'type': 'string', 'enum': list(CHOICES[new])}
        if new == 'fields.periodic.backend':
            rule = {'type':'string','enum':['cached_kneq0','panel_spectral_reference','finite_images'],
                    'description':'One owner of periodic field backend selection; split backends derive exclude_k0.'}
        insert(grouped,new,rule)
    for new,old in PRESERVED.items():
        insert(grouped,new,lookup(legacy,old,schema))
    insert(grouped,'run.restart', {'type':'object','additionalProperties':False,
        'properties':{'from':lookup(legacy,'output.restart_from',schema)},
        'description':'Presence resumes from a checkpoint. Omit from to search output.dir.'})
    species = object_rule()
    for old,new in SPECIES_PATHS.items():
        insert(species,new,schema['$defs']['species']['properties'][old])
    species['properties']['source']['required'] = ['mode']
    species['allOf'] = [{'if': {'not': {'required': ['source']}},
        'then': {'properties': {'sampling': {'properties': {
            'volume_macro_particles_per_batch': {'const': 0}}}}}}]
    insert(grouped,'particles.species', {**copy.deepcopy(schema['$defs']['particles']['properties']['species']),
        'items':{'$ref':'#/$defs/groupedSpecies'}})
    # The same mutual exclusions apply before and after authoring normalization.
    for container,pairs in [('run.batch',[('duration_s','duration_steps')]),
                            ('particles.tracking',[]),
                            ('fields.external',[('electric_v_m','electric_magnitude_v_m')])]:
        rule = lookup(grouped,container,schema)
        if pairs:
            rule['allOf'] = [{'not':{'required':list(pair)}} for pair in pairs]
    distribution=species['properties']['distribution']
    distribution['allOf']=[{'not':{'required':list(pair)}} for pair in
                          [('number_density_cm3','number_density_m3'),('temperature_ev','temperature_k')]]
    descriptions={
        'run':'Batch clock, stopping target, random seed, and restart.',
        'domain':'Finite box geometry and periodic topology.',
        'mesh':'Surface geometry, OBJ input, groups, and templates.',
        'particles':'Boundary actions, reservoir, tracking, species and sampling.',
        'fields':'Field boundary, imposed fields, solver and periodic backend.',
        'sheath':'External sheath zero-current closure, roles, photoelectron source and outer-root refresh.',
        'output':'Files, history, checkpoints and diagnostics.'}
    for key,description in descriptions.items():
        grouped['properties'][key]['description']=description
    definitions={**schema['$defs'], 'legacyConfig':legacy, 'groupedConfig':grouped, 'groupedSpecies':species}
    return {'$schema':schema['$schema'],'title':'BEACH Parameter File',
        'description':'Grouped BEACH authoring input. Released flat input is readable during 1.x; removal is planned for 2.0. Runtime preflight additionally checks physical combinations.',
        'anyOf':[{'$ref':'#/$defs/groupedConfig'},{'$ref':'#/$defs/legacyConfig'}], '$defs':definitions}


def fortran_paths():
    paths={new:old for old,new in PATHS.items()}
    paths.update(PRESERVED)
    paths.update({f'particles.species.{new}':f'particles.species.{old}'
                  for old,new in SPECIES_PATHS.items()})
    preserved={'domain','mesh.groups','mesh.templates','particles.species.boundary','particles.species.inflow'}
    containers={'', 'run.restart','particles.species'}
    for path in paths:
        parts=path.split('.')
        containers.update('.'.join(parts[:i]) for i in range(1,len(parts)))
    lines=['! Generated by tools/generate_config_layout.py; edit beach/config/_layout.py.',
           'module bem_config_layout_paths','  implicit none','  private',
           '  public :: legacy_layout_path, grouped_layout_table','contains',
           '  !> Return the existing reader path for one grouped setting.',
           '  function legacy_layout_path(path) result(legacy)',
           '    character(len=*), intent(in) :: path',
           '    character(len=:), allocatable :: legacy',
           '    select case (path)']
    for new,old in paths.items():
        lines.extend([f"    case ('{new}')",f"      legacy = '{old}'"])
    lines.extend(['    case default',"      legacy = ''",'    end select',
                  '  end function legacy_layout_path', '',
                  '  !> Recognize structural tables, rejecting unknown empty groups.',
                  '  logical function grouped_layout_table(path) result(known)',
                  '    character(len=*), intent(in) :: path','    select case (path)'])
    for path in sorted(containers-preserved):
        lines.extend([f"    case ('{path}')",'      known = .true.'])
    lines.extend(['    case default','      known = .false.','    end select',
                  '  end function grouped_layout_table','end module bem_config_layout_paths',''])
    return '\n'.join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args=parser.parse_args()
    schema=generate_schema(json.loads((ROOT/'schemas/beach.schema.json').read_text()))
    contents=json.dumps(schema,ensure_ascii=False,indent=2)+'\n'
    outputs={ROOT/path:contents for path in ['schemas/beach.schema.json',
        'beach/config/schemas/beach.schema.json','plugins/beach-context/references/schemas/beach.schema.json']}
    outputs[ROOT/'src/config/bem_config_layout_paths.f90']=fortran_paths()
    stale=[]
    for path,text in outputs.items():
        if args.check:
            if not path.exists() or path.read_text()!=text:
                stale.append(str(path.relative_to(ROOT)))
        else:
            path.write_text(text)
    if stale:
        raise SystemExit('Run tools/generate_config_layout.py: '+', '.join(stale))

if __name__=='__main__':
    main()
