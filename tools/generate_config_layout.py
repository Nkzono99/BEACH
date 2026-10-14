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
from beach.config._layout import CHOICES, PATHS, PRESERVED, SPECIES_PATHS, reader_schema


def object_rule():
    return {"type": "object", "additionalProperties": False, "properties": {}}


def insert(rule, path, value, *, structural=True):
    parts = path.split('.')
    for part in parts[:-1]:
        rule = rule.setdefault('properties', {}).setdefault(part, object_rule() if structural
                                             else {'type': 'object', 'properties': {}})
    rule.setdefault('properties', {})[parts[-1]] = copy.deepcopy(value)


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


def map_constraint(rule, paths, schema, choices=None):
    """Project reader constraints onto grouped paths without a public reader schema."""
    rule = dereference(rule, schema)
    choices = choices or {}
    result = {}
    for key, value in rule.items():
        if key == 'properties':
            result.setdefault('properties', {})
            for name, constraint in value.items():
                constraint = copy.deepcopy(constraint)
                if name in choices:
                    inverse = {old: new for new, old in choices[name].items()}
                    if 'const' in constraint:
                        constraint['const'] = inverse[constraint['const']]
                    if 'enum' in constraint:
                        constraint['enum'] = [inverse[value] for value in constraint['enum']]
                insert(result, paths[name], constraint, structural=False)
        elif key == 'required':
            for name in value:
                node = result
                parts = paths[name].split('.')
                for part in parts[:-1]:
                    required = node.setdefault('required', [])
                    if part not in required:
                        required.append(part)
                    node = node.setdefault('properties', {}).setdefault(part, {})
                node.setdefault('required', []).append(parts[-1])
        elif key in {'allOf', 'anyOf', 'oneOf'}:
            mapped = [map_constraint(item, paths, schema, choices) for item in value]
            if key == 'allOf':
                result.setdefault(key, []).extend(mapped)
            else:
                result[key] = mapped
        elif key in {'not', 'if', 'then', 'else'}:
            result[key] = map_constraint(value, paths, schema, choices)
        else:
            result[key] = copy.deepcopy(value)
    return result


def generate_schema(schema):
    reader = reader_schema(schema)
    grouped = object_rule()
    grouped['required'] = ['particles']
    for old,new in PATHS.items():
        try:
            rule = lookup(reader, old, schema)
        except KeyError:  # Online-only APIs are not advertised in released main.
            continue
        if new in CHOICES:
            rule = {'type': 'string', 'enum': list(CHOICES[new])}
        if new == 'fields.periodic.backend':
            rule = {'type':'string','enum':['cached_kneq0','panel_spectral_reference','finite_images'],
                    'description':'One owner of periodic field backend selection; split backends derive exclude_k0.'}
        insert(grouped,new,rule)
    for new,old in PRESERVED.items():
        insert(grouped,new,lookup(reader,old,schema))
    insert(grouped,'run.restart', {'type':'object','additionalProperties':False,
        'properties':{'from':lookup(reader,'output.restart_from',schema)},
        'description':'Presence resumes from a checkpoint. Omit from to search output.dir.'})
    species = object_rule()
    for old,new in SPECIES_PATHS.items():
        insert(species,new,schema['$defs']['species']['properties'][old])
    species['properties']['source']['required'] = ['mode']
    species['allOf'] = [{'if': {'not': {'required': ['source']}},
        'then': {'properties': {'sampling': {'properties': {
            'volume_macro_particles_per_batch': {'const': 0}}}}}}]
    species['allOf'].extend(map_constraint(rule, SPECIES_PATHS, schema)
                            for rule in schema['$defs']['species'].get('allOf', []))
    insert(grouped,'particles.species', {**copy.deepcopy(schema['$defs']['particles']['properties']['species']),
        'items':{'$ref':'#/$defs/groupedSpecies'}})
    distribution=species['properties']['distribution']
    distribution['allOf']=[{'not':{'required':list(pair)}} for pair in
                          [('number_density_cm3','number_density_m3')]]
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
    grouped['allOf'] = [map_constraint(rule, {
        key: PATHS['sim.' + key] for key in schema['$defs']['sim']['properties']
    }, schema) for rule in schema['$defs']['sim'].get('allOf', [])]
    # Restart table presence is the public representation of the reader resume flag.
    grouped['allOf'].append({
        'if': {'required': ['run'], 'properties': {'run': {'required': ['restart']}}},
        'then': {'properties': {'output': {'properties': {'enabled': {'const': True}}}}},
    })
    for definition, container, paths in [
        ('surfaceCurrentModel', 'sheath', {key: PATHS['surface_current_model.' + key].removeprefix('sheath.')
            for key in schema['$defs']['surfaceCurrentModel']['properties']}),
        ('mesh', 'mesh', {key: PATHS.get('mesh.' + key, 'mesh.' + key).removeprefix('mesh.')
            for key in schema['$defs']['mesh']['properties']}),
    ]:
        target = lookup(grouped, container, schema)
        for composition in ('allOf', 'anyOf', 'oneOf'):
            if composition in schema['$defs'][definition]:
                target[composition] = [map_constraint(rule, paths, schema,
                    {"model": CHOICES["sheath.closure"]} if definition == "surfaceCurrentModel" else {})
                    for rule in schema['$defs'][definition][composition]]
    definitions={**schema['$defs'], 'groupedConfig':grouped, 'groupedSpecies':species}
    return {'$schema':schema['$schema'],'title':'BEACH Parameter File',
        'description':'Grouped BEACH input (run/domain/mesh/particles/fields/sheath/output). Runtime preflight additionally checks physical combinations.',
        '$ref':'#/$defs/groupedConfig', '$defs':definitions}


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
           '  public :: reader_layout_path, grouped_layout_table','contains',
           '  !> Return the existing reader path for one grouped setting.',
           '  function reader_layout_path(path) result(reader_path)',
           '    character(len=*), intent(in) :: path',
           '    character(len=:), allocatable :: reader_path',
           '    select case (path)']
    for new,old in paths.items():
        lines.extend([f"    case ('{new}')",f"      reader_path = '{old}'"])
    lines.extend(['    case default',"      reader_path = ''",'    end select',
                  '  end function reader_layout_path', '',
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
