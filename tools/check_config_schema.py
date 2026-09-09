"""Compare the public schema keys with the Fortran TOML dispatch tables.

This is a declaration check, not a substitute for the executable config contract
tests: values, defaults, and physical combinations are checked by those tests.
"""

from __future__ import annotations

import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
TABLES = (
    (None, "apply_toml_document"),
    ("sim", "apply_sim_toml_table"),
    ("domain", "apply_domain_toml_table"),
    ("fieldBoundary", "apply_field_boundary_toml_table"),
    ("particleBoundary", "apply_particle_boundary_toml_table"),
    ("reservoir", "apply_reservoir_toml_table"),
    ("surfaceCurrentModel", "apply_surface_current_model_toml_table"),
    ("periodic2", "apply_periodic2_toml_table"),
    ("particles", "apply_particles_toml_table"),
    ("species", "apply_particles_species_toml_table"),
    ("speciesParticleBoundary", "apply_species_boundary_toml_table"),
    ("speciesBoundaryInflow", "apply_species_boundary_inflow_toml_table"),
    ("mesh", "apply_mesh_toml_table"),
    ("meshGroup", "apply_mesh_group_toml_table"),
    ("template", "apply_template_toml_table"),
    ("output", "apply_output_toml_table"),
)


def dispatched_keys(source: str, procedure: str) -> set[str]:
    start = rf"^\s*(?:module procedure {procedure}\s*$|(?:module )?subroutine {procedure}\s*\()"
    end = rf"^\s*end (?:subroutine|procedure) {procedure}\b"
    bodies = re.findall(start + r"(.*?)" + end, source, re.I | re.M | re.S)
    for body in bodies:
        keys = set()
        depth = 0
        for line in body.splitlines():
            if re.match(r"\s*select case\b", line, re.I):
                depth += 1
            elif re.match(r"\s*end select\b", line, re.I):
                depth -= 1
            elif depth == 1:
                case = re.match(r"\s*case\s*\((.*?)\)", line, re.I)
                if case:
                    keys.update(re.findall(r"'([^']+)'", case.group(1)))
        if keys:
            return keys
    raise ValueError(f"no TOML key dispatch found for {procedure}")


def main() -> None:
    schema = json.loads((ROOT / "schemas/beach.schema.json").read_text())
    source = "\n".join(path.read_text() for path in sorted(
        (ROOT / "src/config/app_config_parser").glob("*.f90")
    ))
    failures = []
    for definition, procedure in TABLES:
        if definition is None:
            rule = schema
        elif definition == "periodic2":
            rule = schema["properties"][definition]
        else:
            rule = schema["$defs"][definition]
        expected = set(rule["properties"])
        actual = dispatched_keys(source, procedure)
        if actual != expected:
            failures.append(
                f"{definition or '<root>'}: schema-only={sorted(expected - actual)}, "
                f"Fortran-only={sorted(actual - expected)}"
            )
    if failures:
        raise SystemExit("Schema/Fortran key mismatch:\n" + "\n".join(failures))
    print(f"Schema/Fortran keys agree in {len(TABLES)} tables")


if __name__ == "__main__":
    main()
