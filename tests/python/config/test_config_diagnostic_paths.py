"""Keep internal reader names out of Fortran configuration diagnostics."""

from pathlib import Path
import re

import pytest

from beach.config._layout import PATHS, SPECIES_PATHS


ROOT = Path(__file__).resolve().parents[3]
LAYOUT_FILES = {"bem_config_layout.f90", "bem_config_layout_paths.f90"}
# Match comments after strings so ! inside a literal remains part of the literal.
TOKENS = re.compile(r"'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\"|![^\n]*")
READER_TABLES = ("sim.", "field_boundary.", "particle_boundary.",
                 "surface_current_model.", "periodic2.", "reservoir.")
REMOVED_MESH_NAMES = {"mesh.mode", "mesh.obj", "mesh.templates", "mesh.obj_path",
                      "mesh.obj_scale", "mesh.obj_rotation", "mesh.obj_offset"}
TABLE_NAMES = {reader for reader, public in PATHS.items()
               if reader != public and reader.startswith(READER_TABLES)}
SPECIES_NAMES = {reader for reader, public in SPECIES_PATHS.items() if reader != public}
RENAMED_LEAVES = {reader for reader in SPECIES_NAMES
                  if reader != SPECIES_PATHS[reader].rsplit(".", 1)[-1]}


def _alternatives(names: set[str]) -> str:
    return "|".join(re.escape(name) for name in sorted(names, key=len, reverse=True))


# Grouped paths can retain a reader leaf (e.g. distribution.number_density_m3).
# Bare names, renamed leaves, and ungrouped species paths are still violations.
READER_PATTERN = re.compile(r"(?<![\w.])(?:" + _alternatives(TABLE_NAMES | SPECIES_NAMES | REMOVED_MESH_NAMES) + r")(?!\w)")
SPECIES_PATTERN = re.compile(
    r"(?<![\w.])particles\.species(?:\[[^]\n]+\])?\.(?P<key>"
    + _alternatives(SPECIES_NAMES) + r")(?!\w)"
)
LEAF_PATTERN = re.compile(r"(?<!\w)(?:" + _alternatives(RENAMED_LEAVES) + r")(?!\w)")


def _diagnostic_reader_names(source: str) -> list[str]:
    source = TOKENS.sub(lambda match: "" if match[0].startswith("!") else match[0], source)
    source = re.sub(r"&\s*\n\s*&?", "", source)
    found = []
    # Semicolons outside literals separate statements too.
    masked = TOKENS.sub(lambda match: " " * len(match[0]), source)
    separators = [match.start() for match in re.finditer(r"[;\n]", masked)]
    start = 0
    statements = []
    for end in [*separators, len(source)]:
        statements.append(source[start:end])
        start = end + 1
    for statement in statements:
        # SELECT CASE values are the reader's dispatch keys, not diagnostics.
        if re.match(r"\s*case\s*\(", statement, re.IGNORECASE):
            continue
        for match in TOKENS.finditer(statement):
            # Literal second arguments of get_value are TOML lookup keys.
            if re.search(r"\bcall\s+get_value\s*\(\s*\w+\s*,\s*$",
                         statement[:match.start()], re.IGNORECASE):
                continue
            value = match[0][1:-1].replace("''", "'").replace('""', '"')
            names = {hit[0] for hit in READER_PATTERN.finditer(value)}
            names.update(hit["key"] for hit in SPECIES_PATTERN.finditer(value))
            names.update(hit[0] for hit in LEAF_PATTERN.finditer(value))
            if value in READER_TABLES:
                names.add(value)
            found.extend(sorted(names))
    return found


def test_fortran_diagnostics_use_public_config_paths() -> None:
    violations = []
    for path in sorted((ROOT / "src").rglob("*.f90")):
        if path.name in LAYOUT_FILES:
            continue
        names = _diagnostic_reader_names(path.read_text(encoding="utf-8"))
        if names:
            violations.append(f"{path.relative_to(ROOT)}: {', '.join(names)}")
    assert not violations, "Reader names in Fortran literals:\n" + "\n".join(violations)


@pytest.mark.parametrize("source, expected", [
    ("error stop 'sim.batch_duration must be finite'", ["sim.batch_duration"]),
    ("error stop 'mesh.templates.surface_side is invalid'", ["mesh.templates"]),
    ("error stop 'mesh.parts[2].path is required'", []),
    ('error stop "particles.species.q_particle must be nonzero"', ["q_particle"]),
    ("call require_toml_success(stat, 'sim.batch_duration')", ["sim.batch_duration"]),
    ("error stop 'bad ''q_particle'' ! value' ! sim.dt", ["q_particle"]),
    ("error stop &\n 'particle_boundary.z_high is invalid'", ["particle_boundary.z_high"]),
    ("case ('q_particle', 'm_particle')", []),
    ("case ('q_particle'); error stop 'q_particle'", ["q_particle"]),
    ("error stop 'particle_boundary.'//face", ["particle_boundary."]),
    ("error stop path//'.q_particle must be nonzero'", ["q_particle"]),
    ("error stop 'particles.species[2].number_density_m3 is invalid'", ["number_density_m3"]),
    ("error stop path//'.distribution.number_density_m3 must be positive'", []),
    ("call get_value(table, 'q_particle', value)", []),
    ("call get_value(table, 'q_particle', value); error stop 'q_particle'", ["q_particle"]),
    ("! error stop 'sim.dt'", []),
    ("error stop 'particles.species[2].charge_c must be nonzero'", []),
    ("error stop 'particles.reservoir.inflow_model is invalid'", []),
    ("write (u, '(a)') 'surface_current_model_phi0_V='", []),
])
def test_diagnostic_scan_distinguishes_reader_keys(source: str, expected: list[str]) -> None:
    assert _diagnostic_reader_names(source) == expected
