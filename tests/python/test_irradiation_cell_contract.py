"""Keep experiment assumptions in the irradiation-cell contract, outside kernels."""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = "src/runtime/sheath/bem_irradiation_cell_contract.f90"
# Exact statements, rather than whole-file exclusions: both are generic six-face mappings.
ALLOWED = {
    ("src/runtime/configuration/bem_app_config_particle_runtime.f90", "name = 'z_high'"):
        "boundary_face_name maps every bem_boundary face index to its name",
    ("src/runtime/configuration/bem_app_config_potential_runtime.f90", "case ('z_high')"):
        "resolve_face_sampling_geometry samples any of the six box faces",
}
PATTERNS = {
    "wall bit": r"\bbtest\s*\(.*?,\s*5(?:_i32)?\s*\)",
    "wall boundary": r"%face_bc\s*\(\s*6(?:_i32)?\s*\)",
    "kinetic wall face": r"%(?:outflow_barrier_face|inflow_kinetic_face)\b.*=\s*6(?:_i32)?\b",
    "wall barrier query": r"\bopen_face_uses_kinetic_barrier\s*\(\s*6(?:_i32)?\s*,",
    "wall coordinate": r"%box_max\s*\(\s*3(?:_i32)?\s*\)",
    "wall face name": r"['\"]z_high['\"]",
    "sheath normal drift": r"%drift_velocity\s*\(\s*3(?:_i32)?\s*\)",
    "fixed periodic tangents": r"%bc_(?:low|high)\s*\(\s*1\s*:\s*2\s*\).*\bbc_periodic\b",
    "lower model literal": r"['\"](?:e_bottom_zero|symmetric_vacuum|legacy_implicit)['\"]",
}


def _statements(source: str):
    """Join Fortran continuations, dropping comments without stripping quoted '!'."""
    current = ""
    start = 0
    for number, line in enumerate(source.splitlines(), 1):
        code = re.match(r"(?:'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\"|[^!'\"])*", line)[0]
        code = code.strip()
        if not code:
            continue
        if not current:
            start = number
        current += " " + code.lstrip("&").rstrip("&").strip()
        if code.endswith("&"):
            continue
        yield start, " ".join(current.split()).lower()
        current = ""
    assert not current, "unfinished Fortran statement"


def _violations(statement: str) -> list[str]:
    return [name for name, pattern in PATTERNS.items() if re.search(pattern, statement, re.I)]


def test_runtime_queries_irradiation_cell_contract():
    failures = []
    allowed_seen = set()
    for directory in ("src/runtime", "src/physics/sheath"):
        for path in sorted((ROOT / directory).rglob("*")):
            if path.suffix.lower() != ".f90":
                continue
            relative = path.relative_to(ROOT).as_posix()
            if relative == CONTRACT:
                continue
            for number, statement in _statements(path.read_text()):
                violations = _violations(statement)
                if not violations:
                    continue
                key = (relative, statement)
                if key in ALLOWED:
                    allowed_seen.add(key)
                else:
                    failures.append(f"{relative}:{number}: {', '.join(violations)}: {statement}")
    assert not failures, "\n".join(failures)
    assert allowed_seen == ALLOWED.keys(), "Remove stale assumption exceptions"
    assert all(ALLOWED.values()), "Each allowed statement needs a reason"


@pytest.mark.parametrize("source", [
    "if (btest(event%face_mask, &\n  5_i32)) return",
    "if (event%face_bc(6_i32) == bc_open) return",
    "result%outflow_barrier_face(photo_idx) = 6_i32",
    "result%inflow_kinetic_face([electron_idx, ion_idx]) = 6_i32",
    "uses = open_face_uses_kinetic_barrier(6_i32, boundary)",
    "reference = app%sim%box_max(3)",
    "call sample_face(mesh, 'z_high', potential)",
    "drift = -app%particle_species(idx)%drift_velocity(3)",
    "uniform = all(app%sim%bc_low(1:2) == bc_periodic)",
    "model = 'e_bottom_zero'",
])
def test_assumption_detector_rejects_wall_literals(source):
    assert any(_violations(statement) for _, statement in _statements(source))


def test_assumption_detector_keeps_dimensions_generic_faces_and_output_keys():
    source = """
    real(dp) :: box_max(3)
    do face = 1_i32, 6_i32
    coordinate = sim%box_max(axis) ! Wall used to be sim%box_max(3)
    write (u, '(a)') 'top_reference_z_high=,field_periodic_bc_high_z=z_high'
    """
    assert not any(_violations(statement) for _, statement in _statements(source))


def test_generic_kernels_do_not_depend_on_experiment_contract():
    directories = ("src/mesh", "src/particles", "src/physics")
    for directory in directories:
        for path in (ROOT / directory).rglob("*"):
            if path.suffix.lower() != ".f90" or "sheath" in path.parts:
                continue
            for number, statement in _statements(path.read_text()):
                assert not re.search(r"\buse\s+bem_irradiation_cell_contract\b", statement), f"{path}:{number}"
