"""Keep the independent failure diagnostic honest on captured physical endpoints."""

import json
from pathlib import Path
import subprocess
import sys

import pytest


pytest.importorskip("scipy")

REPO_ROOT = Path(__file__).resolve().parents[2]

# The captured runs' physical inputs; geometry is irrelevant to the outer model.
CONFIG = """
[surface_current_model]
response_backend = "zhao_online"
photoelectron_closure = "energy_spectrum"
zhao_branch = "auto"
electron_species = "electron"
ion_species = "ion"
photoelectron_species = "photoelectron"

[[particles.species]]
species_key = "electron"
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
temperature_ev = 10.0
drift_velocity = [0.0, 0.0, -4.0e5]

[[particles.species]]
species_key = "ion"
q_particle = 1.602176634e-19
m_particle = 1.67262192369e-27
number_density_cm3 = 5.0
drift_velocity = [0.0, 0.0, -4.0e5]

[[particles.species]]
species_key = "photoelectron"
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
"""


@pytest.mark.parametrize(
    ("spectrum_name", "displacement_before", "expected_phi"),
    [
        pytest.param(
            "matching_plane_pe_captured_spectrum.csv",
            0.0,
            [5.829194512642798, 5.945909126703211, 6.055444426495559],
            id="three-physical-endpoints",
        ),
        pytest.param(
            "matching_plane_pe_batch14_spectrum.csv",
            1.9174651506096518e-11,
            [6.432956538216116],
            id="batch14-bisection-tolerance-miss",
        ),
        pytest.param(None, 0.0, [], id="zero-source-absence"),
    ],
)
def test_failure_diagnostic_classifies_physical_endpoints(
    tmp_path: Path,
    spectrum_name: str | None,
    displacement_before: float,
    expected_phi: list[float],
) -> None:
    config = tmp_path / "beach.toml"
    config.write_text(CONFIG)
    if spectrum_name is None:
        spectrum = tmp_path / "zero-source.csv"
        spectrum.write_text("energy_low_ev,energy_high_ev,flux_m2_s\n0.0,1.0,0.0\n")
    else:
        spectrum = REPO_ROOT / "tests" / "fixtures" / spectrum_name
    output = tmp_path / "diagnosis"
    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "tools" / "diagnose_matching_plane_failure.py"),
            str(spectrum),
            "--config", str(config),
            "--duration", "2",
            "--displacement-before", str(displacement_before),
            "--output", str(output),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    diagnosis = json.loads(output.with_suffix(".json").read_text())
    physical_roots = [
        candidate for candidate in diagnosis["candidates"]
        if candidate["kind"] == "backward_euler" and candidate["physical"]
    ]
    assert diagnosis["physical_BE_candidates"] == len(expected_phi)
    assert sorted(root["phi_H"] for root in physical_roots) == pytest.approx(
        expected_phi, rel=0.0, abs=1e-6,
    )
    assert diagnosis["absence_certificate"]["complete"] is (not expected_phi)
    if expected_phi:
        assert diagnosis["diagnosis"] == "physical_BE_endpoint_found"
    else:
        assert diagnosis["diagnosis"] == "type_B_absence_certified_in_model_with_roundoff_margin"
        assert diagnosis["absence_certificate"]["flat_state_excluded"] is True
