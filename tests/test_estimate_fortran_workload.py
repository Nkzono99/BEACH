from __future__ import annotations

import pytest

from beach.cli_estimate_fortran_workload import estimate_workload


def test_estimate_workload_resolves_batch_duration_from_target_species1() -> None:
    config = {
        "sim": {
            "batch_count": 3,
            "target_npcls_species1": 300,
            "use_box": True,
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [1.0, 1.0, 1.0],
        },
        "particles": {
            "species": [
                {
                    "source_mode": "reservoir_face",
                    "number_density_m3": 1000.0,
                    "temperature_k": 0.0,
                    "m_particle": 1.0,
                    "w_particle": 10.0,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                    "drift_velocity": [0.0, 0.0, -1.0],
                },
                {
                    "source_mode": "reservoir_face",
                    "number_density_m3": 250.0,
                    "temperature_k": 0.0,
                    "m_particle": 1.0,
                    "w_particle": 5.0,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                    "drift_velocity": [0.0, 0.0, -1.0],
                },
            ]
        },
    }

    result = estimate_workload(config=config, threads=8)

    assert result["resolved_batch_duration"] == pytest.approx(3.0)
    assert result["species_per_batch"] == [[300, 150], [300, 150], [300, 150]]


def test_estimate_workload_rejects_batch_duration_and_target_together() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "batch_duration": 1.0e-6,
            "target_npcls_species1": 10,
            "use_box": True,
        },
        "particles": {
            "species": [
                {
                    "source_mode": "reservoir_face",
                    "number_density_m3": 1.0,
                    "temperature_k": 0.0,
                    "m_particle": 1.0,
                    "w_particle": 1.0,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                    "drift_velocity": [0.0, 0.0, -1.0],
                }
            ]
        },
    }

    with pytest.raises(SystemExit, match="cannot be used together"):
        estimate_workload(config=config, threads=1)


def test_estimate_workload_requires_species1_reservoir_for_target() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "target_npcls_species1": 10,
            "use_box": True,
        },
        "particles": {
            "species": [
                {
                    "source_mode": "volume_seed",
                    "npcls_per_step": 5,
                },
                {
                    "source_mode": "reservoir_face",
                    "number_density_m3": 1.0,
                    "temperature_k": 0.0,
                    "m_particle": 1.0,
                    "w_particle": 1.0,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                    "drift_velocity": [0.0, 0.0, -1.0],
                },
            ]
        },
    }

    with pytest.raises(SystemExit, match="species\\[1\\].source_mode"):
        estimate_workload(config=config, threads=1)


def test_estimate_workload_keeps_legacy_batch_duration_behavior() -> None:
    config = {
        "sim": {
            "batch_count": 3,
            "batch_duration": 0.5,
            "use_box": True,
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [1.0, 1.0, 1.0],
        },
        "particles": {
            "species": [
                {
                    "source_mode": "reservoir_face",
                    "number_density_m3": 100.0,
                    "temperature_k": 0.0,
                    "m_particle": 1.0,
                    "w_particle": 10.0,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                    "drift_velocity": [0.0, 0.0, -1.0],
                }
            ]
        },
    }

    result = estimate_workload(config=config, threads=2)

    assert result["resolved_batch_duration"] == pytest.approx(0.5)
    assert result["species_per_batch"] == [[5], [5], [5]]
