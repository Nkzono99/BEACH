from __future__ import annotations

import pytest

from beach.cli_estimate_fortran_workload import estimate_workload


def test_estimate_workload_resolves_batch_duration_from_step_and_species_targets() -> None:
    config = {
        "sim": {
            "batch_count": 3,
            "dt": 1.0,
            "batch_duration_step": 3.0,
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
                    "target_macro_particles_per_batch": 300,
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
                    "target_macro_particles_per_batch": 150,
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


def test_estimate_workload_supports_species_target_minus_one_following_species1_w() -> None:
    config = {
        "sim": {
            "batch_count": 3,
            "dt": 1.0,
            "batch_duration_step": 3.0,
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
                    "target_macro_particles_per_batch": 300,
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
                    "target_macro_particles_per_batch": -1,
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
    assert result["species_per_batch"] == [[300, 75], [300, 75], [300, 75]]


def test_estimate_workload_rejects_batch_duration_and_step_together() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "batch_duration": 1.0e-6,
            "batch_duration_step": 10.0,
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


def test_estimate_workload_rejects_removed_target_npcls_species1() -> None:
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
                    "npcls_per_step": 1,
                },
            ]
        },
    }

    with pytest.raises(SystemExit, match="was removed"):
        estimate_workload(config=config, threads=1)


def test_estimate_workload_rejects_w_and_target_together_for_reservoir() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "batch_duration": 1.0,
            "use_box": True,
        },
        "particles": {
            "species": [
                {
                    "source_mode": "reservoir_face",
                    "number_density_m3": 100.0,
                    "temperature_k": 0.0,
                    "m_particle": 1.0,
                    "w_particle": 10.0,
                    "target_macro_particles_per_batch": 10,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                    "drift_velocity": [0.0, 0.0, -1.0],
                }
            ]
        },
    }

    with pytest.raises(SystemExit, match="does not allow both"):
        estimate_workload(config=config, threads=1)


def test_estimate_workload_requires_w_or_target_for_reservoir() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "batch_duration": 1.0,
            "use_box": True,
        },
        "particles": {
            "species": [
                {
                    "source_mode": "reservoir_face",
                    "number_density_m3": 100.0,
                    "temperature_k": 0.0,
                    "m_particle": 1.0,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                    "drift_velocity": [0.0, 0.0, -1.0],
                }
            ]
        },
    }

    with pytest.raises(SystemExit, match="requires either w_particle or target"):
        estimate_workload(config=config, threads=1)


def test_estimate_workload_rejects_target_for_volume_seed() -> None:
    config = {
        "sim": {
            "batch_count": 1,
        },
        "particles": {
            "species": [
                {
                    "source_mode": "volume_seed",
                    "npcls_per_step": 5,
                    "target_macro_particles_per_batch": 10,
                }
            ]
        },
    }

    with pytest.raises(SystemExit, match="only valid for reservoir_face"):
        estimate_workload(config=config, threads=1)


def test_estimate_workload_rejects_minus_one_for_species1() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "batch_duration": 1.0,
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
                    "target_macro_particles_per_batch": -1,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                    "drift_velocity": [0.0, 0.0, -1.0],
                }
            ]
        },
    }

    with pytest.raises(SystemExit, match="cannot be -1"):
        estimate_workload(config=config, threads=1)


def test_estimate_workload_rejects_minus_one_if_species1_is_not_reservoir() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "batch_duration": 1.0,
            "use_box": True,
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [1.0, 1.0, 1.0],
        },
        "particles": {
            "species": [
                {
                    "source_mode": "volume_seed",
                    "npcls_per_step": 10,
                },
                {
                    "source_mode": "reservoir_face",
                    "number_density_m3": 100.0,
                    "temperature_k": 0.0,
                    "m_particle": 1.0,
                    "target_macro_particles_per_batch": -1,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                    "drift_velocity": [0.0, 0.0, -1.0],
                },
            ]
        },
    }

    with pytest.raises(SystemExit, match='source_mode=\"reservoir_face\"'):
        estimate_workload(config=config, threads=1)


def test_estimate_workload_keeps_batch_duration_behavior_with_manual_w() -> None:
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
