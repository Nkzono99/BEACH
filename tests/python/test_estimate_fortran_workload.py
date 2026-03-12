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


def test_estimate_workload_rejects_unknown_sim_key() -> None:
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

    with pytest.raises(SystemExit, match=r"Unknown key in \[sim\]"):
        estimate_workload(config=config, threads=1)


def test_estimate_workload_accepts_treecode_sim_keys() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "field_solver": "treecode",
            "tree_theta": 0.5,
            "tree_leaf_max": 16,
            "tree_min_nelem": 256,
            "use_box": True,
        },
        "particles": {
            "species": [
                {
                    "source_mode": "volume_seed",
                    "npcls_per_step": 3,
                },
            ]
        },
    }

    result = estimate_workload(config=config, threads=1)
    assert result["batch_totals"] == [3]


def test_estimate_workload_accepts_periodic_field_sim_keys() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "field_solver": "fmm",
            "field_bc_mode": "periodic2",
            "field_periodic_image_layers": 2,
            "field_periodic_far_correction": "ewald_like",
            "field_periodic_ewald_alpha": 1.2,
            "field_periodic_ewald_layers": 6,
            "use_box": True,
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [1.0, 1.0, 1.0],
            "bc_x_low": "periodic",
            "bc_x_high": "periodic",
            "bc_y_low": "periodic",
            "bc_y_high": "periodic",
            "bc_z_low": "open",
            "bc_z_high": "open",
        },
        "particles": {
            "species": [
                {
                    "source_mode": "volume_seed",
                    "npcls_per_step": 3,
                },
            ]
        },
    }

    result = estimate_workload(config=config, threads=1)
    assert result["batch_totals"] == [3]


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


def test_estimate_workload_supports_photo_raycast_as_upper_bound() -> None:
    config = {
        "sim": {
            "batch_count": 2,
            "batch_duration": 1.0e-6,
            "use_box": True,
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [1.0, 1.0, 1.0],
        },
        "particles": {
            "species": [
                {
                    "source_mode": "photo_raycast",
                    "emit_current_density_a_m2": 1.0e-3,
                    "rays_per_batch": 25,
                    "q_particle": -1.0,
                    "m_particle": 1.0,
                    "temperature_k": 0.0,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                }
            ]
        },
    }

    result = estimate_workload(config=config, threads=4)

    assert result["resolved_batch_duration"] == pytest.approx(1.0e-6)
    assert result["species_per_batch"] == [[25], [25]]


def test_estimate_workload_rejects_photo_raycast_with_outward_direction() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "batch_duration": 1.0e-6,
            "use_box": True,
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [1.0, 1.0, 1.0],
        },
        "particles": {
            "species": [
                {
                    "source_mode": "photo_raycast",
                    "emit_current_density_a_m2": 1.0e-3,
                    "rays_per_batch": 10,
                    "q_particle": -1.0,
                    "m_particle": 1.0,
                    "temperature_k": 0.0,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                    "ray_direction": [0.0, 0.0, 1.0],
                }
            ]
        },
    }

    with pytest.raises(SystemExit, match="must point inward"):
        estimate_workload(config=config, threads=1)


def test_estimate_workload_requires_positive_batch_duration_for_photo_raycast() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "use_box": True,
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [1.0, 1.0, 1.0],
        },
        "particles": {
            "species": [
                {
                    "source_mode": "photo_raycast",
                    "emit_current_density_a_m2": 1.0e-3,
                    "rays_per_batch": 10,
                    "q_particle": -1.0,
                    "m_particle": 1.0,
                    "temperature_k": 0.0,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                }
            ]
        },
    }

    with pytest.raises(SystemExit, match="batch_duration must be > 0"):
        estimate_workload(config=config, threads=1)


def test_estimate_workload_splits_volume_seed_by_mpi_rank() -> None:
    config = {
        "sim": {"batch_count": 2},
        "particles": {
            "species": [
                {"source_mode": "volume_seed", "npcls_per_step": 10},
            ]
        },
    }

    rank0 = estimate_workload(config=config, threads=2, mpi_ranks=3, mpi_rank=0)
    rank2 = estimate_workload(config=config, threads=2, mpi_ranks=3, mpi_rank=2)

    assert rank0["species_per_batch"] == [[4], [4]]
    assert rank0["batch_totals"] == [4, 4]
    assert rank0["batch_thread_min"] == [2, 2]
    assert rank0["batch_thread_max"] == [2, 2]

    assert rank2["species_per_batch"] == [[3], [3]]
    assert rank2["batch_totals"] == [3, 3]
    assert rank2["batch_thread_min"] == [1, 1]
    assert rank2["batch_thread_max"] == [2, 2]


def test_estimate_workload_splits_photo_raycast_by_mpi_rank() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "batch_duration": 1.0e-6,
            "use_box": True,
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [1.0, 1.0, 1.0],
        },
        "particles": {
            "species": [
                {
                    "source_mode": "photo_raycast",
                    "emit_current_density_a_m2": 1.0e-3,
                    "rays_per_batch": 10,
                    "q_particle": -1.0,
                    "m_particle": 1.0,
                    "temperature_k": 0.0,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                }
            ]
        },
    }

    rank0 = estimate_workload(config=config, threads=1, mpi_ranks=4, mpi_rank=0)
    rank3 = estimate_workload(config=config, threads=1, mpi_ranks=4, mpi_rank=3)

    assert rank0["species_per_batch"] == [[3]]
    assert rank3["species_per_batch"] == [[2]]


def test_estimate_workload_scales_reservoir_face_by_mpi_ranks() -> None:
    config = {
        "sim": {
            "batch_count": 3,
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
                    "w_particle": 10.0,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                    "drift_velocity": [0.0, 0.0, -1.0],
                }
            ]
        },
    }

    result = estimate_workload(config=config, threads=2, mpi_ranks=4, mpi_rank=0)

    assert result["species_per_batch"] == [[2], [3], [2]]
    assert result["batch_totals"] == [2, 3, 2]


def test_estimate_workload_rejects_invalid_mpi_rank_parameters() -> None:
    config = {
        "sim": {"batch_count": 1},
        "particles": {
            "species": [
                {"source_mode": "volume_seed", "npcls_per_step": 1},
            ]
        },
    }

    with pytest.raises(SystemExit, match="mpi_ranks must be > 0"):
        estimate_workload(config=config, threads=1, mpi_ranks=0, mpi_rank=0)
    with pytest.raises(SystemExit, match="mpi_rank must satisfy"):
        estimate_workload(config=config, threads=1, mpi_ranks=2, mpi_rank=2)
