from __future__ import annotations

import pytest

from beach.cli_estimate_fortran_workload import (
    completed_batches_from_resume_config,
    estimate_workload,
    read_macro_residuals,
)


def _fractional_reservoir_config() -> dict[str, object]:
    return {
        "sim": {
            "batch_count": 4,
            "batch_duration": 1.0,
            "use_box": True,
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [1.0, 1.0, 1.0],
        },
        "particles": {
            "species": [
                {
                    "source_mode": "reservoir_face",
                    "number_density_m3": 1.0,
                    "temperature_k": 0.0,
                    "m_particle": 1.0,
                    "w_particle": 4.0,
                    "inject_face": "z_low",
                    "pos_low": [0.0, 0.0, 0.0],
                    "pos_high": [1.0, 1.0, 0.0],
                    "drift_velocity": [0.0, 0.0, 1.0],
                }
            ]
        },
    }


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


def test_estimate_workload_resume_counts_only_remaining_batches() -> None:
    config = {
        "sim": {
            "batch_count": 5,
            "use_box": True,
        },
        "particles": {
            "species": [
                {
                    "species_key": "electron",
                    "source_mode": "volume_seed",
                    "npcls_per_step": 10,
                },
            ]
        },
    }

    result = estimate_workload(config=config, threads=2, completed_batches=2)

    assert result["target_batch_count"] == 5
    assert result["completed_batches"] == 2
    assert result["batch_count"] == 3
    assert result["batch_totals"] == [10, 10, 10]
    assert result["total_particles"] == 30


def test_estimate_workload_resume_all_batches_completed_is_zero_work() -> None:
    config = {
        "sim": {
            "batch_count": 5,
            "use_box": True,
        },
        "particles": {
            "species": [
                {
                    "source_mode": "volume_seed",
                    "npcls_per_step": 10,
                },
            ]
        },
    }

    result = estimate_workload(config=config, threads=2, completed_batches=5)

    assert result["batch_count"] == 0
    assert result["batch_totals"] == []
    assert result["total_particles"] == 0


def test_estimate_workload_rejects_resume_target_before_checkpoint() -> None:
    config = {
        "sim": {
            "batch_count": 2,
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

    with pytest.raises(SystemExit, match="completed checkpoint batches"):
        estimate_workload(config=config, threads=1, completed_batches=3)


def test_completed_batches_from_resume_config_reads_summary(tmp_path) -> None:
    out_dir = tmp_path / "outputs" / "latest"
    out_dir.mkdir(parents=True)
    (out_dir / "summary.txt").write_text(
        "mesh_nelem=1\nbatches=4\nlast_rel_change=0.0\n",
        encoding="utf-8",
    )
    config = {"output": {"resume": True, "dir": str(out_dir)}}

    assert completed_batches_from_resume_config(config) == 4


def test_completed_batches_from_resume_config_reads_restart_from(tmp_path) -> None:
    out_dir = tmp_path / "outputs" / "new"
    restart_dir = tmp_path / "outputs" / "parent"
    out_dir.mkdir(parents=True)
    restart_dir.mkdir(parents=True)
    (out_dir / "summary.txt").write_text(
        "mesh_nelem=1\nbatches=1\nlast_rel_change=0.0\n",
        encoding="utf-8",
    )
    (restart_dir / "summary.txt").write_text(
        "mesh_nelem=1\nbatches=6\nlast_rel_change=0.0\n",
        encoding="utf-8",
    )
    config = {
        "output": {
            "resume": True,
            "dir": str(out_dir),
            "restart_from": str(restart_dir),
        }
    }

    assert completed_batches_from_resume_config(config) == 6


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


def test_estimate_workload_rejects_removed_reserved_sim_key() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "use_hybrid": True,
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

    with pytest.raises(SystemExit, match=r"Unknown key in \[sim\]: use_hybrid"):
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
            "field_periodic_far_correction": "cached_kneq0",
            "field_periodic_ewald_layers": 6,
            "field_periodic_cache_dir": ".cache-test",
            "field_periodic_generation_tolerance": 1.0e-7,
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


def test_estimate_workload_accepts_normalization_and_sheath_sim_keys() -> None:
    config = {
        "sim": {
            "batch_count": 1,
            "field_normalization": "length",
            "field_length_scale": 2.0,
            "sheath_injection_model": "none",
            "sheath_alpha_deg": 60.0,
            "sheath_photoelectron_ref_density_cm3": 64.0,
            "sheath_reference_coordinate": 0.02,
            "sheath_electron_drift_mode": "normal",
            "sheath_ion_drift_mode": "normal",
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
                    "photo_escape_model": "boltzmann_cutoff",
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


def test_estimate_workload_splits_global_reservoir_count_by_mpi_rank() -> None:
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

    assert result["species_per_batch"] == [[3], [3], [3]]
    assert result["global_species_per_batch"] == [[10], [10], [10]]
    assert result["batch_totals"] == [3, 3, 3]
    assert result["global_batch_totals"] == [10, 10, 10]
    assert result["local_reservoir_particles"] == 9
    assert result["global_reservoir_particles"] == 30


@pytest.mark.parametrize("mpi_ranks", [1, 2, 4])
def test_estimate_workload_global_fractional_sequence_is_mpi_size_independent(
    mpi_ranks: int,
) -> None:
    results = [
        estimate_workload(
            config=_fractional_reservoir_config(),
            threads=1,
            mpi_ranks=mpi_ranks,
            mpi_rank=rank,
        )
        for rank in range(mpi_ranks)
    ]

    expected_global = [[0], [0], [0], [1]]
    assert all(result["global_species_per_batch"] == expected_global for result in results)
    assert [
        sum(result["species_per_batch"][batch_idx][0] for result in results)
        for batch_idx in range(4)
    ] == [0, 0, 0, 1]
    assert all(result["global_reservoir_particles"] == 1 for result in results)
    assert sum(result["local_reservoir_particles"] for result in results) == 1


def test_estimate_workload_resumed_global_reservoir_sequence_matches_uninterrupted() -> None:
    config = _fractional_reservoir_config()
    uninterrupted = estimate_workload(
        config=config, threads=1, mpi_ranks=4, mpi_rank=0
    )
    resumed = estimate_workload(
        config=config,
        threads=1,
        initial_residuals=[0.5],
        mpi_ranks=4,
        mpi_rank=0,
        completed_batches=2,
    )

    assert resumed["species_per_batch"] == uninterrupted["species_per_batch"][2:]
    assert resumed["global_species_per_batch"] == uninterrupted[
        "global_species_per_batch"
    ][2:]
    assert resumed["final_residuals"] == uninterrupted["final_residuals"]


def test_read_macro_residuals_rejects_legacy_rank_path(tmp_path) -> None:
    path = tmp_path / "macro_residuals_rank00000.csv"
    path.write_text("species_idx,residual\n1,0.5\n", encoding="utf-8")

    with pytest.raises(SystemExit, match="legacy rank-local"):
        read_macro_residuals(path, 1)


def test_estimate_workload_supports_velocity_grid_particle_flux() -> None:
    config = {
        "sim": {
            "batch_count": 2,
            "batch_duration": 1.0,
            "use_box": True,
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [1.0, 1.0, 1.0],
            "e0_abs": 5.0,
            "e0_phi_z_deg": 90.0,
        },
        "particles": {
            "species": [
                {
                    "source_mode": "reservoir_face",
                    "velocity_distribution": "grid",
                    "velocity_grid_path": "vgrid.csv",
                    "particle_flux_m2_s": 10.0,
                    "q_particle": -1.0,
                    "m_particle": 1.0,
                    "target_macro_particles_per_batch": 5,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                }
            ]
        },
    }

    result = estimate_workload(config=config, threads=2)

    assert result["species_per_batch"] == [[5], [5]]
    assert result["batch_totals"] == [5, 5]


def test_estimate_workload_supports_velocity_grid_current_density() -> None:
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
                    "velocity_distribution": "grid",
                    "velocity_grid_path": "vgrid.csv",
                    "current_density_a_m2": 8.0,
                    "q_particle": -2.0,
                    "m_particle": 1.0,
                    "target_macro_particles_per_batch": 4,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                }
            ]
        },
    }

    result = estimate_workload(config=config, threads=2)

    assert result["species_per_batch"] == [[4]]


def test_estimate_workload_rejects_velocity_grid_flux_and_current_together() -> None:
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
                    "velocity_distribution": "grid",
                    "velocity_grid_path": "vgrid.csv",
                    "particle_flux_m2_s": 10.0,
                    "current_density_a_m2": 8.0,
                    "q_particle": -2.0,
                    "m_particle": 1.0,
                    "target_macro_particles_per_batch": 4,
                    "inject_face": "z_high",
                    "pos_low": [0.0, 0.0, 1.0],
                    "pos_high": [1.0, 1.0, 1.0],
                }
            ]
        },
    }

    with pytest.raises(SystemExit, match="either particle_flux_m2_s or current_density_a_m2"):
        estimate_workload(config=config, threads=1)


def test_estimate_workload_rejects_mixed_external_e_field_forms() -> None:
    config = {
        "sim": {"batch_count": 1, "e0": [1.0, 0.0, 0.0], "e0_abs": 1.0},
        "particles": {
            "species": [
                {"source_mode": "volume_seed", "npcls_per_step": 1},
            ]
        },
    }

    with pytest.raises(SystemExit, match="sim.e0 cannot be combined"):
        estimate_workload(config=config, threads=1)


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
