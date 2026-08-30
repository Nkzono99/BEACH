from __future__ import annotations

import pytest

from beach.fortran_results.periodic import (
    Periodic2Config,
    coerce_periodic2,
    periodic2_from_sim,
)


def test_periodic2_config_normalizes_documented_mapping_contract() -> None:
    raw = {
        "axes": [0, 1],
        "lengths": [2.0, 3.0],
        "box_min": [-1.0, -2.0, -3.0],
        "image_layers": 2,
        "far_correction": "auto",
    }

    config = coerce_periodic2(raw)
    explicit_origins = coerce_periodic2({**raw, "origins": [4.0, 5.0]})

    assert isinstance(config, Periodic2Config)
    assert config == ((0, 1), (2.0, 3.0), (-1.0, -2.0), 2, "none", 0.0, 4)
    assert explicit_origins is not None
    assert explicit_origins.origins == (4.0, 5.0)


def test_periodic2_config_revalidates_legacy_tuple_and_cached_policy() -> None:
    cached = ((0, 1), (2.0, 3.0), (0.0, 0.0), 1, "cached_kneq0", 0.0, 4)

    with pytest.raises(ValueError, match="far_correction"):
        coerce_periodic2(cached)
    accepted = coerce_periodic2(cached, allow_cached_kneq0=True)
    assert isinstance(accepted, Periodic2Config)
    assert accepted.far_correction == "cached_kneq0"


def test_periodic2_from_sim_translates_only_periodic_geometry() -> None:
    sim = {
        "field_bc_mode": "periodic2",
        "box_min": [1.0, 2.0, 3.0],
        "box_max": [5.0, 8.0, 10.0],
        "bc_x_low": "periodic",
        "bc_x_high": "periodic",
        "bc_y_low": "periodic",
        "bc_y_high": "periodic",
        "bc_z_low": "open",
        "bc_z_high": "open",
    }

    mapping = periodic2_from_sim(sim)

    assert mapping == {
        "axes": (0, 1),
        "lengths": (4.0, 6.0),
        "origins": (1.0, 2.0),
        "image_layers": 1,
        "far_correction": "none",
        "ewald_alpha": 0.0,
        "ewald_layers": 4,
    }
    assert periodic2_from_sim({"field_bc_mode": "free"}) is None


def test_periodic2_from_sim_rejects_removed_root_oracle() -> None:
    sim = {
        "field_bc_mode": "periodic2",
        "box_min": [0.0, 0.0, -1.0],
        "box_max": [1.0, 1.0, 1.0],
        "bc_x_low": "periodic",
        "bc_x_high": "periodic",
        "bc_y_low": "periodic",
        "bc_y_high": "periodic",
        "bc_z_low": "open",
        "bc_z_high": "open",
        "field_periodic_far_correction": "m2l_root_oracle",
    }

    with pytest.raises(ValueError, match='was removed; use "none"'):
        periodic2_from_sim(sim)
