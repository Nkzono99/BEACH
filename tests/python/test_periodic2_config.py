from __future__ import annotations

import copy

import pytest

from beach.fortran_results.periodic import (
    Periodic2Config,
    coerce_periodic2,
    periodic2_from_sim,
)


def test_periodic2_config_is_typed_and_legacy_tuple_compatible() -> None:
    raw = {
        "axes": [0, 1],
        "lengths": [2.0, 3.0],
        "box_min": [-1.0, -2.0, -3.0],
        "image_layers": 2,
        "far_correction": "auto",
    }
    original = copy.deepcopy(raw)

    config = coerce_periodic2(raw)

    assert isinstance(config, Periodic2Config)
    assert config == ((0, 1), (2.0, 3.0), (-1.0, -2.0), 2, "none", 0.0, 4)
    assert config.axes == (0, 1)
    assert config.image_layers == 2
    assert raw == original
    with pytest.raises(AttributeError):
        config.axes = (1, 2)  # type: ignore[misc]


def test_periodic2_config_revalidates_legacy_tuple_and_cached_policy() -> None:
    cached = ((0, 1), (2.0, 3.0), (0.0, 0.0), 1, "cached_kneq0", 0.0, 4)

    with pytest.raises(ValueError, match="far_correction"):
        coerce_periodic2(cached)
    assert coerce_periodic2(cached, allow_cached_kneq0=True) == cached


def test_periodic2_from_sim_keeps_split_periodic_table_out_of_scope() -> None:
    sim = {
        "field_bc_mode": "periodic2",
        "box_min": [1.0, 2.0, 3.0],
        "box_max": [5.0, 8.0, 10.0],
        "bc_x_low": "periodic",
        "bc_x_high": "periodic",
        "bc_y_low": "open",
        "bc_y_high": "open",
        "bc_z_low": "periodic",
        "bc_z_high": "periodic",
    }

    mapping = periodic2_from_sim(sim)
    config = coerce_periodic2(mapping)

    assert config is not None
    assert config.axes == (0, 2)
    assert config.lengths == (4.0, 7.0)
    assert config.origins == (1.0, 3.0)
    assert config.far_correction == "none"
    assert "nonzero_mode_backend" not in mapping  # type: ignore[operator]


def test_periodic2_from_free_sim_is_none() -> None:
    assert periodic2_from_sim({"field_bc_mode": "free"}) is None
