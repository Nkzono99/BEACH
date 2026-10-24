from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from beach import FieldKernel, FieldKernelOptions

pytestmark = pytest.mark.field_kernel_cache


def test_field_kernel_periodic_cache_reuses_only_matching_identity(
    tmp_path: Path, field_kernel_library: Path,
) -> None:
    lib = field_kernel_library
    source_centers = np.array(
        [
            [0.30, 0.40, -0.20],
            [1.50, 0.50, 0.16],
            [0.60, 1.60, 0.24],
            [1.64, 1.44, -0.12],
        ],
        dtype=float,
    )
    source_triangles = np.stack(
        (
            source_centers + [-0.025, -0.025, 0.0],
            source_centers + [0.025, -0.025, 0.0],
            source_centers + [-0.025, 0.025, 0.0],
        ),
        axis=1,
    )
    source_q = np.array([1.0e-9, -2.0e-9, 1.5e-9, -0.5e-9], dtype=float)
    options = FieldKernelOptions(
        theta=0.45,
        leaf_max=2,
        order=2,
        periodic2=((0, 1), (2.0, 2.0), (0.0, 0.0), 1, "cached_kneq0", 0.0, 4),
        box_min=(0.0, 0.0, -1.0),
        box_max=(2.0, 2.0, 1.0),
        periodic_cache_dir=str(tmp_path / "cache-\N{LATIN SMALL LETTER E WITH ACUTE}"),
        periodic_generation_tolerance=2.5e-9,
    )

    with FieldKernel(
        source_triangles, source_q, options=options, library_path=lib
    ) as cold_kernel:
        cold = cold_kernel.diagnostics()
    with FieldKernel(
        source_triangles, source_q, options=options, library_path=lib
    ) as warm_kernel:
        warm = warm_kernel.diagnostics()
    changed_options = replace(options, periodic_generation_tolerance=5.0e-9)
    with FieldKernel(
        source_triangles, source_q, options=changed_options, library_path=lib
    ) as changed_kernel:
        changed = changed_kernel.diagnostics()

    assert cold.periodic_cache_hit is False
    assert cold.periodic_operator_build_count == 1
    assert cold.periodic_cache_path is not None
    assert cold.periodic_cache_path.parent == tmp_path / "cache-\N{LATIN SMALL LETTER E WITH ACUTE}"
    assert cold.periodic_cache_path.is_file()
    assert cold.periodic_cache_fingerprint
    assert warm.periodic_cache_hit is True
    assert warm.periodic_operator_build_count == 0
    assert warm.periodic_cache_fingerprint == cold.periodic_cache_fingerprint
    assert warm.periodic_cache_path == cold.periodic_cache_path
    assert changed.periodic_cache_hit is False
    assert changed.periodic_operator_build_count == 1
    assert changed.periodic_cache_fingerprint
    assert changed.periodic_cache_fingerprint != cold.periodic_cache_fingerprint
    assert changed.periodic_cache_path is not None
    assert changed.periodic_cache_path != cold.periodic_cache_path
    assert changed.periodic_cache_path.is_file()
