from __future__ import annotations

import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def _read(relative_path: str) -> str:
    return (REPO_ROOT / relative_path).read_text(encoding="utf-8")


def _normalized(relative_path: str) -> str:
    return " ".join(_read(relative_path).split())


def test_fortran_far_correction_defaults_and_normalization_match() -> None:
    sim_types = _normalized("src/core/bem_types.f90")
    app_defaults = _normalized("src/config/bem_app_config_types.f90")
    solver_config = _normalized(
        "src/physics/field_solver/bem_field_solver_config.f90"
    )
    fmm_types = _normalized(
        "src/physics/field_solver/fmm/internal/common/bem_coulomb_fmm_types.f90"
    )
    core_plan_ops = _normalized(
        "src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_plan_ops.f90"
    )

    assert "field_periodic_far_correction = 'none'" in sim_types
    assert "field_periodic_far_correction = 'none'" in app_defaults
    assert "periodic_far_correction = 'none'" in fmm_types
    assert "case ('auto') self%periodic_far_correction = 'none'" in solver_config
    assert (
        "case ('auto') plan%options%periodic_far_correction = 'none'"
        in core_plan_ops
    )


def test_far_correction_schema_copies_and_metadata_match() -> None:
    schema_paths = (
        "schemas/beach.schema.json",
        "beach/config/schemas/beach.schema.json",
        "plugins/beach-context/references/schemas/beach.schema.json",
    )
    schema_texts = [_read(path) for path in schema_paths]

    assert schema_texts[0] == schema_texts[1] == schema_texts[2]
    schema = json.loads(schema_texts[0])
    far = schema["$defs"]["sim"]["properties"][
        "field_periodic_far_correction"
    ]
    assert far["default"] == "none"
    assert far["enum"] == ["auto", "none", "m2l_root_oracle", "cached_kneq0"]
    assert "cached_kneq0" in far["description"]
    assert "high-cost diagnostic" in far["description"]


def test_spec_and_canonical_docs_state_the_staged_contract() -> None:
    spec = _read("SPEC.md")
    parameters_ja = _read("docs/Parameters.md")
    fmm_core_ja = _read("docs/FMMCore.md")
    fmm_core_en = _read("docs/FMMCore.en.md")

    assert '`field_periodic_far_correction="none"`' in spec
    assert "`auto`" in spec and "互換用" in spec and "`none` に正規化" in spec
    assert (
        "`m2l_root_oracle`" in spec
        and "明示" in spec
        and "高コスト診断" in spec
    )
    assert "有限画像" in spec
    assert "cached_kneq0" in spec and "generator version" in spec
    assert "production" in spec and "O(log n)" in spec

    assert "`periodic2` の遠方補正は `auto` を既定" not in spec
    assert "内部的に `m2l_root_oracle` に正規化" not in spec
    assert "`auto` の alias ではありません" not in spec

    assert (
        '`field_periodic_far_correction="auto"` は互換用に受理され'
        in parameters_ja
    )
    assert "現在は `none` と同じ扱い" in parameters_ja
    assert "明示 opt-in の高コスト診断モード" in fmm_core_ja
    assert "production" in fmm_core_ja and "`cached_kneq0`" in fmm_core_ja
    assert "Infinite-periodic production runs use `cached_kneq0`" in fmm_core_en


def test_plugin_references_match_canonical_files_without_stale_contract() -> None:
    mirrors = {
        "plugins/beach-context/references/SPEC.md": "SPEC.md",
        "plugins/beach-context/references/fortran_parameter_file.md": (
            "docs/Parameters.md"
        ),
        "plugins/beach-context/references/fortran_fmm_core.md": "docs/FMMCore.md",
        "plugins/beach-context/references/periodic_zero_mode_outer_plasma.md": (
            "docs/PeriodicZeroModeOuterPlasma.md"
        ),
        "plugins/beach-context/references/sheath_reservoir_boundary.md": (
            "docs/SheathReservoirBoundary.md"
        ),
        "plugins/beach-context/references/examples/beach.toml": (
            "examples/beach.toml"
        ),
        "plugins/beach-context/references/examples/periodic2_basic/beach.toml": (
            "examples/periodic2_basic/beach.toml"
        ),
        "plugins/beach-context/references/agent-user-guide.md": (
            "docs/agent-user-guide.md"
        ),
        "plugins/beach-context/references/python_postprocess_api.md": (
            "docs/PythonPostprocessAPI.md"
        ),
    }
    stale_sentinels = (
        "`periodic2` の遠方補正は `auto` を既定",
        "内部的に `m2l_root_oracle` に正規化",
        "`auto` の alias ではありません",
        '"auto"            — periodic2 + fmm では内部的に m2l_root_oracle へ正規化',
        '`periodic2` の `auto` は `m2l_root_oracle` に正規化されます',
    )

    for mirror_path, canonical_path in mirrors.items():
        mirror = _read(mirror_path)
        assert mirror == _read(canonical_path), (
            f"{mirror_path} must mirror {canonical_path} byte-for-byte"
        )
        for sentinel in stale_sentinels:
            assert sentinel not in mirror, (
                f"{mirror_path} contains stale far-correction text: {sentinel}"
            )


def test_field_kernel_real_cache_tests_are_opt_in_only() -> None:
    makefile = _read("Makefile")
    cache_test = _read("tests/python/test_field_kernel_cache.py")

    assert (
        "test_field_kernel_cache_c"
        in makefile.split("FORTRAN_FAR_CORRECTION_TARGETS ?=", 1)[1]
    )
    assert (
        "test_field_kernel_cache_c"
        not in makefile.split("FORTRAN_L2_TARGETS ?=", 1)[1].split(
            "FORTRAN_L3_TARGETS ?=", 1
        )[0]
    )
    assert 'os.environ.get("BEACH_RUN_FIELD_KERNEL_CACHE_TESTS") != "1"' in cache_test
    assert "allow_module_level=True" in cache_test
    assert makefile.count("BEACH_RUN_FIELD_KERNEL_CACHE_TESTS=1") == 1
