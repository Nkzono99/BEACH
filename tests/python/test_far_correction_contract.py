from __future__ import annotations

import json
import os
import runpy
import subprocess
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
GROUPED_BACKENDS = frozenset({"finite_images", "cached_kneq0", "panel_spectral_reference"})


def _read(relative_path: str) -> str:
    return (REPO_ROOT / relative_path).read_text(encoding="utf-8")


def _line_starting_with(text: str, prefix: str) -> str:
    for line in text.splitlines():
        if line.startswith(prefix):
            return line
    raise AssertionError(f"missing line starting with {prefix!r}")


def test_far_correction_schema_exposes_public_backends() -> None:
    schema = json.loads(_read("schemas/beach.schema.json"))
    backend = schema["$defs"]["groupedConfig"]["properties"]["fields"][
        "properties"
    ]["periodic"]["properties"]["backend"]

    assert backend["type"] == "string"
    assert set(backend["enum"]) == GROUPED_BACKENDS


def test_canonical_docs_state_far_correction_semantics() -> None:
    spec = _read("SPEC.md")
    spec_rows = {
        mode: _line_starting_with(spec, f"- `{mode}`:")
        for mode in GROUPED_BACKENDS
    }
    assert "有限" in spec_rows["finite_images"] and "image" in spec_rows["finite_images"]
    assert "direct" in spec_rows["panel_spectral_reference"]
    assert "無限周期" in spec_rows["cached_kneq0"]
    assert "Fourier" in spec_rows["cached_kneq0"]

    contracts = {
        "docs/developer/PeriodicFarCorrection.md": {
            "rows": {
                "finite_images": ("有限画像",),
                "cached_kneq0": ("無限周期", "演算子"),
                "panel_spectral_reference": ("Fourier", "検証"),
            },
        },
    }
    for path, contract in contracts.items():
        text = _read(path)
        assert set(contract["rows"]) == GROUPED_BACKENDS
        for mode, fragments in contract["rows"].items():
            row = _line_starting_with(text, f"| `{mode}` |")
            assert all(fragment in row for fragment in fragments), (path, mode)

        assert "物理的な平均場はこの演算子の外で一度だけ加えます" in text
        assert "periodic-kneq0-generator-v12" in text
        assert "BEACH_CACHE_DIR" in text


def test_bundled_plugin_references_match_canonical_sources() -> None:
    mirrors = {
        "plugins/beach-context/references/SPEC.md": "SPEC.md",
        "plugins/beach-context/references/fortran_parameter_file.md": (
            "docs/reference/Parameters.md"
        ),
        "plugins/beach-context/references/fortran_fmm_core.md": "docs/developer/FMMCore.md",
        "plugins/beach-context/references/examples/beach.toml": (
            "examples/beach.toml"
        ),
        "plugins/beach-context/references/examples/periodic2_basic.toml": (
            "examples/periodic2_basic.toml"
        ),
        "plugins/beach-context/references/agent-user-guide.md": (
            "docs/agent-user-guide.md"
        ),
        "plugins/beach-context/references/python_postprocess_api.md": (
            "docs/reference/PythonPostprocessAPI.md"
        ),
    }

    for source in REPO_ROOT.joinpath("examples").glob("*.toml"):
        mirrors[f"plugins/beach-context/references/examples/{source.name}"] = f"examples/{source.name}"

    for mirror_path, canonical_path in mirrors.items():
        mirror = (REPO_ROOT / mirror_path).read_bytes()
        canonical = (REPO_ROOT / canonical_path).read_bytes()
        assert mirror == canonical, (
            f"{mirror_path} must mirror {canonical_path} byte-for-byte"
        )


def test_field_kernel_cache_tests_are_opt_in(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    environment = os.environ.copy()
    for name in (
        "FORTRAN_FAR_CORRECTION_TARGETS",
        "FORTRAN_L1_TARGETS",
        "FORTRAN_L2_TARGETS",
        "FORTRAN_L3_TARGETS",
        "MAKEFLAGS",
        "MAKEOVERRIDES",
        "MFLAGS",
    ):
        environment.pop(name, None)
    result = subprocess.run(
        ["make", "-pn", "--no-print-directory", "-f", "Makefile"],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr

    far_targets = set(
        _line_starting_with(
            result.stdout, "FORTRAN_FAR_CORRECTION_TARGETS = "
        ).split(" = ", maxsplit=1)[1].split()
    )
    assert "test_field_kernel_cache_c" in far_targets
    for tier in ("L1", "L2", "L3"):
        tier_targets = set(
            _line_starting_with(
                result.stdout, f"FORTRAN_{tier}_TARGETS = "
            ).split(" = ", maxsplit=1)[1].split()
        )
        assert "test_field_kernel_cache_c" not in tier_targets

    opt_in_assignment = _line_starting_with(
        result.stdout,
        "test-field-kernel-cache: BEACH_RUN_FIELD_KERNEL_CACHE_TESTS ",
    )
    assert opt_in_assignment.rsplit("=", maxsplit=1)[1].strip() == "1"

    default_routes = subprocess.run(
        [
            "make",
            "-n",
            "--trace",
            "--no-print-directory",
            "-f",
            "Makefile",
            "test-l1",
            "test-l2",
            "test-l3",
            "test-physics-release",
        ],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert default_routes.returncode == 0, default_routes.stderr
    assert "target 'test-field-kernel-cache'" not in default_routes.stdout

    monkeypatch.delenv("BEACH_RUN_FIELD_KERNEL_CACHE_TESTS", raising=False)
    with pytest.raises(pytest.skip.Exception, match="opt-in"):
        runpy.run_path(str(REPO_ROOT / "tests/python/test_field_kernel_cache.py"))
