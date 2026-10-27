from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]


def _load_golden():
    spec = importlib.util.spec_from_file_location(
        "beach_golden_regression", ROOT / "tools" / "golden_regression.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


golden = _load_golden()

BASE_CHARGES = "elem,charge\n1,1.0\n2,-2.5e-12\n"


def _tree(root: Path, files: dict[str, str]) -> Path:
    for rel, text in files.items():
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    return root


def _compare(tmp_path: Path, base: dict[str, str], new: dict[str, str], rtol: float) -> int:
    return golden.compare_trees(
        _tree(tmp_path / "base", base), _tree(tmp_path / "new", new), rtol
    )


def test_identical_trees_pass_ignoring_volatile_lines(tmp_path: Path) -> None:
    base = {"g1/summary.txt": "batches=3\nelapsed_s=1.0\n", "g1/charges.csv": BASE_CHARGES}
    new = {"g1/summary.txt": "batches=3\nelapsed_s=9.0\n", "g1/charges.csv": BASE_CHARGES}

    assert _compare(tmp_path, base, new, rtol=0.0) == 0


@pytest.mark.parametrize("rtol", [0.0, 1.0e-3])
@pytest.mark.parametrize("token", ["NaN", "Infinity", "-Infinity", "nan", "inf"])
def test_finite_to_nonfinite_fails(tmp_path: Path, token: str, rtol: float) -> None:
    new_charges = BASE_CHARGES.replace("1,1.0", f"1,{token}")

    status = _compare(
        tmp_path, {"g1/charges.csv": BASE_CHARGES}, {"g1/charges.csv": new_charges}, rtol
    )

    assert status == 1


def test_same_nonfinite_token_matches_but_different_one_fails(tmp_path: Path) -> None:
    nan_charges = BASE_CHARGES.replace("1,1.0", "1,NaN")
    within_rtol = nan_charges.replace("-2.5e-12", "-2.5000001e-12")
    assert (
        _compare(tmp_path / "same", {"g1/c.csv": nan_charges}, {"g1/c.csv": within_rtol}, 1.0e-6)
        == 0
    )

    inf_charges = nan_charges.replace("NaN", "Infinity")
    assert (
        _compare(tmp_path / "other", {"g1/c.csv": nan_charges}, {"g1/c.csv": inf_charges}, 1.0e-6)
        == 1
    )


def test_file_only_on_new_side_fails(tmp_path: Path) -> None:
    base = {"g1/charges.csv": BASE_CHARGES}
    new = {"g1/charges.csv": BASE_CHARGES, "g1/extra.txt": "x\n"}

    assert _compare(tmp_path, base, new, rtol=0.0) == 1


def test_case_only_on_new_side_fails(tmp_path: Path) -> None:
    base = {"g1/charges.csv": BASE_CHARGES}
    new = {"g1/charges.csv": BASE_CHARGES, "g2/charges.csv": BASE_CHARGES}

    assert _compare(tmp_path, base, new, rtol=0.0) == 1


def test_missing_file_and_empty_trees_fail(tmp_path: Path) -> None:
    base = {"g1/charges.csv": BASE_CHARGES, "g1/summary.txt": "batches=3\n"}
    assert _compare(tmp_path / "missing", base, {"g1/charges.csv": BASE_CHARGES}, 0.0) == 1

    (tmp_path / "empty" / "base").mkdir(parents=True)
    (tmp_path / "empty" / "new").mkdir(parents=True)
    assert golden.compare_trees(tmp_path / "empty" / "base", tmp_path / "empty" / "new", 0.0) == 1


def test_zero_rtol_requires_identical_text(tmp_path: Path) -> None:
    new_charges = BASE_CHARGES.replace("1,1.0", "1,1.00")

    assert _compare(tmp_path / "zero", {"g1/c.csv": BASE_CHARGES}, {"g1/c.csv": new_charges}, 0.0) == 1
    assert (
        _compare(tmp_path / "numeric", {"g1/c.csv": BASE_CHARGES}, {"g1/c.csv": new_charges}, 1.0e-14)
        == 0
    )


@pytest.mark.parametrize(
    ("new_value", "expected"),
    [("1.0000000001", 0), ("1.001", 1)],
)
def test_positive_rtol_compares_numbers_against_the_threshold(
    tmp_path: Path, new_value: str, expected: int
) -> None:
    new_charges = BASE_CHARGES.replace("1,1.0", f"1,{new_value}")

    status = _compare(tmp_path, {"g1/c.csv": BASE_CHARGES}, {"g1/c.csv": new_charges}, 1.0e-6)

    assert status == expected


def test_positive_rtol_fails_on_text_difference(tmp_path: Path) -> None:
    base = {"g1/summary.txt": "status=done\n"}
    new = {"g1/summary.txt": "status=failed\n"}

    assert _compare(tmp_path, base, new, rtol=1.0e-6) == 1
