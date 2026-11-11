from __future__ import annotations

from pathlib import Path

import tools.generate_m2l_order4_kernel as generator


def test_generated_m2l_order4_kernel_is_current() -> None:
    path = Path(__file__).resolve().parents[3] / generator.OUTPUT
    assert path.read_text(encoding="utf-8") == generator.generate(), (
        "run tools/generate_m2l_order4_kernel.py"
    )
