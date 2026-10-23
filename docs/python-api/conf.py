"""Build the Python reference from the checkout without loading a native library."""

from __future__ import annotations

import os
import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("MPLBACKEND", "Agg")

project = "BEACH Python API"
with (ROOT / "pyproject.toml").open("rb") as stream:
    release = tomllib.load(stream)["project"]["version"]
language = "ja"
extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon", "sphinx.ext.viewcode"]
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_type_aliases = {
    "Axes": "matplotlib.axes.Axes",
    "Axes3D": "mpl_toolkits.mplot3d.axes3d.Axes3D",
    "Figure": "matplotlib.figure.Figure",
    "FieldKernel": "beach.native.kernel.FieldKernel",
    "NDArray": "numpy.typing.NDArray",
}
napoleon_google_docstring = False
napoleon_numpy_docstring = True
# Dataclass attributes are indexed by autodoc; keep their NumPy descriptions inline.
napoleon_use_ivar = True
html_theme = "alabaster"
html_title = f"BEACH {release} Python API"
html_theme_options = {"description": "BEACH の Python 公開 API"}
html_baseurl = "https://nkzono99.github.io/BEACH/python/"
