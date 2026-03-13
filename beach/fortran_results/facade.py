"""Object-oriented facade for one Fortran output directory."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Literal, Mapping, TYPE_CHECKING

import numpy as np

from .animation import animate_history_mesh
from .coulomb import calc_coulomb
from .io import load_fortran_result
from .plotting import (
    plot_charge_mesh,
    plot_charges,
    plot_mesh_source_boxplot,
    plot_potential_mesh,
    plot_potential_slices,
)
from .potential import (
    compute_potential_mesh,
    compute_potential_points,
    compute_potential_slices,
)
from .selection import _build_mesh_selection, _mesh_ids_or_default
from .types import FortranRunResult, MeshSelection

if TYPE_CHECKING:
    from matplotlib.animation import FuncAnimation


class Beach:
    """High-level facade for one Fortran output directory.

    Parameters
    ----------
    output_dir : str or pathlib.Path, default "outputs/latest"
        Directory containing Fortran output files.
    """

    def __init__(
        self,
        output_dir: str | Path = "outputs/latest",
    ) -> None:
        self.output_dir = Path(output_dir)
        self._result: FortranRunResult | None = None

    @property
    def result(self) -> FortranRunResult:
        """Return the loaded run result.

        Returns
        -------
        FortranRunResult
            Parsed result for ``output_dir``.
        """

        if self._result is None:
            self._result = load_fortran_result(self.output_dir)
        return self._result

    def reload(self) -> FortranRunResult:
        """Reload output files from disk.

        Returns
        -------
        FortranRunResult
            Refreshed parsed result.
        """

        self._result = load_fortran_result(self.output_dir)
        return self._result

    @property
    def mesh_ids(self) -> tuple[int, ...]:
        """Return available mesh ids.

        Returns
        -------
        tuple of int
            Sorted unique mesh identifiers.
        """

        ids = np.unique(_mesh_ids_or_default(self.result))
        return tuple(int(v) for v in ids)

    def get_mesh(self, *mesh_ids: int, step: int | None = -1):
        """Return mesh selections for requested ids.

        Parameters
        ----------
        *mesh_ids : int
            One or more mesh identifiers.
        step : int or None, default -1
            History step used to select charges. ``-1`` selects latest history,
            ``None`` uses final charges.

        Returns
        -------
        MeshSelection or tuple[MeshSelection, ...]
            Single selection when one id is given, tuple otherwise.

        Raises
        ------
        ValueError
            If no mesh ids are provided or a mesh id is invalid.
        """

        if len(mesh_ids) == 0:
            raise ValueError("at least one mesh id must be provided.")
        selections = tuple(
            _build_mesh_selection(self.result, (int(mesh_id),), step=step)
            for mesh_id in mesh_ids
        )
        if len(selections) == 1:
            return selections[0]
        return selections

    def get_mesh_charge(self, *mesh_ids: int, step: int | None = -1):
        """Return per-element charges for requested mesh ids.

        Parameters
        ----------
        *mesh_ids : int
            One or more mesh identifiers.
        step : int or None, default -1
            History step used to select charges.

        Returns
        -------
        numpy.ndarray or tuple[numpy.ndarray, ...]
            Charge array for each requested mesh.
        """

        selection = self.get_mesh(*mesh_ids, step=step)
        if isinstance(selection, tuple):
            return tuple(mesh.charges.copy() for mesh in selection)
        return selection.charges.copy()

    def calc_coulomb(
        self,
        target: int | MeshSelection | Iterable[int | MeshSelection],
        source: int | MeshSelection | Iterable[int | MeshSelection],
        *,
        step: int | None = -1,
        softening: float = 0.0,
        torque_origin: Literal[
            "target_center",
            "source_center",
            "origin",
            "group_a_center",
            "group_b_center",
        ] = "target_center",
    ):
        """Compute Coulomb force/torque from source acting on target.

        Parameters
        ----------
        target : int, MeshSelection, or iterable of those
            Target mesh group (group A).
        source : int, MeshSelection, or iterable of those
            Source mesh group (group B).
        step : int or None, default -1
            History step used to select charges.
        softening : float, default 0.0
            Softening length in meters.
        torque_origin : {"target_center", "source_center", "origin", "group_a_center", "group_b_center"}, default "target_center"
            Torque reference point.

        Returns
        -------
        CoulombInteraction
            Aggregated force/torque summary.
        """

        return calc_coulomb(
            self.result,
            target,
            source,
            step=step,
            softening=softening,
            torque_origin=torque_origin,
        )

    def plot_mesh(self, *, cmap: str = "coolwarm"):
        """Plot a 3D mesh colored by surface charge density.

        Parameters
        ----------
        cmap : str, default "coolwarm"
            Matplotlib colormap name.

        Returns
        -------
        tuple
            ``(figure, axes)`` from matplotlib.
        """

        return plot_charge_mesh(self.result, cmap=cmap)

    def plot_bar(self):
        """Plot per-element charges as a bar chart.

        Returns
        -------
        tuple
            ``(figure, axes)`` from matplotlib.
        """

        return plot_charges(self.result)

    def compute_potential(
        self,
        *,
        softening: float | None = None,
        self_term: str = "auto",
        periodic2: Mapping[str, object] | None = None,
        reference_point: Iterable[float] | str | None = None,
    ) -> np.ndarray:
        """Compute potential values at triangle centroids.

        Parameters
        ----------
        softening : float or None, default None
            Softening length in meters.
        self_term : {"auto", "area_equivalent", "exclude", "softened_point"}, default "auto"
            Self-interaction model.
        periodic2 : mapping or None, default None
            Two-axis periodic setting. See
            :func:`beach.fortran_results.compute_potential_mesh`.
            ``None`` の場合は出力ディレクトリ近傍の ``beach.toml`` から自動判定。
        reference_point : iterable of float, {"species1_injection_center"}, or None, default None
            基準電位を差し引く参照点。

        Returns
        -------
        numpy.ndarray
            Potential values in volts with shape ``(mesh_nelem,)``.
        """

        return compute_potential_mesh(
            self.result,
            softening=softening,
            self_term=self_term,
            periodic2=periodic2,
            reference_point=reference_point,
        )

    def compute_potential_points(
        self,
        points: np.ndarray,
        *,
        softening: float | None = None,
        chunk_size: int = 2048,
        periodic2: Mapping[str, object] | None = None,
        reference_point: Iterable[float] | str | None = None,
    ) -> np.ndarray:
        """Compute potential values at arbitrary 3D points.

        Parameters
        ----------
        points : numpy.ndarray
            Sampling points with shape ``(n_points, 3)``.
        softening : float or None, default None
            Softening length in meters.
        chunk_size : int, default 2048
            Number of points processed per chunk.
        periodic2 : mapping or None, default None
            Two-axis periodic setting. See
            :func:`beach.fortran_results.compute_potential_points`.
            ``None`` の場合は出力ディレクトリ近傍の ``beach.toml`` から自動判定。
        reference_point : iterable of float, {"species1_injection_center"}, or None, default None
            基準電位を差し引く参照点。

        Returns
        -------
        numpy.ndarray
            Potential values in volts with shape ``(n_points,)``.
        """

        return compute_potential_points(
            self.result,
            points,
            softening=softening,
            chunk_size=chunk_size,
            periodic2=periodic2,
            reference_point=reference_point,
        )

    def compute_potential_slices(
        self,
        *,
        box_min: Iterable[float],
        box_max: Iterable[float],
        grid_n: int = 200,
        xy_z: float | None = None,
        yz_x: float | None = None,
        xz_y: float | None = None,
        softening: float | None = None,
        chunk_size: int = 2048,
        periodic2: Mapping[str, object] | None = None,
        reference_point: Iterable[float] | str | None = None,
    ):
        """Compute potential slices on XY/YZ/XZ planes.

        Parameters
        ----------
        box_min : iterable of float
            Lower simulation-box corner ``[x, y, z]``.
        box_max : iterable of float
            Upper simulation-box corner ``[x, y, z]``.
        grid_n : int, default 200
            Grid size per slice axis.
        xy_z : float or None, default None
            Z coordinate of XY slice.
        yz_x : float or None, default None
            X coordinate of YZ slice.
        xz_y : float or None, default None
            Y coordinate of XZ slice.
        softening : float or None, default None
            Softening length in meters.
        chunk_size : int, default 2048
            Sampling chunk size.
        periodic2 : mapping or None, default None
            Two-axis periodic setting. See
            :func:`beach.fortran_results.compute_potential_points`.
            ``None`` の場合は出力ディレクトリ近傍の ``beach.toml`` から自動判定。
        reference_point : iterable of float, {"species1_injection_center"}, or None, default None
            基準電位を差し引く参照点。

        Returns
        -------
        dict[str, PotentialSlice2D]
            Slice data keyed by ``"xy"``, ``"yz"``, and ``"xz"``.
        """

        return compute_potential_slices(
            self.result,
            box_min=box_min,
            box_max=box_max,
            grid_n=grid_n,
            xy_z=xy_z,
            yz_x=yz_x,
            xz_y=xz_y,
            softening=softening,
            chunk_size=chunk_size,
            periodic2=periodic2,
            reference_point=reference_point,
        )

    def plot_potential(
        self,
        *,
        softening: float | None = None,
        self_term: str = "auto",
        cmap: str = "viridis",
        periodic2: Mapping[str, object] | None = None,
        reference_point: Iterable[float] | str | None = "species1_injection_center",
    ):
        """Plot a 3D mesh colored by reconstructed electric potential.

        Parameters
        ----------
        softening : float or None, default None
            Softening length in meters.
        self_term : {"auto", "area_equivalent", "exclude", "softened_point"}, default "auto"
            Self-interaction model.
        cmap : str, default "viridis"
            Matplotlib colormap name.
        periodic2 : mapping or None, default None
            Two-axis periodic setting. See
            :func:`beach.fortran_results.compute_potential_mesh`.
            ``None`` の場合は出力ディレクトリ近傍の ``beach.toml`` から自動判定。
        reference_point : iterable of float, {"species1_injection_center"}, or None, default "species1_injection_center"
            基準電位を差し引く参照点。

        Returns
        -------
        tuple
            ``(figure, axes)`` from matplotlib.
        """

        return plot_potential_mesh(
            self.result,
            softening=softening,
            self_term=self_term,
            cmap=cmap,
            periodic2=periodic2,
            reference_point=reference_point,
        )

    def plot_potential_slices(
        self,
        *,
        box_min: Iterable[float],
        box_max: Iterable[float],
        grid_n: int = 200,
        xy_z: float | None = None,
        yz_x: float | None = None,
        xz_y: float | None = None,
        softening: float | None = None,
        chunk_size: int = 2048,
        cmap: str = "viridis",
        vmin: float | None = None,
        vmax: float | None = None,
        periodic2: Mapping[str, object] | None = None,
        reference_point: Iterable[float] | str | None = "species1_injection_center",
    ):
        """Plot XY/YZ/XZ potential slices with a shared color scale.

        Parameters
        ----------
        box_min : iterable of float
            Lower simulation-box corner ``[x, y, z]``.
        box_max : iterable of float
            Upper simulation-box corner ``[x, y, z]``.
        grid_n : int, default 200
            Grid size per slice axis.
        xy_z : float or None, default None
            Z coordinate of XY slice.
        yz_x : float or None, default None
            X coordinate of YZ slice.
        xz_y : float or None, default None
            Y coordinate of XZ slice.
        softening : float or None, default None
            Softening length in meters.
        chunk_size : int, default 2048
            Sampling chunk size.
        cmap : str, default "viridis"
            Matplotlib colormap name.
        vmin : float or None, default None
            Lower color limit.
        vmax : float or None, default None
            Upper color limit.
        periodic2 : mapping or None, default None
            Two-axis periodic setting. See
            :func:`beach.fortran_results.compute_potential_points`.
            ``None`` の場合は出力ディレクトリ近傍の ``beach.toml`` から自動判定。
        reference_point : iterable of float, {"species1_injection_center"}, or None, default "species1_injection_center"
            基準電位を差し引く参照点。

        Returns
        -------
        tuple
            ``(figure, axes_array)`` from matplotlib.
        """

        return plot_potential_slices(
            self.result,
            box_min=box_min,
            box_max=box_max,
            grid_n=grid_n,
            xy_z=xy_z,
            yz_x=yz_x,
            xz_y=xz_y,
            softening=softening,
            chunk_size=chunk_size,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            periodic2=periodic2,
            reference_point=reference_point,
        )

    def plot_mesh_source_boxplot(
        self,
        *,
        quantity: str = "charge",
        step: int | None = -1,
        softening: float | None = None,
        self_term: str = "auto",
        showfliers: bool = True,
    ):
        """Plot area-weighted boxplots per mesh source.

        Parameters
        ----------
        quantity : {"charge", "potential"}, default "charge"
            Quantity used in boxplot values.
        step : int or None, default -1
            History batch step used for charge snapshot.
            ``None`` uses final charges from ``charges.csv``.
        softening : float or None, default None
            Softening length in meters (potential mode only).
        self_term : {"auto", "area_equivalent", "exclude", "softened_point"}, default "auto"
            Self-interaction model (potential mode only).
        showfliers : bool, default True
            Whether outlier markers are rendered.

        Returns
        -------
        tuple
            ``(figure, axes)`` from matplotlib.
        """

        return plot_mesh_source_boxplot(
            self.result,
            quantity=quantity,
            step=step,
            softening=softening,
            self_term=self_term,
            showfliers=showfliers,
        )

    def animate_mesh(
        self,
        output_path: str | Path | None = None,
        *,
        quantity: Literal["charge", "potential"] = "charge",
        fps: int = 10,
        frame_stride: int = 1,
        total_frames: int | None = None,
        cmap: str | None = None,
        softening: float | None = None,
        self_term: str = "auto",
        periodic2: Mapping[str, object] | None = None,
        reference_point: Iterable[float] | str | None = "species1_injection_center",
    ) -> Path | FuncAnimation:
        """Animate charge or potential history on the 3D surface mesh.

        Parameters
        ----------
        output_path : str or pathlib.Path or None, default None
            GIF output path. If ``None``, return a live animation object.
        quantity : {"charge", "potential"}, default "charge"
            Quantity used for per-face coloring.
        fps : int, default 10
            Frames per second.
        frame_stride : int, default 1
            Use every ``frame_stride``-th snapshot.
        total_frames : int or None, default None
            Number of evenly sampled frames.
        cmap : str or None, default None
            Matplotlib colormap name.
        softening : float or None, default None
            Softening length in meters (potential mode).
        self_term : {"auto", "area_equivalent", "exclude", "softened_point"}, default "auto"
            Potential self-term model (potential mode).
        periodic2 : mapping or None, default None
            Two-axis periodic setting for potential mode. ``None`` の場合は
            出力ディレクトリ近傍の ``beach.toml`` から自動判定。
        reference_point : iterable of float, {"species1_injection_center"}, or None, default "species1_injection_center"
            基準電位を差し引く参照点。

        Returns
        -------
        pathlib.Path or matplotlib.animation.FuncAnimation
            Output GIF path or in-memory animation.
        """

        return animate_history_mesh(
            self.result,
            output_path=output_path,
            quantity=quantity,
            fps=fps,
            frame_stride=frame_stride,
            total_frames=total_frames,
            cmap=cmap,
            softening=softening,
            self_term=self_term,
            periodic2=periodic2,
            reference_point=reference_point,
        )
