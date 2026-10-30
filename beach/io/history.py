"""Lazy per-element batch histories for saved charge and potential snapshots."""

from __future__ import annotations

from math import isfinite
from pathlib import Path

import numpy as np


class BatchHistory:
    """Lazy accessor for per-element CSV values with batch-step indexing.

    Parameters
    ----------
    path : pathlib.Path
        Path to ``charge_history.csv`` or ``potential_history.csv``.
    mesh_nelem : int
        Number of mesh elements. Used to shape per-step vectors.
    value_column : str, default "charge_C"
        Value column to read, for example ``"charge_C"`` or ``"potential_V"``.
        Histories without processed-particle/relative-change metadata expose
        ``None`` for those two properties.
    """

    def __init__(
        self, path: Path, *, mesh_nelem: int, value_column: str = "charge_C"
    ) -> None:
        self.path = Path(path)
        self.value_column = value_column
        self._columns: dict[str, int] = {}
        self.mesh_nelem = int(mesh_nelem)
        self._indexed = False
        self._history_matrix: np.ndarray | None = None
        self._batch_indices = np.empty(0, dtype=np.int64)
        self._processed_particles_by_batch: np.ndarray | None = None
        self._rel_change_by_batch: np.ndarray | None = None
        self._step_offsets: dict[int, tuple[int, int]] = {}
        self._step_columns: dict[int, int] | None = None
        self._step_cache: dict[int, np.ndarray] = {}

    @classmethod
    def from_arrays(
        cls,
        path: Path,
        *,
        mesh_nelem: int,
        history: np.ndarray,
        batch_indices: np.ndarray,
        processed_particles_by_batch: np.ndarray | None = None,
        rel_change_by_batch: np.ndarray | None = None,
        value_column: str = "charge_C",
    ) -> BatchHistory:
        """Construct an in-memory history object without reading CSV from disk.

        Parameters
        ----------
        path : pathlib.Path
            Virtual history path metadata.
        mesh_nelem : int
            Number of mesh elements.
        history : numpy.ndarray
            Charge history matrix with shape ``(mesh_nelem, n_snapshots)``.
        processed_particles_by_batch : numpy.ndarray
            Processed-particle counts for each snapshot.
        rel_change_by_batch : numpy.ndarray
            Relative-change values for each snapshot.
        batch_indices : numpy.ndarray
            Batch index values for each snapshot.

        Returns
        -------
        BatchHistory
            Initialized in-memory history object.
        """

        obj = cls(path, mesh_nelem=mesh_nelem, value_column=value_column)
        obj._indexed = True
        obj._history_matrix = history
        obj._batch_indices = batch_indices.astype(np.int64, copy=False)
        obj._processed_particles_by_batch = (
            None
            if processed_particles_by_batch is None
            else processed_particles_by_batch.astype(np.int64, copy=False)
        )
        obj._rel_change_by_batch = (
            None
            if rel_change_by_batch is None
            else rel_change_by_batch.astype(float, copy=False)
        )
        return obj

    @property
    def has_data(self) -> bool:
        """Return whether history snapshots exist.

        Returns
        -------
        bool
            ``True`` when at least one snapshot is available.
        """

        if self._history_matrix is not None:
            return self._history_matrix.size > 0
        self._ensure_index()
        return self._batch_indices.size > 0

    @property
    def batch_indices(self) -> np.ndarray:
        """Return batch indices corresponding to available history snapshots.

        Returns
        -------
        numpy.ndarray
            Batch indices with shape ``(n_snapshots,)``.
        """

        self._ensure_index()
        return self._batch_indices

    @property
    def processed_particles_by_batch(self) -> np.ndarray | None:
        """Return processed-particle counts aligned with ``batch_indices``.

        Returns
        -------
        numpy.ndarray
            Processed-particle counts with shape ``(n_snapshots,)``.
        """

        self._ensure_index()
        return self._processed_particles_by_batch

    @property
    def rel_change_by_batch(self) -> np.ndarray | None:
        """Return per-batch relative charge-change metrics from history.

        Returns
        -------
        numpy.ndarray
            Relative-change values with shape ``(n_snapshots,)``.
        """

        self._ensure_index()
        return self._rel_change_by_batch

    def __len__(self) -> int:
        return int(self.batch_indices.size)

    def __getitem__(self, step: int) -> np.ndarray:
        return self.get_step(step)

    def as_array(self) -> np.ndarray:
        """Materialize full history as ``(mesh_nelem, n_snapshots)`` array.

        Returns
        -------
        numpy.ndarray
            Full value history matrix.
        """

        if self._history_matrix is not None:
            return self._history_matrix

        self._ensure_index()
        n_snapshots = self._batch_indices.size
        history = np.empty((self.mesh_nelem, n_snapshots), dtype=float)
        for col, batch in enumerate(self._batch_indices):
            history[:, col] = self._load_step_from_file(int(batch))
            # The matrix now owns this snapshot; do not retain a second copy.
            self._step_cache.pop(int(batch))
        self._history_matrix = history
        return history

    def get_step(self, step: int) -> np.ndarray:
        """Return per-element values at one batch step.

        Parameters
        ----------
        step : int
            Batch step to read. ``-1`` selects the latest available step.

        Returns
        -------
        numpy.ndarray
            Per-element value array with shape ``(mesh_nelem,)``.

        Raises
        ------
        ValueError
            If history is missing/empty, or the requested step is unavailable.
        """

        self._ensure_index()
        if self._batch_indices.size == 0:
            raise ValueError(
                f"{self.path.name} is not found or empty. Enable history output and rerun."
            )

        request = int(step)
        if request == -1:
            request = int(self._batch_indices[-1])
        col = self._column_for_step(request)
        if self._history_matrix is not None:
            return self._history_matrix[:, col]
        return self._load_step_from_file(request)

    def _parse_row(self, line: str) -> tuple[int, int | None, float | None, int, float]:
        parts = line.split(",")
        if len(parts) != len(self._columns):
            raise ValueError(
                f"invalid {self.path.name} row; expected {len(self._columns)} columns."
            )
        columns = self._columns
        return (
            int(parts[columns["batch"]]),
            int(parts[columns["processed_particles"]])
            if "processed_particles" in columns
            else None,
            float(parts[columns["rel_change"]]) if "rel_change" in columns else None,
            int(parts[columns["elem_idx"]]) - 1,
            float(parts[columns[self.value_column]]),
        )

    def _column_for_step(self, step: int) -> int:
        if self._step_columns is None:
            columns: dict[int, int] = {}
            for column, batch in enumerate(self._batch_indices):
                columns.setdefault(int(batch), column)
            self._step_columns = columns
        column = self._step_columns.get(step)
        if column is None:
            available = [int(v) for v in self._batch_indices]
            raise ValueError(
                f"step={step} is not found in history. available={available}"
            )
        return column

    def _select_columns(self, columns: np.ndarray) -> np.ndarray:
        """Return the ``(mesh_nelem, len(columns))`` values of selected snapshots.

        Unlike :meth:`get_step`, snapshots read here from the CSV are not kept in
        the step cache, so selecting a few frames of a long history holds only
        those frames. Cached or materialized snapshots, including edits made to
        arrays returned by :meth:`get_step`, are used as they are.
        """

        self._ensure_index()
        if self._history_matrix is not None:
            return self._history_matrix[:, columns]
        selected = np.empty((self.mesh_nelem, len(columns)), dtype=float)
        for position, column in enumerate(columns):
            step = int(self._batch_indices[int(column)])
            cached = self._step_cache.get(step)
            selected[:, position] = (
                self._read_step(step) if cached is None else cached
            )
        return selected

    def _ensure_index(self) -> None:
        if self._indexed:
            return

        step_offsets: dict[int, tuple[int, int]] = {}
        batches: list[int] = []
        processed: list[int | None] = []
        rel_change: list[float | None] = []
        current_batch: int | None = None
        current_start = 0
        current_processed: int | None = None
        current_rel_change: float | None = None

        if self.path.exists() and self.path.stat().st_size > 0:
            seen_elem = np.zeros(self.mesh_nelem, dtype=bool)
            with self.path.open("r", encoding="utf-8") as stream:
                header = stream.readline().strip().split(",")
                self._columns = {name.strip(): i for i, name in enumerate(header)}
                if not {"batch", "elem_idx", self.value_column}.issubset(self._columns):
                    raise ValueError(f"{self.path.name} is missing required columns.")
                while True:
                    pos = stream.tell()
                    line = stream.readline()
                    if line == "":
                        break
                    row = line.strip()
                    if not row:
                        continue
                    batch, processed_count, rel_value, elem_idx, charge = (
                        self._parse_row(row)
                    )

                    if current_batch is None or batch != current_batch:
                        if current_batch is not None:
                            self._validate_complete_batch(current_batch, seen_elem)
                            step_offsets[current_batch] = (current_start, pos)
                            if batch <= current_batch:
                                raise ValueError(
                                    f"{self.path.name} batch={batch} is not greater than "
                                    f"previous batch={current_batch}."
                                )

                        current_batch = batch
                        current_start = pos
                        current_processed = processed_count
                        current_rel_change = rel_value
                        seen_elem.fill(False)
                        batches.append(batch)
                        processed.append(processed_count)
                        rel_change.append(rel_value)

                    elem_number = elem_idx + 1
                    if elem_idx < 0 or elem_idx >= self.mesh_nelem:
                        raise ValueError(
                            f"{self.path.name} batch={batch} has out-of-range "
                            f"elem_idx={elem_number}; expected 1..{self.mesh_nelem}."
                        )
                    if seen_elem[elem_idx]:
                        raise ValueError(
                            f"{self.path.name} batch={batch} has duplicate elem_idx={elem_number}."
                        )
                    if not isfinite(charge):
                        raise ValueError(
                            f"{self.path.name} batch={batch} elem_idx={elem_number} "
                            f"has non-finite {self.value_column}."
                        )
                    if rel_value is not None and not isfinite(rel_value):
                        raise ValueError(
                            f"{self.path.name} batch={batch} elem_idx={elem_number} "
                            "has non-finite rel_change."
                        )
                    if processed_count != current_processed:
                        raise ValueError(
                            f"{self.path.name} batch={batch} has inconsistent "
                            f"processed_particles at elem_idx={elem_number}; "
                            f"expected {current_processed}, got {processed_count}."
                        )
                    if rel_value != current_rel_change:
                        raise ValueError(
                            f"{self.path.name} batch={batch} has inconsistent "
                            f"rel_change at elem_idx={elem_number}; "
                            f"expected {current_rel_change}, got {rel_value}."
                        )
                    seen_elem[elem_idx] = True

                end_pos = stream.tell()

            if current_batch is not None:
                self._validate_complete_batch(current_batch, seen_elem)
                step_offsets[current_batch] = (current_start, end_pos)

        self._step_offsets = step_offsets
        self._step_columns = None
        self._step_cache = {}
        self._batch_indices = np.asarray(batches, dtype=np.int64)
        self._processed_particles_by_batch = (
            np.asarray(processed, dtype=np.int64)
            if "processed_particles" in self._columns
            else None
        )
        self._rel_change_by_batch = (
            np.asarray(rel_change, dtype=float)
            if "rel_change" in self._columns
            else None
        )
        self._indexed = True

    def _validate_complete_batch(self, batch: int, seen_elem: np.ndarray) -> None:
        missing = ~seen_elem
        missing_count = int(np.count_nonzero(missing))
        if missing_count == 0:
            return
        first_missing = int(np.argmax(missing)) + 1
        raise ValueError(
            f"{self.path.name} batch={batch} is incomplete: "
            f"missing elem_idx={first_missing} ({missing_count} missing)."
        )

    def _load_step_from_file(self, step: int) -> np.ndarray:
        if step in self._step_cache:
            return self._step_cache[step]
        charges = self._read_step(step)
        self._step_cache[step] = charges
        return charges

    def _read_step(self, step: int) -> np.ndarray:
        start, end = self._step_offsets[step]
        charges = np.empty(self.mesh_nelem, dtype=float)
        with self.path.open("r", encoding="utf-8") as stream:
            stream.seek(start)
            while stream.tell() < end:
                line = stream.readline()
                if line == "":
                    break
                row = line.strip()
                if not row:
                    continue
                batch, _, _, elem_idx, charge = self._parse_row(row)
                if batch != step:
                    continue
                if elem_idx < 0 or elem_idx >= self.mesh_nelem:
                    raise ValueError(
                        f"elem_idx out of range in {self.path.name}: {elem_idx + 1}"
                    )
                charges[elem_idx] = charge
        return charges
