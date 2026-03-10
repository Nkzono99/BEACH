"""History accessor utilities for charge snapshots."""

from __future__ import annotations

from pathlib import Path

import numpy as np


class FortranChargeHistory:
    """Lazy accessor for ``charge_history.csv`` with batch-step indexing.

    Parameters
    ----------
    path : pathlib.Path
        Path to ``charge_history.csv``.
    mesh_nelem : int
        Number of mesh elements. Used to shape per-step charge vectors.
    """

    def __init__(self, path: Path, *, mesh_nelem: int) -> None:
        self.path = Path(path)
        self.mesh_nelem = int(mesh_nelem)
        self._indexed = False
        self._history_matrix: np.ndarray | None = None
        self._batch_indices = np.empty(0, dtype=np.int64)
        self._processed_particles_by_batch = np.empty(0, dtype=np.int64)
        self._rel_change_by_batch = np.empty(0, dtype=float)
        self._step_offsets: dict[int, tuple[int, int]] = {}
        self._step_cache: dict[int, np.ndarray] = {}

    @classmethod
    def from_arrays(
        cls,
        path: Path,
        *,
        mesh_nelem: int,
        history: np.ndarray,
        processed_particles_by_batch: np.ndarray,
        rel_change_by_batch: np.ndarray,
        batch_indices: np.ndarray,
    ) -> FortranChargeHistory:
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
        FortranChargeHistory
            Initialized in-memory history object.
        """

        obj = cls(path, mesh_nelem=mesh_nelem)
        obj._indexed = True
        obj._history_matrix = history
        obj._batch_indices = batch_indices.astype(np.int64, copy=False)
        obj._processed_particles_by_batch = processed_particles_by_batch.astype(
            np.int64, copy=False
        )
        obj._rel_change_by_batch = rel_change_by_batch.astype(float, copy=False)
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
    def processed_particles_by_batch(self) -> np.ndarray:
        """Return processed-particle counts aligned with ``batch_indices``.

        Returns
        -------
        numpy.ndarray
            Processed-particle counts with shape ``(n_snapshots,)``.
        """

        self._ensure_index()
        return self._processed_particles_by_batch

    @property
    def rel_change_by_batch(self) -> np.ndarray:
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
            Full charge history matrix.
        """

        if self._history_matrix is not None:
            return self._history_matrix

        self._ensure_index()
        n_snapshots = self._batch_indices.size
        history = np.zeros((self.mesh_nelem, n_snapshots), dtype=float)
        for col, batch in enumerate(self._batch_indices):
            history[:, col] = self._load_step_from_file(int(batch))
        self._history_matrix = history
        return history

    def get_step(self, step: int) -> np.ndarray:
        """Return per-element charges at one batch step.

        Parameters
        ----------
        step : int
            Batch step to read. ``-1`` selects the latest available step.

        Returns
        -------
        numpy.ndarray
            Per-element charge array with shape ``(mesh_nelem,)``.

        Raises
        ------
        ValueError
            If history is missing/empty, or the requested step is unavailable.
        """

        self._ensure_index()
        if self._batch_indices.size == 0:
            raise ValueError(
                "charge_history.csv is not found or empty. Enable history output and rerun."
            )

        request = int(step)
        if request == -1:
            request = int(self._batch_indices[-1])
        col = self._column_for_step(request)
        if self._history_matrix is not None:
            return self._history_matrix[:, col]
        return self._load_step_from_file(request)

    @staticmethod
    def _parse_row(
        line: str,
    ) -> tuple[int, int, float, int, float]:
        parts = [part.strip() for part in line.split(",", 4)]
        if len(parts) != 5:
            raise ValueError(
                "invalid charge_history.csv row; expected 5 columns: "
                "batch,processed_particles,rel_change,elem_idx,charge_C"
            )
        return (
            int(parts[0]),
            int(parts[1]),
            float(parts[2]),
            int(parts[3]) - 1,
            float(parts[4]),
        )

    def _column_for_step(self, step: int) -> int:
        cols = np.flatnonzero(self._batch_indices == step)
        if cols.size == 0:
            available = [int(v) for v in self._batch_indices]
            raise ValueError(f"step={step} is not found in history. available={available}")
        return int(cols[0])

    def _ensure_index(self) -> None:
        if self._indexed:
            return

        self._indexed = True
        self._step_offsets = {}
        self._step_cache = {}
        self._batch_indices = np.empty(0, dtype=np.int64)
        self._processed_particles_by_batch = np.empty(0, dtype=np.int64)
        self._rel_change_by_batch = np.empty(0, dtype=float)

        if not self.path.exists():
            return
        if self.path.stat().st_size == 0:
            return

        batches: list[int] = []
        processed: list[int] = []
        rel_change: list[float] = []
        seen_batches: set[int] = set()
        current_batch: int | None = None
        current_start: int | None = None

        with self.path.open("r", encoding="utf-8") as stream:
            stream.readline()  # header
            while True:
                pos = stream.tell()
                line = stream.readline()
                if line == "":
                    break
                row = line.strip()
                if not row:
                    continue
                batch, processed_count, rel_value, _, _ = self._parse_row(row)

                if current_batch is None:
                    current_batch = batch
                    current_start = pos
                    seen_batches.add(batch)
                    batches.append(batch)
                    processed.append(processed_count)
                    rel_change.append(rel_value)
                    continue

                if batch == current_batch:
                    continue

                if current_start is None:
                    raise RuntimeError("internal history index state is inconsistent.")
                self._step_offsets[current_batch] = (current_start, pos)

                if batch in seen_batches:
                    raise ValueError(
                        "charge_history.csv must group rows contiguously by batch."
                    )
                seen_batches.add(batch)
                current_batch = batch
                current_start = pos
                batches.append(batch)
                processed.append(processed_count)
                rel_change.append(rel_value)

            end_pos = stream.tell()

        if current_batch is not None and current_start is not None:
            self._step_offsets[current_batch] = (current_start, end_pos)

        self._batch_indices = np.asarray(batches, dtype=np.int64)
        self._processed_particles_by_batch = np.asarray(processed, dtype=np.int64)
        self._rel_change_by_batch = np.asarray(rel_change, dtype=float)

    def _load_step_from_file(self, step: int) -> np.ndarray:
        if step in self._step_cache:
            return self._step_cache[step]

        offsets = self._step_offsets.get(step)
        if offsets is None:
            self._column_for_step(step)
            raise RuntimeError("history step offset is unavailable after index lookup.")

        start, end = offsets
        charges = np.zeros(self.mesh_nelem, dtype=float)
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
                        f"elem_idx out of range in charge_history.csv: {elem_idx + 1}"
                    )
                charges[elem_idx] = charge

        self._step_cache[step] = charges
        return charges
