"""Read per-batch secondary emission diagnostics."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from ._summary_values import _parse_nonnegative_finite_float, _parse_nonnegative_int
from .types import SecondaryEmissionHistoryEntry

_HEADER = ("batch", "incident_species", "hits", "emitted_count", "emitted_charge_c",
           "mean_yield", "generation_capped_count")


def load_secondary_emission_history(path: Path) -> tuple[SecondaryEmissionHistoryEntry, ...] | None:
    if not path.exists():
        return None
    entries = []
    previous_batch = 0
    species_in_batch: set[str] = set()
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != _HEADER:
            raise ValueError("secondary_emission_history.csv header does not match the 7-column contract.")
        for line, row in enumerate(reader, 2):
            if None in row or any(row.get(key) is None for key in _HEADER):
                raise ValueError(f"secondary_emission_history.csv line {line} must have exactly 7 columns.")
            batch = _parse_nonnegative_int(row["batch"], key="batch")
            species = row["incident_species"]
            if batch <= 0 or batch < previous_batch or not species:
                raise ValueError("secondary_emission_history.csv requires positive ordered batches and incident species.")
            if batch != previous_batch:
                species_in_batch.clear()
            if species in species_in_batch:
                raise ValueError("secondary_emission_history.csv has duplicate batch/incident_species rows.")
            species_in_batch.add(species)
            previous_batch = batch
            # Emitted particle charge is negative; the reaction charge is its opposite.
            charge = float(row["emitted_charge_c"])
            if not math.isfinite(charge) or charge > 0:
                raise ValueError("secondary_emission_history.csv emitted_charge_c must be finite and nonpositive.")
            entries.append(SecondaryEmissionHistoryEntry(
                batch=batch, incident_species=species,
                hits=_parse_nonnegative_int(row["hits"], key="hits"),
                emitted_count=_parse_nonnegative_int(row["emitted_count"], key="emitted_count"),
                emitted_charge_c=charge,
                mean_yield=_parse_nonnegative_finite_float(row["mean_yield"], key="mean_yield"),
                generation_capped_count=_parse_nonnegative_int(row["generation_capped_count"], key="generation_capped_count"),
            ))
    return tuple(entries)
