"""Compare the native finite survey with the independently supplied #38 archive.

Run on a compute node. The archive is an external validation input, not a
runtime dependency or a replacement for BEACH's contract tests.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import hashlib
import subprocess
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, SymLogNorm
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("scanner", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--reuse-native", action="store_true", help="read the existing native_roots.csv instead of rerunning the scanner")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    prefix = "lunar-sheath-solver/results/selection/default-"
    with zipfile.ZipFile(args.archive) as archive:
        queries = list(csv.DictReader(io.StringIO(archive.read(prefix + "queries.csv").decode())))
        reference = defaultdict(list)
        for line in archive.read(prefix + "roots.jsonl").decode().splitlines():
            root = json.loads(line)
            reference[root["query_id"]].append(root)
    inputs = "".join(
        f"{row['query_id']} {row['mach']} {row['temperature_ratio']} "
        f"{row['emission_flux']} {row['wall_field']} 481 1e8\n" for row in queries
    )
    native_path = args.output / "native_roots.csv"
    if not args.reuse_native:
        completed = subprocess.run([str(args.scanner.resolve())], input=inputs, text=True, capture_output=True, check=True)
        native_path.write_text(completed.stdout)
    actual = defaultdict(list)
    statuses = {}
    for row in csv.DictReader(io.StringIO(native_path.read_text())):
        index = int(row["query_id"])
        statuses[index] = int(row["status"])
        if int(row["root_count"]):
            actual[index].append(row)
    differences = []
    max_error = 0.0
    for query in queries:
        index = int(query["query_id"])
        expected = reference[index]
        found = actual[index]
        expected_status = {"solutions": 0, "no_physical_solution": 2, "unresolved": 3}[query["status"]]
        if statuses.get(index) != expected_status:
            differences.append({"query": query, "native_status": statuses.get(index)})
        expected_types = sorted("B0" if r["branch"] == "B" and r["phi_wall"] == 0 else r["branch"] for r in expected)
        if sorted(r["branch"] for r in found) != expected_types:
            differences.append({"query": query, "native": found, "reference": expected})
            continue
        for r in expected:
            kind = "B0" if r["branch"] == "B" and r["phi_wall"] == 0 else r["branch"]
            best = min((s for s in found if s["branch"] == kind), key=lambda s: abs(float(s["phi_h"])-r["phi_wall"]))
            error = max(abs(float(best[key])-r[ref]) / max(1.0, abs(r[ref])) for key, ref in (
                ("phi_h", "phi_wall"), ("phi_min", "phi_min"), ("amplitude", "electron_amplitude"),
                ("q", "escaping_photo_flux"), ("current", "current")))
            max_error = max(max_error, error)
            if error > 2e-6:
                differences.append({"query": query, "error": error, "native": best, "reference": r})
    report = {
        "queries": len(queries), "native_status_counts": dict(Counter(statuses.values())),
        "native_roots": sum(map(len, actual.values())), "differences": differences,
        "source_archive_sha256": hashlib.sha256(args.archive.read_bytes()).hexdigest(),
        "native_csv_sha256": hashlib.sha256(native_path.read_bytes()).hexdigest(),
        "maximum_scaled_root_error": max_error,
        "note": "Finite numerical surveys; root agreement does not establish dynamic stability or completeness.",
    }
    (args.output / "comparison.json").write_text(json.dumps(report, indent=2))
    # The first 980 independent queries are the issue's rectangular E-G slice.
    section = queries[:980]
    fields = sorted({float(r["wall_field"]) for r in section})
    emissions = sorted({float(r["emission_flux"]) for r in section})
    classes = np.full((len(emissions), len(fields)), np.nan)
    counts = classes.copy()
    potential = classes.copy()
    maximum = classes.copy()
    types = classes.copy()
    deep = classes.copy()
    type_sets = sorted({"+".join(sorted(r["branch"] for r in actual[int(row["query_id"])])) or "none" for row in section})
    for row in section:
        i = emissions.index(float(row["emission_flux"]))
        j = fields.index(float(row["wall_field"]))
        index = int(row["query_id"])
        classes[i, j] = statuses[index]
        counts[i, j] = len(actual[index])
        types[i, j] = type_sets.index("+".join(sorted(r["branch"] for r in actual[index])) or "none")
        deep[i, j] = any(float(r["phi_min"]) < -1e4 for r in actual[index])
        if actual[index]:
            potential[i, j] = min(float(r["phi_h"]) for r in actual[index])
            maximum[i, j] = max(float(r["phi_h"]) for r in actual[index])
    report["EG_slice"] = dict(status_counts=dict(Counter(statuses[int(row["query_id"])] for row in section)),
                              roots=sum(len(actual[int(row["query_id"])]) for row in section))
    (args.output / "comparison.json").write_text(json.dumps(report, indent=2))
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=True)
    for panel, (ax, data, title) in enumerate(zip(axes.flat, (classes, types, counts, potential, maximum, deep), (
        "Classification", "Detected Type set", "Detected root count", "Minimum potential / Te", "Maximum potential / Te", "Depth > 10000 Te"))):
        kwargs = {}
        ticks = None
        labels = None
        if panel in (3, 4):
            kwargs = dict(norm=SymLogNorm(linthresh=0.01), cmap="coolwarm")
        elif panel == 0:
            kwargs = dict(norm=BoundaryNorm([-0.5, 1, 2.5, 3.5], 3), cmap=ListedColormap(["#548c9c", "#d4d5da", "#dc9a37"]))
            ticks, labels = [0, 2, 3], ["solutions", "excluded", "unresolved"]
        elif panel == 1:
            kwargs = dict(norm=BoundaryNorm(np.arange(len(type_sets) + 1) - 0.5, len(type_sets)),
                          cmap=plt.get_cmap("tab10", len(type_sets)))
            ticks, labels = range(len(type_sets)), type_sets
        elif panel == 2:
            ticks = [0, 1, 2]
        else:
            ticks, labels = [0, 1], ["no", "yes"]
        plot = ax.imshow(data, origin="lower", aspect="auto", interpolation="nearest", **kwargs)
        bar = fig.colorbar(plot, ax=ax, shrink=0.8, ticks=ticks)
        if labels is not None:
            bar.set_ticklabels(labels)
        ax.set_title(title, fontsize=10)
        ax.set_xticks([0, 17, 34], [f"{fields[j]:g}" for j in [0, 17, 34]])
        ax.set_yticks([0, len(emissions)-1], [f"{emissions[j]:g}" for j in [0, len(emissions)-1]])
        ax.set_xlabel("E_H (sample index spacing)")
        ax.set_ylabel("G (sample index spacing)")
    fig.savefig(args.output / "existence.png", dpi=170)
    print(json.dumps({k: v for k, v in report.items() if k != "differences"}, indent=2))
    print(f"difference_count={len(differences)}")
    if differences:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
