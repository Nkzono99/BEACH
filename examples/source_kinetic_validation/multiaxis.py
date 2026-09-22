"""Compare finite M-E and tau-G slices against the supplied Python solver.

Put the extracted handoff's src directory on PYTHONPATH. This validation-only
script requires scipy and matplotlib; neither is a BEACH runtime dependency.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import numpy as np

from lunar_sheath.extended import solve_extended


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scanner", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    fields = np.r_[-np.geomspace(0.03, 2, 5)[::-1], np.geomspace(0.03, 2, 5)]
    machs = np.geomspace(0.5, 20, 9)
    taus = np.geomspace(0.03, 2, 9)
    emissions = np.geomspace(0.005, 5, 10)
    queries = [(m, 0.2, 0.3, e) for m in machs for e in fields]
    queries += [(10, t, g, 0.1) for t in taus for g in emissions]
    # The zero-drift Type A fixture is shared with the legacy model because
    # both models use the same connected orbits on either side of a minimum.
    queries += [(10, 0.2, 3.3, 2.440435545625)]
    native = defaultdict(list)
    statuses = {}
    for points in (481, 961):
        inputs = "".join(f"{i} {m} {t} {g} {e} {points} 1e8\n" for i, (m, t, g, e) in enumerate(queries))
        output = subprocess.run([str(args.scanner.resolve())], input=inputs, text=True,
                                capture_output=True, check=True).stdout
        (args.output / f"native-{points}.csv").write_text(output)
        for row in csv.DictReader(io.StringIO(output)):
            index = int(row["query_id"])
            statuses[points, index] = int(row["status"])
            if int(row["root_count"]):
                native[points, index].append(row)
    differences, refinements, reference = [], [], []
    maximum_error = 0.0
    for i, (m, t, g, e) in enumerate(queries):
        result = solve_extended(m, t, g, e, electron_model="source_kinetic")
        reference.append(dict(query_id=i, parameters=[m, t, g, e], **result.to_dict()))
        expected = sorted(r.branch for r in result.roots)
        found = native[481, i]
        if expected != sorted(r["branch"] for r in found):
            differences.append(dict(query_id=i, expected=expected, native=found))
        else:
            for root in result.roots:
                best = min((r for r in found if r["branch"] == root.branch),
                           key=lambda r: abs(float(r["phi_h"]) - root.phi_wall))
                error = max(abs(float(best[key]) - value) / max(1.0, abs(value)) for key, value in (
                    ("phi_h", root.phi_wall), ("phi_min", root.phi_min),
                    ("amplitude", root.electron_amplitude), ("q", root.escaping_photo_flux),
                    ("current", root.current)))
                maximum_error = max(maximum_error, error)
                if error > 2e-6:
                    differences.append(dict(query_id=i, error=error))
        expected_status = {"solutions": 0, "no_physical_solution": 2, "unresolved": 3}[result.status]
        if statuses[481, i] != expected_status:
            differences.append(dict(query_id=i, expected_status=expected_status, native_status=statuses[481, i]))
        if statuses[481, i] != statuses[961, i] or sorted(r["branch"] for r in found) != sorted(
            r["branch"] for r in native[961, i]
        ):
            refinements.append(i)
    (args.output / "reference.json").write_text(json.dumps(reference, indent=2))
    report = dict(queries=len(queries), differences=differences, grid_dependent_queries=refinements,
                  maximum_scaled_root_error=maximum_error,
                  note="Finite slices; sampling does not define bifurcation curves or certify completeness.")
    (args.output / "comparison.json").write_text(json.dumps(report, indent=2))
    for name, indices, xs, ys, xlabel, ylabel in (
        ("ME", range(90), fields, machs, "E_H", "M"),
        ("tauG", range(90, 180), emissions, taus, "G", "tau"),
    ):
        sets = sorted({"+".join(sorted(r["branch"] for r in native[481, i])) or "none" for i in indices})
        data = [[], [], [], [], [], []]
        for i in indices:
            roots = native[481, i]
            data[0].append(statuses[481, i])
            data[1].append(sets.index("+".join(sorted(r["branch"] for r in roots)) or "none"))
            data[2].append(len(roots))
            potentials = [float(r["phi_h"]) for r in roots]
            data[3].append(min(potentials, default=np.nan))
            data[4].append(max(potentials, default=np.nan))
            data[5].append(i in refinements)
        fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
        for panel, (ax, values, title) in enumerate(zip(axes.flat, data, (
            "Status (0: roots, 2: excluded, 3: unresolved)", "Type set", "Root count",
            "Minimum phi_H / Te", "Maximum phi_H / Te", "Root count/Type changes at 961 points",
        ))):
            kw = dict(norm=SymLogNorm(linthresh=0.01), cmap="coolwarm") if panel in (3, 4) else {}
            plot = ax.imshow(np.array(values).reshape(len(ys), len(xs)), origin="lower", aspect="auto", **kw)
            bar = fig.colorbar(plot, ax=ax, shrink=0.8)
            if panel == 1:
                bar.set_ticks(range(len(sets)), labels=sets)
            ax.set_xticks(range(len(xs)), [f"{x:.2g}" for x in xs], rotation=45)
            ax.set_yticks(range(len(ys)), [f"{y:.2g}" for y in ys])
            ax.set(xlabel=xlabel + " (sample spacing)", ylabel=ylabel + " (sample spacing)", title=title)
        fig.savefig(args.output / f"{name}.png", dpi=170)
        plt.close(fig)
    print(json.dumps(report, indent=2))
    if differences:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
