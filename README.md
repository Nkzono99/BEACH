# bemtracer

Boundary Element Method (BEM) surface-charging + test-particle prototype in Python.

This repository implements:
- Triangular boundary mesh with per-element charge (insulator accumulation mode)
- Test particle pusher (Boris, E+B ready)
- Segment–triangle collision (first-hit)
- Batched charge deposition (`npcls_per_step`) with electric-field cache rebuild
- Hybrid electric field: far-field centroid point-charge + near-field triangle subdivision correction (difference form)

## Quickstart

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
source .venv/bin/activate
pip install -e .
python examples/minimal_run.py
```

## Files
- `SPEC.md` : specification (v0.1)
- `src/bem_charge_tracer/core.py` : main implementation
- `examples/minimal_run.py` : minimal runnable example

## Notes
- v0.1 focuses on correctness and extensibility. Performance optimizations (BVH/Octree/FMM) are out of scope for now.
- Offset evaluation (δ-shift) is intentionally *not* used in v0.1.
