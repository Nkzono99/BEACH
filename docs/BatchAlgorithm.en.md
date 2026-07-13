title: Batch coupling algorithm

Lang: [English](BatchAlgorithm.en.md) | [日本語](BatchAlgorithm.md)

# Batch coupling algorithm

BEACH couples the field produced by surface charge with particle transport in batches. This page describes the
update order from initialization through charge commit and the numerical consequences of that ordering.

## Initialization

At startup, BEACH approximately follows this sequence:

1. Read `beach.toml`, normalize high-level notation, and validate physical constraints.
2. Build the triangular mesh and surface models.
3. Initialize particle species and injection state.
4. Prepare the field solver and periodic operators.
5. Set surface charge, statistics, and RNG state for a new run or from a checkpoint.
6. Prepare output and history writers.

A resume stops if the model, mesh, or species fingerprint differs from the checkpoint. See
[Run a simulation](Execution.en.html#resume-a-run).

## Batch loop

The following sequence repeats until `sim.batch_count` is reached.

```text
q_elem(batch start)
        │
        ├─ update field-solver state
        ├─ update required potential and boundary profiles
        └─ prepare particle sources
                    ↓
             generate and track particles
                    ↓
       count absorbed, escaped, and unresolved outcomes
                    ↓
               reduce delta_q_elem
                    ↓
       q_elem ← q_elem + delta_q_elem
                    ↓
       update surface models, statistics, and history
```

### 1. Update the field

The direct, treecode, or FMM solver state is updated from `q_elem` at the start of the batch. Under `periodic2`,
near images, the cached nonzero operator, zero mode, and outer model are combined according to the selected
contract.

### 2. Generate particles

Volume-seed, reservoir, photoelectron, and other sources are evaluated for each species. A non-integer reservoir
macro-particle count is retained as a residual for the next batch.

### 3. Track particles

Each particle advances for at most `sim.max_step` steps. Every step evaluates the field, updates velocity and
position, and selects the first box-boundary or triangle event.

- Surface collision: absorb the particle and add its charge to the hit element's `delta_q_elem`.
- Open boundary: decide escape or return according to the boundary model.
- `max_step`: record `survived_max_step` without reclassifying it as absorption or escape.

See [Particle tracking and collision](ParticleTrackingCollision.en.html) for integration and collision details.

### 4. Commit charge

Particle tracking does not modify `q_elem`. It accumulates element-local `delta_q_elem`, which is committed only
after all particles are processed. Particles in one batch therefore do not interact through charge deposited
within that same batch.

An insulator retains deposited charge on the hit element. For a floating conductor, charge is redistributed
between conductor elements after commit to approximate an equipotential condition.

### 5. Update statistics and history

Batch statistics and history are updated from the committed state. `tol_rel` is a monitoring metric for surface
charge change; it is not an early-stop condition. A normal run advances to the configured `sim.batch_count`.

## Parallel-execution contract

OpenMP processes particles in parallel and reduces thread-local statistics and charge deltas at the end of the
batch. MPI distributes particles across ranks and collectively combines element charge deltas and statistics.

Both paths preserve these rules:

- Build the field from surface charge at the start of the batch.
- Do not update shared `q_elem` during particle processing.
- Give every rank the same surface charge after commit.
- Keep injection, absorption, escape, and unresolved outcomes distinct in the charge ledger.

## Distinct time controls

| Value | Role |
| --- | --- |
| `sim.dt` | One particle-trajectory step |
| `sim.max_step` | Maximum steps tracked for one particle |
| `batch_duration` | Physical time and particle supply represented by one batch |
| `sim.batch_count` | Cumulative number of batches to run |

Reducing `dt` and varying `batch_duration` test different convergence axes. See
[`batch_duration` stability and steady value](BatchDurationStability.en.html) for the latter.
