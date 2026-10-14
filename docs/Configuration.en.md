title: Create and Validate beach.toml

Lang: [English](Configuration.en.md) | [日本語](Configuration.md)

# Create and Validate `beach.toml`

The grouped layout is the current authoring format. See [configuration groups and migration](GroupedConfiguration.en.html) for new keys and conversion examples. The flat notation below remains supported during 1.x and is scheduled for removal at 2.0.


This document describes the directly edited `beach.toml` file and the `beachx config` helper commands.
See [Design a Simulation Case](ConfigurationRecipes.en.html) for choosing meshes, particle sources, boundary conditions,
and the rest of the physical setup.

- The Fortran runtime `beach` reads `beach.toml` directly.
- `beachx config init` creates the official multi-particle, 20-batch tutorial `beach.toml`.
- See [Input Parameters Reference](Parameters.en.html) for every key, including coordinate and placement helpers evaluated while loading.

## 1. Basic Flow

```bash
mkdir beach-tutorial
cd beach-tutorial

beachx config init beach.toml
$EDITOR beach.toml
beachx lint beach.toml
beach beach.toml
```

Run `beachx lint` before a simulation and confirm `status=ok` before starting `beach`.

Write `box_origin` / `box_size`, `inject_region_mode`, `mesh.groups`, and related keys directly in TOML. See
[Coordinate and placement helper parameters](Parameters.en.html#coordinate-and-placement-helper-parameters) for the coordinates
or dimensions they calculate and any explicit values they replace.

## 2. Commands

### 2.1 `init`

Create a new `beach.toml`. The command fails if the file already exists.

```bash
beachx config init
beachx config init run.toml
beachx config init --force
```

The generated file is identical to
[`examples/tutorial_insulator.toml`](https://github.com/Nkzono99/BEACH/blob/main/examples/tutorial_insulator.toml).
It is the official beginner case that launches 200 `volume_seed` macro-electrons per batch toward an insulating plane
and follows the charge distribution and its feedback for 20 batches. It uses the easier-to-interpret
`field_solver="direct"` and `[field_boundary] mode="free"`. It does not include periodic boundaries, an ion species,
or `photo_raycast`.

### 2.2 `lint`

Run TOML parsing, the packaged JSON Schema (keys, types, enumerations, ranges), and the coordinate and placement
normalization, then run the Fortran preflight (combination constraints and derived values) through `beach --check-config`.
Success prints `checks=toml,schema,preflight` and `status=ok`. When the preflight rejects the input, lint prints the Fortran diagnostic unchanged.

```bash
beachx lint beach.toml
beachx lint run.toml --schema schemas/beach.schema.json
beachx lint run.toml --beach /path/to/beach
```

`--schema` adds constraints to the packaged BEACH schema. The supplied schema applies before and after normalization;
it cannot disable or relax the normal BEACH validation rules.
The `beach` executable is found through `--beach`, the `BEACH_CONFIG_CHECK_EXE` environment variable, and `PATH`, in that order; lint fails when none is available.

### 2.3 `beach --check-config`

The Fortran preflight that `lint` calls. It applies the same defaults, loading, and normalization as a run, then the
combination constraints and derived values. The configuration path is required. Ordinary use only needs `beachx lint`.

```bash
beach --check-config beach.toml
```

Success returns exit code 0 and prints the following. Detected input errors return a nonzero exit code.

```text
config=beach.toml
checks=toml,semantic
status=ok
```

The check does not start a simulation or create result files. It does not validate the contents of external OBJ files,
response tables, or checkpoints, or the numerical and physical validity of a run. Continue with the normal
`beach beach.toml` command to check external-file loading and model initialization.

### 2.4 `diff`

Semantically compare two configs. By default, coordinate and placement helpers are converted to physical coordinates and sizes
before comparison.

```bash
beachx config diff left.toml right.toml
beachx config diff --raw left.toml right.toml
```

## 3. Schema

Place a `#:schema` directive at the beginning of `beach.toml` to enable completion and type checking in
VS Code extensions such as Even Better TOML / Taplo.

```toml
#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.schema.json
```

For a local checkout:

```toml
#:schema ../schemas/beach.schema.json
```

The BEACH Fortran parser does not accept ordinary `key = value` entries before the first section, so use the
comment directive rather than `"$schema" = "..."`.

## 4. Common Mistakes

### 4.1 Top-Level Key Placement

Settings belong under the [public TOML sections](Parameters.en.html#toml-hierarchy-and-section-list).
Ordinary keys before the first section, or unknown top-level sections, fail validation or Fortran loading.

### 4.2 Specifying the Same Coordinate Two Ways

Writing the same coordinate in two forms, such as `box_origin` / `box_size` together with `box_min` / `box_max`, fails
validation. `size_mode="box_fraction"` and group scaling intentionally replace corresponding dimensions with calculated values;
the affected keys are listed in [Input Parameters Reference](Parameters.en.html#coordinate-and-placement-helper-parameters).

### 4.3 Checking Before Running

Run `beachx lint beach.toml` before a simulation, then start `beach beach.toml` after it succeeds.
Lint calls `beach --check-config` after the schema check, so the combination diagnostics also come from lint.
Numbers must be finite, and integer fields must use
integers. The Fortran reader reports strings exceeding their destination length as errors instead of truncating them.
