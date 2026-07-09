title: beachx config / High-Level Notation Guide

Lang: [English](Configuration.en.md) | [日本語](Configuration.md)

# `beachx config` / High-Level Notation Guide

This document describes the directly edited `beach.toml` file and the `beachx config` helper commands.

- The Fortran runtime `beach` reads `beach.toml` directly and normalizes high-level notation while loading it.
- `beachx config init` creates a small runnable `beach.toml`.
- See [Input Parameters Reference](Parameters.en.html) for the final keys read by the Fortran runtime.

## 1. Basic Flow

```bash
mkdir run_periodic2
cd run_periodic2

beachx config init
$EDITOR beach.toml
beachx lint beach.toml
beach beach.toml
```

No separate expansion file is needed. The Fortran parser normalizes `box_origin` / `box_size`,
`inject_region_mode`, `mesh.groups`, and related authoring keys before validation.

## 2. Commands

### 2.1 `init`

Create a new `beach.toml`. The command fails if the file already exists.

```bash
beachx config init
beachx config init run.toml
beachx config init --force
```

The initial file is a small validation case with two-periodic-axis FMM, electron/ion volume seeds,
`photo_raycast` electron emission, a plane mesh, and standard output settings.

### 2.2 `lint`

Run TOML parsing, JSON Schema validation, high-level notation checks, and BEACH-specific constraints together.

```bash
beachx lint beach.toml
beachx lint run.toml --schema schemas/beach.schema.json
```

### 2.3 `validate`

Read `beach.toml` and validate high-level notation consistency plus known constraints of the final settings.
Use `beachx lint` when JSON Schema validation should also be included.

```bash
beachx config validate
beachx config validate run.toml
```

### 2.4 `diff`

Semantically compare two configs. By default, high-level notation is normalized before comparison.

```bash
beachx config diff left.toml right.toml
beachx config diff --raw left.toml right.toml
```

## 3. High-Level Notation

High-level notation is a convenience layer for writing coordinate intent directly in `beach.toml`.
The Fortran parser converts it to ordinary numerical keys while loading the file.

### 3.1 Simulation Box

Use `sim.box_origin` and `sim.box_size` to derive `box_min` / `box_max`.

```toml
[sim]
use_box = true
box_origin = [0.0, 0.0, -1.0]
box_size = [1.0, 1.0, 2.0]
```

Resolved while loading:

```toml
[sim]
use_box = true
box_min = [0.0, 0.0, -1.0]
box_max = [1.0, 1.0, 1.0]
```

### 3.2 Injection Region

For `reservoir_face` and `photo_raycast`, the injection region can be expressed as fractions on the face.

```toml
[[particles.species]]
source_mode = "reservoir_face"
inject_face = "z_high"
inject_region_mode = "face_fraction"
uv_low = [0.25, 0.25]
uv_high = [0.75, 0.75]
```

These values are resolved into `pos_low` / `pos_high` while loading.

### 3.3 Mesh Placement

`mesh.templates` supports anchor placement relative to the simulation box.

```toml
[[mesh.templates]]
kind = "plane"
size_mode = "box_fraction"
size_frac = [1.0, 1.0]
placement_mode = "box_anchor"
anchor = "z_low_face_center"
offset_frac = [0.0, 0.0, 0.02]
nx = 20
ny = 20
```

The loader resolves this to `size_x` / `size_y` / `center` and other final template keys.

### 3.4 Group Placement

`mesh.groups` is a table for shared origins and scales across multiple templates.

```toml
[mesh.groups.cavity_unit]
placement_mode = "box_anchor"
anchor = "box_center"
scale_from = "box_x"
scale_factor = 0.5

[[mesh.templates]]
group = "cavity_unit"
kind = "sphere"
radius = 0.2
center = [0.0, 0.0, 0.0]
```

While loading, `mesh.groups`, `group`, `scale_from`, and related keys determine each template's
physical coordinates and dimensions.

## 4. Schema

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

## 5. Common Mistakes

### 5.1 Top-Level Key Placement

Runtime settings belong under `sim`, `particles`, `mesh`, and `output`.
Ordinary keys before the first section, or unknown top-level sections, fail validation or Fortran loading.

### 5.2 Mixing High-Level and Runtime Keys

Using high-level keys and their runtime equivalents for the same value, such as `box_origin` / `box_size`
with `box_min` / `box_max`, fails validation. Use one style for a given quantity.

### 5.3 Checking Before Running

Run `beachx lint beach.toml` to check TOML parsing, JSON Schema, high-level notation consistency,
and known BEACH constraints together.
