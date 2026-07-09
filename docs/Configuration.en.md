title: beachx config / High-Level Notation Guide

Lang: [English](Configuration.en.md) | [日本語](Configuration.md)

# `beachx config` / High-Level Notation Guide

This document describes the directly edited `beach.toml` file and the `beachx config` helper commands.

- The Fortran runtime `beach` reads a TOML file expanded to final runtime keys.
- `beachx config init` creates a small runnable `beach.toml`.
- `beachx config render` expands high-level notation in `beach.toml` into final runtime keys.
- See [Input Parameters Reference](Parameters.en.html) for the final keys read by the Fortran runtime.

## 1. Basic Flow

```bash
mkdir run_periodic2
cd run_periodic2

beachx config init
$EDITOR beach.toml
beachx lint beach.toml
beachx config render beach.toml --output beach.rendered.toml
beach beach.rendered.toml
```

When no output path is specified, `render` overwrites the input file. For first-time use and shared examples,
prefer `--output` so the source config and rendered runtime config remain separate. Use `--stdout` to inspect
the rendered result first.

```bash
beachx config render --stdout
beachx config render beach.toml --output beach.rendered.toml
```

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

### 2.4 `render`

Expand high-level notation into final keys consumed by the Fortran runtime.

```bash
beachx config render
beachx config render run.toml --output run.rendered.toml
```

### 2.5 `diff`

Semantically compare two configs. By default, high-level notation is rendered before comparison.

```bash
beachx config diff left.toml right.toml
beachx config diff --raw left.toml right.toml
```

## 3. High-Level Notation

High-level notation is a convenience layer for writing coordinate intent directly in `beach.toml`.
After `beachx config render`, it is converted to ordinary numerical keys.

### 3.1 Simulation Box

Use `sim.box_origin` and `sim.box_size` to derive `box_min` / `box_max`.

```toml
[sim]
use_box = true
box_origin = [0.0, 0.0, -1.0]
box_size = [1.0, 1.0, 2.0]
```

After rendering:

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

Rendering expands these values into `pos_low` / `pos_high`.

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

Rendering expands this to `size_x` / `size_y` / `center` and other final template keys.

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

After rendering, `mesh.groups`, `group`, `scale_from`, and other high-level keys are removed and each
template is written with physical coordinates and dimensions.

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

### 5.2 High-Level Keys Disappear After Rendering

This is expected. Keys such as `box_origin`, `box_size`, `inject_region_mode`, and `mesh.groups` are not final
runtime keys read by `beach`.

### 5.3 Avoiding In-Place Overwrite

Pass `--stdout` or `--output` to `render`.

```bash
beachx config render --stdout
beachx config render beach.toml --output beach.rendered.toml
```
