title: Edit configuration

Lang: [English](Configuration.en.md) | [日本語](Configuration.md)

# Edit configuration

This document describes the directly edited `beach.toml` file and the `beachx config` helper commands.

- The Fortran runtime `beach` reads `beach.toml` directly.
- `beachx config init` creates a small runnable `beach.toml`.
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
It is the nonperiodic official beginner case that launches one `volume_seed` electron toward an insulating plane with
`field_solver="direct"` and `field_bc_mode="free"`. It does not include FMM, periodic boundaries, an ion species, or
`photo_raycast`.

### 2.2 `lint`

Run TOML parsing, JSON Schema validation, coordinate and placement combination checks, and BEACH-specific constraints together.

```bash
beachx lint beach.toml
beachx lint run.toml --schema schemas/beach.schema.json
```

### 2.3 `validate`

Read `beach.toml` and validate coordinate and placement combinations plus known constraints of the final settings.
Use `beachx lint` when JSON Schema validation should also be included.

```bash
beachx config validate
beachx config validate run.toml
```

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

Runtime settings belong under `sim`, `particles`, `mesh`, and `output`.
Ordinary keys before the first section, or unknown top-level sections, fail validation or Fortran loading.

### 4.2 Specifying the Same Coordinate Two Ways

Writing the same coordinate in two forms, such as `box_origin` / `box_size` together with `box_min` / `box_max`, fails
validation. `size_mode="box_fraction"` and group scaling intentionally replace corresponding dimensions with calculated values;
the affected keys are listed in [Input Parameters Reference](Parameters.en.html#coordinate-and-placement-helper-parameters).

### 4.3 Checking Before Running

Run `beachx lint beach.toml` to check TOML parsing, JSON Schema, coordinate and placement combinations, and known BEACH
constraints together.
