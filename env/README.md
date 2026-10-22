# Build Profiles

`install.sh` can load host-specific build presets from this directory.

- `common.env`: shared defaults
- `camphor.env`: optimized Intel MPI + OpenMP flags for camphor (`-fp-model precise`, bit-reproducible against the test builds)
- `camphor-fast.env`: opt-in camphor build with `-fp-model fast=1`; about 1.3x faster per batch because the near-field
  kernel can use 4-lane SVML `log`/`atan2`. Results differ from `camphor.env` at the 1e-13 level and are not
  bit-reproducible against it, so use it for production runs only and keep `camphor.env` for tests and golden
  comparisons. It sets `FPM_NO_PROFILE=1` because fpm would otherwise append `-fp-model precise` after `FFLAGS`.
- `profile.template.env`: new profile template

Profile selection:

- `BUILD_PROFILE=auto` (default): detect by hostname (`camphor* -> camphor`, `fn* -> fugaku`, otherwise `generic`)
- `BUILD_PROFILE=generic`: portable fallback build
- `BUILD_PROFILE=<name>`: load `env/<name>.env`

Variables in each `*.env` file should use `: "${VAR:=...}"` style so user-provided env vars can still override.

`FPM_NO_PROFILE=1` makes `build.sh` pass only `FFLAGS` to fpm (no `--profile`). fpm appends its profile flags
after `--flag`, so without it a profile's `-O3 -fp-model precise` overrides any `-fp-model` given in `FFLAGS`;
a profile that disables fpm's flags must spell out `-O3` and the rest itself (see `camphor-fast.env`).

To add a new host profile:

1. `cp env/profile.template.env env/<profile>.env`
2. Fill `FC`, `FFLAGS`, module settings (`USE_MODULES`, `MODULES_TO_LOAD`, `MODULES_TO_UNLOAD`)
3. Run with `BUILD_PROFILE=<profile> ./install.sh` or `make install INSTALL_PROFILE=<profile>`
