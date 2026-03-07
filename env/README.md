# Build Profiles

`install.sh` can load host-specific build presets from this directory.

- `common.env`: shared defaults
- `camphor.env`: optimized Intel MPI + OpenMP flags for camphor
- `profile.template.env`: new profile template

Profile selection:

- `BUILD_PROFILE=auto` (default): detect by hostname (`camphor* -> camphor`, `fn* -> fugaku`, otherwise `generic`)
- `BUILD_PROFILE=generic`: portable fallback build
- `BUILD_PROFILE=<name>`: load `env/<name>.env`

Variables in each `*.env` file should use `: "${VAR:=...}"` style so user-provided env vars can still override.

To add a new host profile:

1. `cp env/profile.template.env env/<profile>.env`
2. Fill `FC`, `FFLAGS`, module settings (`USE_MODULES`, `MODULES_TO_LOAD`, `MODULES_TO_UNLOAD`)
3. Run with `BUILD_PROFILE=<profile> ./install.sh` or `make install INSTALL_PROFILE=<profile>`
