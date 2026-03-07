# Build Profiles

`install.sh` can load host-specific build presets from this directory.

- `common.env`: shared defaults
- `camphor.env`: optimized Intel MPI + OpenMP flags for camphor

Profile selection:

- `BUILD_PROFILE=auto` (default): detect by hostname (`camphor* -> camphor`, `fn* -> fugaku`, otherwise `generic`)
- `BUILD_PROFILE=generic`: portable fallback build
- `BUILD_PROFILE=<name>`: load `env/<name>.env`

Variables in each `*.env` file should use `: "${VAR:=...}"` style so user-provided env vars can still override.
