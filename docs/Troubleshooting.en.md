title: Troubleshooting

Lang: [English](Troubleshooting.en.md) | [日本語](Troubleshooting.md)

# Troubleshooting

| Symptom | Check |
| --- | --- |
| `gfortran` or `fpm` missing | Verify prerequisites; load the compiler module on HPC systems |
| `beach` or `beachx` missing | Inspect `python -m site --user-base` and add its `bin` directory to `PATH` |
| lint passes but runtime stops | Read the model-combination error; Fortran physics validation can be stricter than schema validation |
| no output directory | Check `write_files=true`, `output.dir`, and process exit status |
| many `survived_max_step` particles | Revisit `dt`, `max_step`, box size, and injection velocity; do not relabel unresolved particles |
| empty or huge history | Check `history_stride`; increase it and enable potential history only when needed |
| restart rejected | Compare model/mesh/species fingerprints, `restart_from`, and cumulative `batch_count` |
| unsupported conductor/periodic combination or dielectric input | Consult the support matrix; dielectric polarization is not implemented |
| slow far correction | Reuse a `cached_kneq0` warm cache and check how often cold operator generation occurs |

For a reproducible issue, provide the configuration, BEACH version, compiler/MPI, rank/thread counts,
complete error, and a minimal mesh. Do not loosen numerical tolerances to bypass a physically inapplicable model.
