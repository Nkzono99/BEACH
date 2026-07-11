title: Installation

Lang: [English](Installation.en.md) | [日本語](Installation.md)

# Installation

Installing the BEACH Python package also builds the Fortran `beach` executable.

## Requirements

| Item | Requirement |
| --- | --- |
| OS | Linux is the primary tested environment; use site compiler/MPI modules on HPC systems |
| Python | 3.10 or newer |
| Build tool | `make` |
| Fortran | `gfortran`, Intel Fortran, or another fpm-compatible compiler |
| Fortran package manager | `fpm` |

```bash
python --version
make --version
gfortran --version
fpm --version
```

## Install from PyPI

```bash
python -m pip install -U pip setuptools wheel
python -m pip install beach-bem
beach --version
beachx --help
```

If the commands are not found, inspect the user base and add its `bin` directory to `PATH`:

```bash
python -m site --user-base
export PATH="$HOME/.local/bin:$PATH"
```

## Development version

```bash
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"
python -m pip install -e . --no-build-isolation  # from a checkout
```

## Upgrade or remove

```bash
python -m pip install -U beach-bem
python -m pip uninstall beach-bem
```

See [Troubleshooting](Troubleshooting.en.html) for compiler, fpm, PATH, and runtime failures.
Continue with the [Ten-Minute Tutorial](Tutorial.en.html).
