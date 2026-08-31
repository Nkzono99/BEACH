title: Installation

Lang: [English](Installation.en.md) | [日本語](Installation.md)

# Installation

Installing the BEACH Python package also installs `beachx` and builds the Fortran `beach` executable. A normal
`pip install` automatically installs the Fortran Package Manager (`fpm`) into the isolated build environment.

## Requirements

| Item | Requirement |
| --- | --- |
| OS | Linux is the primary tested environment; use site compiler/MPI modules on HPC systems |
| Python | 3.10 or newer |
| Source retrieval | `git` |
| Build tool | `make` |
| Fortran | `gfortran`, Intel Fortran, or another fpm-compatible compiler |

```bash
python --version
git --version
make --version
gfortran --version
```

## Install the version described by this site

This site documents the configuration and output behavior on the GitHub `main` branch. Install the current source to
reproduce the official 20-batch tutorial and its expected values.

```bash
python -m pip install -U pip setuptools wheel
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"
beach --version
beachx --help
```

If the commands are not found, inspect the user base and add its `bin` directory to `PATH`:

```bash
python -m site --user-base
export PATH="$HOME/.local/bin:$PATH"
```

## Stable release from PyPI

```bash
python -m pip install beach-bem
```

The PyPI package is a stable snapshot. Its `beachx config init` tutorial may lag behind the `main` branch documented by
this site. Do not combine an older generated case with the 20-batch expected values; use the GitHub installation above
when following this site.

## Work from an editable checkout

```bash
python -m pip install -e . --no-build-isolation
```

`--no-build-isolation` disables automatic installation of the build dependencies declared in
`pyproject.toml`. Install `fpm` on `PATH` before using this mode or invoking `make` directly.

```bash
python -m pip install fpm
fpm --version
```

Current source builds also provide the `beach-zhao-response` helper used to create matching-plane response tables.
PyPI distributions may not include that auxiliary executable. Install the GitHub version above when you
need it, then continue with the
[matching-plane numerical and response-table reference](MatchingPlaneReference.en.html#build-a-response-table-for-the-table-backend).

## Upgrade or remove

```bash
python -m pip install -U beach-bem
python -m pip uninstall beach-bem
```

See [Troubleshooting](Troubleshooting.en.html) for compiler, fpm, PATH, and runtime failures.
Continue with the [Ten-Minute Tutorial](Tutorial.en.html).
