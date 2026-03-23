.PHONY: \
	install install-local install-auto install-generic install-camphor install-camphor-local \
	install-intel install-intel-local \
	build run test \
	fmt-fortran fmt-check-fortran install-hooks \
	build-mpi run-mpi test-mpi \
	docs-fortran docs-clean

PYTHON ?= $(shell if command -v python >/dev/null 2>&1; then echo python; \
	elif command -v python3.11 >/dev/null 2>&1; then echo python3.11; \
	elif command -v python3.12 >/dev/null 2>&1; then echo python3.12; \
	else echo python3; fi)
FPM ?= fpm
FORD ?= $(PYTHON) -m ford
FORD_CONFIG ?= preprocessor = "$(PYTHON) -W ignore::RuntimeWarning -m pcpp.pcmd -D__GFORTRAN__ --passthru-comments"
PROFILE ?= release
CONFIG ?= beach.toml
OPENMP_FLAG ?= -fopenmp
INSTALL_PROFILE ?= auto
FPRETTIFY ?= fprettify
PRE_COMMIT ?= pre-commit
DOCS_PROJECT_FILE ?= ford.md
DOCS_OUTPUT_DIR ?= build/ford-docs
FORTRAN_DEP_MAP_MD ?= docs/fortran_dependency_map.md
FORTRAN_DEP_MAP_DOT ?= build/fortran_module_dependencies.dot
FORTRAN_DEP_MAP_SVG ?= docs/media/fortran_module_dependencies.svg

MPI_FC ?= mpiifort
MPI_OPENMP_FLAG ?= -qopenmp
MPI_CPP_FLAG ?= -fpp -DUSE_MPI
MPI_NP ?= 2
MPI_RUNNER ?= mpirun -n $(MPI_NP)

install:
	BUILD_PROFILE=$(INSTALL_PROFILE) ./install.sh

install-local:
	BUILD_PROFILE=$(INSTALL_PROFILE) PREFIX=$(PWD)/.local ./install.sh

install-auto:
	BUILD_PROFILE=auto ./install.sh

install-generic:
	BUILD_PROFILE=generic ./install.sh

install-camphor:
	BUILD_PROFILE=camphor ./install.sh

install-camphor-local:
	BUILD_PROFILE=camphor PREFIX=$(PWD)/.local-intel ./install.sh

install-intel:
	$(MAKE) install-camphor

install-intel-local:
	$(MAKE) install-camphor-local

build:
	$(FPM) build --profile $(PROFILE) --flag "$(OPENMP_FLAG)"

run:
	$(FPM) run --profile $(PROFILE) --flag "$(OPENMP_FLAG)" -- $(CONFIG)

test:
	$(FPM) test --profile debug --flag "$(OPENMP_FLAG)"

fmt-fortran:
	find src app tests/fortran -type f \( -name '*.f90' -o -name '*.F90' \) -exec $(FPRETTIFY) -i 2 {} +

fmt-check-fortran:
	find src app tests/fortran -type f \( -name '*.f90' -o -name '*.F90' \) -exec $(FPRETTIFY) -i 2 --silent {} +
	git diff --exit-code -- src app tests/fortran

install-hooks:
	$(PRE_COMMIT) install

build-mpi:
	FPM_FC=$(MPI_FC) $(FPM) build --profile $(PROFILE) --flag "$(MPI_CPP_FLAG) $(MPI_OPENMP_FLAG)"

run-mpi:
	FPM_FC=$(MPI_FC) $(FPM) run --profile $(PROFILE) --flag "$(MPI_CPP_FLAG) $(MPI_OPENMP_FLAG)" \
		--runner "$(MPI_RUNNER)" -- $(CONFIG)

test-mpi:
	FPM_FC=$(MPI_FC) $(FPM) test --target test_mpi_hybrid --profile debug \
		--flag "$(MPI_CPP_FLAG) $(MPI_OPENMP_FLAG)" --runner "$(MPI_RUNNER)"

docs-fortran:
	$(PYTHON) tools/generate_fortran_dependency_report.py \
		--markdown $(FORTRAN_DEP_MAP_MD) \
		--dot $(FORTRAN_DEP_MAP_DOT) \
		--svg $(FORTRAN_DEP_MAP_SVG)
	$(FORD) $(DOCS_PROJECT_FILE) --output_dir $(DOCS_OUTPUT_DIR) --config '$(FORD_CONFIG)'

docs-clean:
	rm -rf $(DOCS_OUTPUT_DIR)
