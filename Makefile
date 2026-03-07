.PHONY: \
	install install-local install-auto install-generic install-camphor install-camphor-local \
	install-intel install-intel-local \
	build run test \
	build-mpi run-mpi test-mpi

FPM ?= fpm
PROFILE ?= release
CONFIG ?= beach.toml
OPENMP_FLAG ?= -fopenmp
INSTALL_PROFILE ?= auto

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

build-mpi:
	FPM_FC=$(MPI_FC) $(FPM) build --profile $(PROFILE) --flag "$(MPI_CPP_FLAG) $(MPI_OPENMP_FLAG)"

run-mpi:
	FPM_FC=$(MPI_FC) $(FPM) run --profile $(PROFILE) --flag "$(MPI_CPP_FLAG) $(MPI_OPENMP_FLAG)" \
		--runner "$(MPI_RUNNER)" -- $(CONFIG)

test-mpi:
	FPM_FC=$(MPI_FC) $(FPM) test --target test_mpi_hybrid --profile debug \
		--flag "$(MPI_CPP_FLAG) $(MPI_OPENMP_FLAG)" --runner "$(MPI_RUNNER)"
