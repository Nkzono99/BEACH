.PHONY: \
	install-intel install-intel-local \
	build run test \
	build-mpi run-mpi test-mpi

FPM ?= fpm
PROFILE ?= release
CONFIG ?= beach.toml
OPENMP_FLAG ?= -fopenmp

MPI_FC ?= mpiifort
MPI_OPENMP_FLAG ?= -qopenmp
MPI_CPP_FLAG ?= -fpp -DUSE_MPI
MPI_NP ?= 2
MPI_RUNNER ?= mpirun -n $(MPI_NP)

install-intel:
	./install.sh

install-intel-local:
	PREFIX=$(PWD)/.local-intel ./install.sh

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
