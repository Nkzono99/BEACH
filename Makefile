.PHONY: \
	install install-local install-auto install-generic install-camphor install-camphor-local \
	install-intel install-intel-local \
	build check run \
	static-check schema-check \
	test test-l0 test-l1 test-l2 test-l3 test-heavy test-full \
	test-fortran test-fortran-light test-fortran-contract test-fortran-heavy \
	test-python test-quick test-ci test-local \
	build-kernel \
	fmt-fortran fmt-check-fortran install-hooks \
	build-mpi run-mpi test-mpi \
	docs-fortran docs-clean

.NOTPARALLEL: test test-l0 test-l1 test-l2 test-l3 test-heavy test-full test-quick test-ci test-local \
	test-fortran test-fortran-light test-fortran-contract test-fortran-heavy test-mpi

.DEFAULT_GOAL := install

PYTHON ?= $(shell if command -v python >/dev/null 2>&1; then echo python; \
	elif command -v python3.11 >/dev/null 2>&1; then echo python3.11; \
	elif command -v python3.12 >/dev/null 2>&1; then echo python3.12; \
	else echo python3; fi)
FPM ?= fpm
BUILD_SH ?= ./build.sh
FORD ?= $(PYTHON) -m ford
FORD_CONFIG ?= preprocessor = "$(PYTHON) -W ignore::RuntimeWarning -m pcpp.pcmd -D__GFORTRAN__ --passthru-comments"
PROFILE ?= release
CONFIG ?= beach.toml
OPENMP_FLAG ?= -fopenmp
VERSION_MODE ?=
BUILD_VERSION_MODE ?= $(if $(VERSION_MODE),$(VERSION_MODE),git)
CHECK_VERSION_MODE ?= $(if $(VERSION_MODE),$(VERSION_MODE),dev)
RUN_VERSION_MODE ?= $(if $(VERSION_MODE),$(VERSION_MODE),dev)
SCHEMAS ?= schemas/beach.schema.json schemas/beach.case.schema.json schemas/beach.preset.schema.json
FORTRAN_L1_TARGETS ?= \
	test_version \
	test_app_config_parser \
	test_boundary \
	test_restart \
	test_reservoir_injection \
	test_external_field_velocity_grid \
	test_dynamics_basic \
	test_dynamics_field_solver \
	test_templates_importers_runtime \
	test_simulator \
	test_injection_sampling \
	test_sheath_injection_model \
	test_sheath_model_core \
	test_performance_profile \
	test_output_writer_io \
	test_output_writer_potential
FORTRAN_L2_TARGETS ?= \
	test_field_kernel_c
FORTRAN_L3_TARGETS ?= \
	test_dynamics_fmm \
	test_coulomb_fmm_core_basic \
	test_coulomb_fmm_core_periodic \
	test_periodic2_flat_oracle_diag
KERNEL_FC ?= gfortran
KERNEL_LIB ?= build/libbeach_field_kernel.so
KERNEL_FPM_FLAG ?= $(OPENMP_FLAG) -fPIC
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

define run_fortran_targets
	@set -eu; \
	for target in $(1); do \
		echo "==> fpm test --target $$target"; \
		BEACH_VERSION_MODE=$(CHECK_VERSION_MODE) FPM=$(FPM) FPM_ACTION=test \
			FPM_PROFILE=debug FPM_FFLAGS="$(OPENMP_FLAG)" $(BUILD_SH) --target "$$target"; \
	done
endef

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
	BEACH_VERSION_MODE=$(BUILD_VERSION_MODE) FPM=$(FPM) FPM_ACTION=build \
		FPM_PROFILE=$(PROFILE) FPM_FFLAGS="$(OPENMP_FLAG)" $(BUILD_SH)

check:
	BEACH_VERSION_MODE=$(CHECK_VERSION_MODE) FPM=$(FPM) FPM_ACTION=build \
		FPM_PROFILE=debug FPM_FFLAGS="$(OPENMP_FLAG)" $(BUILD_SH)

build-kernel:
	BEACH_VERSION_MODE=$(BUILD_VERSION_MODE) FPM=$(FPM) FPM_ACTION=build \
		FPM_PROFILE=$(PROFILE) FPM_FFLAGS="$(KERNEL_FPM_FLAG)" $(BUILD_SH)
	@set -eu; \
	lib=$$(find build -name libbeach_fortran.a -printf '%T@ %p\n' | sort -nr | awk 'NR==1 {print $$2}'); \
	if [ -z "$$lib" ]; then echo "libbeach_fortran.a not found; run fpm build first." >&2; exit 1; fi; \
	mkdir -p "$$(dirname "$(KERNEL_LIB)")"; \
	$(KERNEL_FC) -shared -o "$(KERNEL_LIB)" \
		-Wl,--whole-archive "$$lib" -Wl,--no-whole-archive $(OPENMP_FLAG); \
	echo "built $(KERNEL_LIB)"

run:
	BEACH_VERSION_MODE=$(RUN_VERSION_MODE) FPM=$(FPM) FPM_ACTION=run \
		FPM_PROFILE=$(PROFILE) FPM_FFLAGS="$(OPENMP_FLAG)" $(BUILD_SH) -- $(CONFIG)

static-check:
	git diff --check

schema-check:
	@set -eu; \
	for schema in $(SCHEMAS); do \
		echo "==> validate $$schema"; \
		$(PYTHON) -m json.tool "$$schema" >/dev/null; \
	done

test-l0: static-check schema-check check

test: test-l1

test-l1: test-l0 test-python test-fortran-light

test-l2: test-l1 test-fortran-contract

test-l3: test-l2 test-fortran-heavy

test-heavy: test-fortran-heavy

test-full:
	BEACH_VERSION_MODE=$(CHECK_VERSION_MODE) FPM=$(FPM) FPM_ACTION=test \
		FPM_PROFILE=debug FPM_FFLAGS="$(OPENMP_FLAG)" $(BUILD_SH)

test-fortran: test-fortran-light

test-fortran-light:
	$(call run_fortran_targets,$(FORTRAN_L1_TARGETS))

test-fortran-contract:
	$(call run_fortran_targets,$(FORTRAN_L2_TARGETS))

test-fortran-heavy:
	$(call run_fortran_targets,$(FORTRAN_L3_TARGETS))

test-python:
	$(PYTHON) -m pytest -q

test-quick: test-l1

test-ci: test-l2

test-local: test-l1

fmt-fortran:
	find src app tests/fortran -type f \( -name '*.f90' -o -name '*.F90' \) -exec $(FPRETTIFY) -i 2 {} +

fmt-check-fortran:
	find src app tests/fortran -type f \( -name '*.f90' -o -name '*.F90' \) -exec $(FPRETTIFY) -i 2 --silent {} +
	git diff --exit-code -- src app tests/fortran

install-hooks:
	$(PRE_COMMIT) install

build-mpi:
	BEACH_VERSION_MODE=$(BUILD_VERSION_MODE) FPM=$(FPM) FPM_ACTION=build \
		FPM_PROFILE=$(PROFILE) FPM_FC=$(MPI_FC) FPM_FFLAGS="$(MPI_CPP_FLAG) $(MPI_OPENMP_FLAG)" $(BUILD_SH)

run-mpi:
	BEACH_VERSION_MODE=$(RUN_VERSION_MODE) FPM=$(FPM) FPM_ACTION=run \
		FPM_PROFILE=$(PROFILE) FPM_FC=$(MPI_FC) FPM_FFLAGS="$(MPI_CPP_FLAG) $(MPI_OPENMP_FLAG)" \
		$(BUILD_SH) --runner "$(MPI_RUNNER)" -- $(CONFIG)

test-mpi:
	BEACH_VERSION_MODE=$(CHECK_VERSION_MODE) FPM=$(FPM) FPM_ACTION=test \
		FPM_PROFILE=debug FPM_FC=$(MPI_FC) FPM_FFLAGS="$(MPI_CPP_FLAG) $(MPI_OPENMP_FLAG)" \
		$(BUILD_SH) --target test_mpi_hybrid --runner "$(MPI_RUNNER)"

docs-fortran:
	$(PYTHON) tools/generate_fortran_dependency_report.py \
		--markdown $(FORTRAN_DEP_MAP_MD) \
		--dot $(FORTRAN_DEP_MAP_DOT) \
		--svg $(FORTRAN_DEP_MAP_SVG)
	$(FORD) $(DOCS_PROJECT_FILE) --output_dir $(DOCS_OUTPUT_DIR) --config '$(FORD_CONFIG)'

docs-clean:
	rm -rf $(DOCS_OUTPUT_DIR)
