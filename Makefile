.PHONY: \
	install install-local install-auto install-generic install-camphor install-camphor-local \
	install-intel install-intel-local \
	build check run \
	static-check schema-check \
	test test-l0 test-l1 test-l2 test-l3 test-heavy test-full test-physics-release \
	test-fortran test-fortran-light test-fortran-contract test-fortran-heavy test-fortran-release-correctness \
	test-fortran-far-correction test-fortran-far-correction-diagnostics test-fortran-benchmark \
	test-field-kernel-cache \
	test-python test-quick test-ci test-local \
	build-kernel \
	fmt-fortran fmt-check-fortran install-hooks \
	build-mpi run-mpi test-mpi test-mpi-periodic-cache \
	docs-deps docs-starlight docs-fortran docs-pages docs-clean

.NOTPARALLEL: test test-l0 test-l1 test-l2 test-l3 test-heavy test-full test-quick test-ci test-local \
	test-fortran test-fortran-light test-fortran-contract test-fortran-heavy test-fortran-release-correctness \
	test-fortran-far-correction test-fortran-far-correction-diagnostics test-fortran-benchmark \
	test-field-kernel-cache \
	test-mpi test-mpi-periodic-cache

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
INTEL_TEST_QUIET_FLAG = $(if $(filter ifort ifx mpiifort mpiifx,$(notdir $(FPM_FC))),-check noarg_temp_created)
FORTRAN_TEST_FLAGS ?= $(strip $(OPENMP_FLAG) $(INTEL_TEST_QUIET_FLAG))
VERSION_MODE ?=
BUILD_VERSION_MODE ?= $(if $(VERSION_MODE),$(VERSION_MODE),git)
CHECK_VERSION_MODE ?= $(if $(VERSION_MODE),$(VERSION_MODE),dev)
RUN_VERSION_MODE ?= $(if $(VERSION_MODE),$(VERSION_MODE),dev)
SCHEMAS ?= schemas/beach.schema.json
PHYSICS_RELEASE_MANIFEST ?= build/physics-release/manifest.txt
FORTRAN_L1_TARGETS ?= \
	test_version \
	test_app_config_parser \
	test_physics_config_types \
	test_charge_ledger \
	test_model_fingerprint \
	test_panel_moments \
	test_panel_kernel \
	test_panel_geometry_near \
	test_periodic2_cache_codec \
	test_periodic_zero_mode \
	test_outer_plasma_linear \
	test_outer_plasma_kinetic_core \
	test_outer_plasma_grid \
	test_outer_plasma_kinetic_runtime \
	test_outer_plasma_local_mean \
	test_periodic2_nonzero_tail \
	test_outer_plasma_unified \
	test_electrostatic_unified \
	test_electrostatic_snapshot \
	test_outer_coupler \
	test_outer_plasma_interface \
	test_outer_plasma_orbit \
	test_outer_plasma_photoelectron \
	test_interface_particle_buffer \
	test_boundary \
	test_restart \
	test_reservoir_injection \
	test_external_field_velocity_grid \
	test_dynamics_basic \
	test_particle_stepper \
	test_dynamics_field_solver \
	test_surface_models \
	test_templates_importers_runtime \
	test_simulator \
	test_injection_sampling \
	test_sheath_injection_model \
	test_sheath_model_core \
	test_performance_profile \
	test_output_writer_io \
	test_output_writer_potential
FORTRAN_L2_TARGETS ?= \
	test_field_kernel_c \
	test_periodic_zero_mode_c
FORTRAN_L3_TARGETS ?= \
	test_dynamics_fmm \
	test_coulomb_fmm_core_basic \
	test_panel_near_correction \
	test_coulomb_fmm_core_panel \
	test_dynamics_panel_fmm \
	test_outer_plasma_kinetic
FORTRAN_RELEASE_CONVERGENCE_TARGETS ?= \
	test_dynamics_basic \
	test_electrostatic_unified \
	test_outer_plasma_orbit
FORTRAN_FAR_CORRECTION_TARGETS ?= \
	test_field_kernel_cache_c \
	test_periodic2_operator_cache \
	test_periodic2_infinite_operator \
	test_periodic2_cached_snapshot \
	test_coulomb_fmm_core_periodic \
	test_periodic_nonzero_reference
FORTRAN_FAR_CORRECTION_DIAGNOSTIC_TARGETS ?= \
	test_periodic2_flat_oracle_diag
FORTRAN_BENCHMARK_TARGETS ?= \
	benchmark_periodic2_runtime
KERNEL_FC ?= gfortran
KERNEL_LIB ?= build/libbeach_field_kernel.so
KERNEL_FPM_FLAG ?= $(OPENMP_FLAG) -fPIC
INSTALL_PROFILE ?= auto
FPRETTIFY ?= fprettify
PRE_COMMIT ?= pre-commit
DOCS_PROJECT_FILE ?= ford.md
DOCS_OUTPUT_DIR ?= build/ford-docs
DOCS_SITE_DIR ?= docs-site
DOCS_SITE_OUTPUT_DIR ?= build/starlight-site
DOCS_SITE_FORTRAN_DIR ?= $(DOCS_SITE_OUTPUT_DIR)/fortran
FORTRAN_DEP_MAP_MD ?= docs/FortranDependencyMap.md
FORTRAN_DEP_MAP_EN_MD ?= docs/FortranDependencyMap.en.md
FORTRAN_DEP_MAP_DOT ?= build/fortran_module_dependencies.dot
FORTRAN_DEP_MAP_SVG ?= docs/media/fortran_module_dependencies.svg

MPI_FC ?= mpiifort
MPI_OPENMP_FLAG ?= -qopenmp
MPI_CPP_FLAG ?= -fpp -DUSE_MPI
MPI_TEST_QUIET_FLAG = $(if $(filter ifort ifx mpiifort mpiifx,$(notdir $(MPI_FC))),-check noarg_temp_created)
MPI_TEST_FLAGS ?= $(strip $(MPI_CPP_FLAG) $(MPI_OPENMP_FLAG) $(MPI_TEST_QUIET_FLAG))
MPI_NP ?= 2
MPI_RUNNER ?= mpirun -n $(MPI_NP)

define run_fortran_targets
	@BEACH_VERSION_MODE=$(CHECK_VERSION_MODE) FPM=$(FPM) FPM_ACTION=$(3) FPM_TARGET_KIND=$(4) \
		FPM_PROFILE=$(2) \
		FPM_FFLAGS="$(if $(filter release,$(2)),$(OPENMP_FLAG),$(FORTRAN_TEST_FLAGS))" \
		BUILD_SH="$(BUILD_SH)" \
		tools/run_fortran_targets.sh $(1)
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
		FPM_PROFILE=$(PROFILE) FPM_FC="$(KERNEL_FC)" FPM_FFLAGS="$(KERNEL_FPM_FLAG)" $(BUILD_SH)
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

test-physics-release:
	tools/run_physics_release_gate.sh --manifest "$(PHYSICS_RELEASE_MANIFEST)"

test-full:
	BEACH_VERSION_MODE=$(CHECK_VERSION_MODE) FPM=$(FPM) FPM_ACTION=test \
		FPM_PROFILE=debug FPM_FFLAGS="$(FORTRAN_TEST_FLAGS)" $(BUILD_SH)

test-fortran: test-fortran-light

test-fortran-light:
	$(call run_fortran_targets,$(FORTRAN_L1_TARGETS),debug,test)

test-fortran-contract:
	$(call run_fortran_targets,$(FORTRAN_L2_TARGETS),debug,test)

test-fortran-heavy:
	$(call run_fortran_targets,$(FORTRAN_L3_TARGETS),debug,test)

test-fortran-release-correctness:
	$(call run_fortran_targets,$(FORTRAN_RELEASE_CONVERGENCE_TARGETS) $(FORTRAN_L3_TARGETS),debug,test)

test-fortran-far-correction:
	$(call run_fortran_targets,$(FORTRAN_FAR_CORRECTION_TARGETS),debug,test)

test-fortran-far-correction-diagnostics:
	$(call run_fortran_targets,$(FORTRAN_FAR_CORRECTION_DIAGNOSTIC_TARGETS),debug,test)

test-fortran-benchmark:
	$(call run_fortran_targets,$(FORTRAN_BENCHMARK_TARGETS),release,run,example)

test-field-kernel-cache: export BEACH_RUN_FIELD_KERNEL_CACHE_TESTS=1
test-field-kernel-cache: build-kernel
	$(call run_fortran_targets,test_field_kernel_cache_c,debug,test)
	$(PYTHON) -m pytest -q tests/python/test_field_kernel_cache.py
	BEACH_FIELD_KERNEL_LIB="$(abspath $(KERNEL_LIB))" \
		$(PYTHON) -m pytest -q \
		tests/python/test_periodic_force_oracle.py::test_validation_tool_native_plane_oracles_match_receipt_contract

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
		FPM_PROFILE=debug FPM_FC=$(MPI_FC) FPM_FFLAGS="$(MPI_TEST_FLAGS)" \
		$(BUILD_SH) --target test_mpi_hybrid --runner "$(MPI_RUNNER)"

test-mpi-periodic-cache:
	BEACH_VERSION_MODE=$(CHECK_VERSION_MODE) FPM=$(FPM) FPM_ACTION=test \
		FPM_PROFILE=debug FPM_FC=$(MPI_FC) FPM_FFLAGS="$(MPI_TEST_FLAGS)" \
		$(BUILD_SH) --target test_periodic2_operator_cache_mpi --runner "$(MPI_RUNNER)"

docs-deps:
	npm --prefix $(DOCS_SITE_DIR) ci

docs-starlight:
	$(PYTHON) tools/generate_fortran_dependency_report.py \
		--markdown $(FORTRAN_DEP_MAP_MD) \
		--markdown-en $(FORTRAN_DEP_MAP_EN_MD) \
		--dot $(FORTRAN_DEP_MAP_DOT) \
		--svg $(FORTRAN_DEP_MAP_SVG)
	$(PYTHON) tools/sync_starlight_docs.py
	npm --prefix $(DOCS_SITE_DIR) run build

docs-fortran:
	$(PYTHON) tools/generate_fortran_dependency_report.py \
		--markdown $(FORTRAN_DEP_MAP_MD) \
		--markdown-en $(FORTRAN_DEP_MAP_EN_MD) \
		--dot $(FORTRAN_DEP_MAP_DOT) \
		--svg $(FORTRAN_DEP_MAP_SVG)
	$(FORD) $(DOCS_PROJECT_FILE) --output_dir $(DOCS_OUTPUT_DIR) --config '$(FORD_CONFIG)'

docs-pages: docs-fortran
	$(PYTHON) tools/sync_starlight_docs.py
	npm --prefix $(DOCS_SITE_DIR) run build
	rm -rf $(DOCS_SITE_FORTRAN_DIR)
	mkdir -p $(DOCS_SITE_FORTRAN_DIR)
	cp -a $(DOCS_OUTPUT_DIR)/. $(DOCS_SITE_FORTRAN_DIR)/

docs-clean:
	rm -rf $(DOCS_OUTPUT_DIR) $(DOCS_SITE_OUTPUT_DIR) \
		$(DOCS_SITE_DIR)/src/content/docs \
		$(DOCS_SITE_DIR)/public/images \
		$(DOCS_SITE_DIR)/public/media
