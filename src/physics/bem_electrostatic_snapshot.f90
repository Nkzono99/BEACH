module bem_electrostatic_snapshot
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: eps0
  use bem_types, only: mesh_type, sim_config, bc_periodic
  use bem_field_solver, only: field_solver_type
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config, &
                                      outer_plasma_config
  use bem_string_utils, only: lower_ascii
  use bem_periodic_zero_mode_plan, only: periodic_zero_mode_plan_type, periodic_zero_mode_state_type, &
                                         build_periodic_zero_mode_plan, refresh_periodic_zero_mode_state, &
                                         symmetric_vacuum_bottom_field, periodic_zero_mode_ok
  use bem_periodic_zero_mode_eval, only: eval_periodic_zero_mode, zero_mode_trace_plus
  use bem_coulomb_fmm_periodic_nonzero_reference, only: eval_periodic_nonzero_panel_reference
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok
  use bem_outer_plasma_kinetic, only: kinetic_outer_plasma_options_type, solve_outer_plasma_kinetic
  use bem_outer_plasma_zhao, only: zhao_continuation_diagnostics_type, &
                                   zhao_continuation_reason_none, &
                                   write_zhao_continuation_diagnostics
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  use bem_mpi, only: mpi_context, mpi_is_root, mpi_bcast_i32_array, mpi_bcast_real_dp_array, &
                     mpi_world_barrier
  implicit none
  private

  type, public :: electrostatic_diagnostics_type
    logical :: split_periodic_active = .false.
    logical :: applicable = .true.
    integer(i32) :: interface_sample_n = 0_i32
    real(dp) :: interface_potential = 0.0_dp
    real(dp) :: interface_field = 0.0_dp
    real(dp) :: eta_phi_kneq0 = 0.0_dp
    real(dp) :: eta_field_kneq0 = 0.0_dp
    real(dp) :: eta_gap = 0.0_dp
    real(dp) :: eta_local_charge = 0.0_dp
    real(dp) :: gauss_residual = 0.0_dp
    real(dp) :: outer_integrated_charge = 0.0_dp
    integer(i32) :: outer_nonlinear_iterations = 0_i32
    real(dp) :: outer_nonlinear_residual = 0.0_dp
    integer(i32) :: outer_applicability_status = 0_i32
    real(dp) :: outer_infinity_potential = 0.0_dp
    real(dp) :: outer_debye_length = 0.0_dp
    real(dp) :: outer_integrated_charge_per_area = 0.0_dp
    real(dp) :: outer_electron_current_density = 0.0_dp
    real(dp) :: outer_ion_current_density = 0.0_dp
    real(dp) :: outer_photoelectron_current_density = 0.0_dp
    real(dp) :: outer_total_current_density = 0.0_dp
    character(len=32) :: outer_kinetic_closure = 'none'
    character(len=1) :: outer_zhao_branch = ' '
    real(dp) :: outer_zhao_phi0 = 0.0_dp
    real(dp) :: outer_zhao_phi_minimum = 0.0_dp
    real(dp) :: outer_zhao_electron_density_infinity = 0.0_dp
    real(dp) :: outer_photoelectron_source_scale = 1.0_dp
    real(dp) :: outer_photoelectron_population_fraction = 1.0_dp
    real(dp) :: outer_photoelectron_column_per_area = 0.0_dp
    real(dp) :: outer_photoelectron_column_target_per_area = 0.0_dp
    real(dp) :: outer_photoelectron_column_residual_per_area = 0.0_dp
    integer(i64) :: outer_queue_event_count = 0_i64
    real(dp) :: outer_queue_signed_charge = 0.0_dp
    character(len=16) :: outer_queue_fingerprint = ''
    character(len=64) :: status = 'legacy_or_not_applicable'
    integer(i32) :: last_outer_update_batch = -1_i32
    real(dp) :: max_outer_flight_time = 0.0_dp
    real(dp) :: max_frozen_field_ratio = 0.0_dp
    real(dp) :: max_outer_energy_relative_error = 0.0_dp
    logical :: implicit_mean_shadow_diagnostics_available = .false.
    real(dp) :: implicit_mean_last_returned_outer_flight_time_mean = 0.0_dp
    real(dp) :: implicit_mean_last_returning_pe_column_charge_per_area = 0.0_dp
    logical :: periodic_cache_hit = .false.
    integer(i32) :: periodic_operator_build_count = 0_i32
    character(len=128) :: periodic_cache_fingerprint = ''
    character(len=512) :: periodic_cache_path = ''
    real(dp), allocatable :: outer_profile_z(:)
    real(dp), allocatable :: outer_profile_potential(:)
    real(dp), allocatable :: outer_profile_field(:)
    real(dp), allocatable :: outer_profile_charge_density(:)
  end type electrostatic_diagnostics_type

  type, public :: electrostatic_restart_state_type
    integer(i32) :: checkpoint_schema_version = 0_i32
    logical :: outer_ready = .false.
    logical :: outer_profile_complete = .false.
    real(dp) :: outer_interface_potential = 0.0_dp
    real(dp) :: outer_interface_field = 0.0_dp
    integer(i32) :: outer_applicability_status = 0_i32
    integer(i32) :: outer_nonlinear_iterations = 0_i32
    real(dp) :: outer_nonlinear_residual = 0.0_dp
    real(dp) :: outer_infinity_potential = 0.0_dp
    real(dp) :: outer_debye_length = 0.0_dp
    real(dp) :: outer_integrated_charge_per_area = 0.0_dp
    real(dp) :: outer_electron_current_density = 0.0_dp
    real(dp) :: outer_ion_current_density = 0.0_dp
    real(dp) :: outer_photoelectron_current_density = 0.0_dp
    real(dp) :: outer_total_current_density = 0.0_dp
    character(len=1) :: outer_zhao_branch = ' '
    real(dp) :: outer_zhao_phi0 = 0.0_dp
    real(dp) :: outer_zhao_phi_minimum = 0.0_dp
    real(dp) :: outer_zhao_electron_density_infinity = 0.0_dp
    real(dp) :: outer_photoelectron_source_scale = 1.0_dp
    real(dp) :: outer_photoelectron_population_fraction = 1.0_dp
    real(dp) :: outer_photoelectron_column_per_area = 0.0_dp
    real(dp) :: outer_photoelectron_column_target_per_area = 0.0_dp
    real(dp) :: outer_photoelectron_column_residual_per_area = 0.0_dp
    integer(i64) :: outer_queue_event_count = 0_i64
    real(dp) :: outer_queue_signed_charge = 0.0_dp
    character(len=16) :: outer_queue_fingerprint = ''
    logical :: outer_zhao_state_complete = .false.
    logical :: outer_zhao_source_scale_complete = .false.
    logical :: outer_zhao_transient_state_complete = .false.
    logical :: outer_queue_inventory_complete = .false.
    real(dp) :: max_outer_flight_time = 0.0_dp
    real(dp) :: max_frozen_field_ratio = 0.0_dp
    real(dp) :: max_outer_energy_relative_error = 0.0_dp
    logical :: outer_max_diagnostics_complete = .false.
    integer(i32) :: last_outer_update_batch = -1_i32
    real(dp), allocatable :: outer_profile_z(:)
    real(dp), allocatable :: outer_profile_potential(:)
    real(dp), allocatable :: outer_profile_field(:)
    real(dp), allocatable :: outer_profile_charge_density(:)
  end type electrostatic_restart_state_type

  type, public :: electrostatic_snapshot_type
    type(field_solver_type) :: nonzero_solver
    real(dp) :: prescribed_e(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: prescribed_phi_origin(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    logical :: use_zero_mode = .false.
    logical :: use_outer_plasma = .false.
    logical :: use_panel_spectral_reference = .false.
    logical :: use_cached_kneq0 = .false.
    real(dp) :: periodic_length(2) = [0.0_dp, 0.0_dp]
    real(dp) :: periodic_origin(2) = [0.0_dp, 0.0_dp]
    integer(i32) :: reference_mode_layers = 0_i32
    integer(i32) :: panel_quadrature_order = 0_i32
    type(periodic_zero_mode_plan_type) :: zero_plan
    type(periodic_zero_mode_state_type) :: zero_state
    type(outer_plasma_state_type) :: outer
    type(outer_plasma_config) :: outer_options
    type(kinetic_outer_plasma_options_type) :: kinetic_options
    type(mpi_context) :: mpi
    type(periodic2_physics_config) :: periodic_options
    type(electrostatic_diagnostics_type) :: diagnostics
    real(dp) :: mesh_top_z = 0.0_dp
    real(dp) :: gauss_residual = 0.0_dp
  contains
    procedure :: init => init_electrostatic_snapshot
    procedure :: refresh => refresh_electrostatic_snapshot
    procedure :: refresh_mean_only => refresh_electrostatic_snapshot_mean_only
    procedure :: adopt_mean_outer => adopt_electrostatic_snapshot_mean_outer
    procedure :: eval_local_e => eval_snapshot_local_e
    procedure :: eval_local_phi => eval_snapshot_local_phi
    procedure :: eval_local_phi_without_primary_self => eval_snapshot_local_phi_without_primary_self
    procedure :: compute_mesh_potential => compute_snapshot_mesh_potential
    procedure :: get_diagnostics => get_snapshot_diagnostics
    procedure :: restore_outer_state => restore_snapshot_outer_state
    procedure :: export_restart_state => export_snapshot_restart_state
  end type electrostatic_snapshot_type

  interface
    !> 観測点の局所電場を、選択中の非零モードと零モードから評価する。
    module subroutine eval_snapshot_local_e(self, mesh, position, electric_field)
      class(electrostatic_snapshot_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
      real(dp), intent(in) :: position(3)
      real(dp), intent(out) :: electric_field(3)
    end subroutine eval_snapshot_local_e

    !> 観測点の局所電位を、選択中の非零モードと零モードから評価する。
    module subroutine eval_snapshot_local_phi(self, mesh, sim, position, potential)
      class(electrostatic_snapshot_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
      type(sim_config), intent(in) :: sim
      real(dp), intent(in) :: position(3)
      real(dp), intent(out) :: potential
    end subroutine eval_snapshot_local_phi

    !> 指定要素の主自己寄与を除いた重心電位を評価する。
    module subroutine eval_snapshot_local_phi_without_primary_self(self, mesh, sim, element, potential)
      class(electrostatic_snapshot_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
      type(sim_config), intent(in) :: sim
      integer(i32), intent(in) :: element
      real(dp), intent(out) :: potential
    end subroutine eval_snapshot_local_phi_without_primary_self

    !> 全要素重心での電位を評価する。
    module subroutine compute_snapshot_mesh_potential(self, mesh, sim, potential)
      class(electrostatic_snapshot_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
      type(sim_config), intent(in) :: sim
      real(dp), intent(out) :: potential(:)
    end subroutine compute_snapshot_mesh_potential

  end interface

contains

  subroutine init_electrostatic_snapshot( &
    self, mesh, sim, field_config, periodic_config, panel_config, outer_config, kinetic_options, mpi &
    )
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(field_physics_config), intent(in), optional :: field_config
    type(periodic2_physics_config), intent(in), optional :: periodic_config
    type(panel_kernel_config), intent(in), optional :: panel_config
    type(outer_plasma_config), intent(in), optional :: outer_config
    type(kinetic_outer_plasma_options_type), intent(in), optional :: kinetic_options
    type(mpi_context), intent(in), optional :: mpi

    self%prescribed_e = sim%e0
    self%prescribed_phi_origin = 0.0_dp
    self%use_zero_mode = .false.
    self%use_outer_plasma = .false.
    self%use_panel_spectral_reference = .false.
    self%use_cached_kneq0 = .false.
    self%gauss_residual = 0.0_dp
    self%diagnostics = electrostatic_diagnostics_type()
    self%mpi = mpi_context()
    if (present(mpi)) self%mpi = mpi
    if (present(periodic_config) .and. present(panel_config) .and. present(outer_config)) then
      if (trim(lower_ascii(periodic_config%nonzero_mode_backend)) == 'panel_spectral_reference') then
        call init_split_periodic_snapshot(self, mesh, sim, periodic_config, outer_config, kinetic_options)
        return
      else if (trim(lower_ascii(periodic_config%nonzero_mode_backend)) == 'cached_kneq0') then
        if (.not. present(field_config)) error stop 'cached_kneq0 requires typed field configuration.'
        call init_cached_periodic_snapshot( &
          self, mesh, sim, field_config, periodic_config, panel_config, outer_config, kinetic_options &
          )
        return
      end if
    end if
    if (present(outer_config)) then
      if (trim(lower_ascii(outer_config%model)) /= 'none') then
        error stop 'outer plasma requires the split periodic snapshot path.'
      end if
    end if
    if (present(field_config) .and. present(periodic_config) .and. present(panel_config)) then
      call self%nonzero_solver%init(mesh, sim, field_config, periodic_config, panel_config)
    else
      call self%nonzero_solver%init(mesh, sim)
    end if
  end subroutine init_electrostatic_snapshot

  subroutine refresh_electrostatic_snapshot( &
    self, mesh, update_outer, continuation_stage, continuation_batch &
    )
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    logical, intent(in), optional :: update_outer
    character(len=*), intent(in), optional :: continuation_stage
    integer(i32), intent(in), optional :: continuation_batch

    real(dp) :: interface_potential, interface_field, raw_potential
    logical :: refresh_outer

    if (self%use_cached_kneq0 .and. .not. self%use_outer_plasma) then
      call self%nonzero_solver%refresh(mesh)
      call refresh_periodic_zero_mode_state( &
        self%zero_plan, mesh%q_elem, zero_mode_bottom_field(self, mesh%q_elem), &
        self%zero_plan%break_z(1), 0.0_dp, self%zero_state &
        )
      return
    else if (self%use_cached_kneq0) then
      call refresh_cached_kinetic_outer( &
        self, mesh, update_outer, continuation_stage, continuation_batch &
        )
      return
    else if (.not. self%use_panel_spectral_reference) then
      call self%nonzero_solver%refresh(mesh)
      return
    end if
    refresh_outer = .true.
    if (present(update_outer)) refresh_outer = update_outer
    if (.not. self%outer%ready) refresh_outer = .true.
    call refresh_periodic_zero_mode_state( &
      self%zero_plan, mesh%q_elem, zero_mode_bottom_field(self, mesh%q_elem), &
      self%zero_plan%break_z(1), 0.0_dp, self%zero_state &
      )
    call eval_periodic_zero_mode( &
      self%zero_plan, self%zero_state, self%outer_options%interface_z, zero_mode_trace_plus, &
      raw_potential, interface_field &
      )
    if (refresh_outer) then
      call solve_kinetic_collective( &
        self, interface_field, continuation_stage, continuation_batch &
        )
      interface_potential = self%outer%interface_potential
    else
      interface_potential = self%outer%interface_potential
    end if
    call refresh_periodic_zero_mode_state( &
      self%zero_plan, mesh%q_elem, zero_mode_bottom_field(self, mesh%q_elem), self%zero_plan%break_z(1), &
      interface_potential - raw_potential, self%zero_state &
      )
    self%gauss_residual = upper_flux_charge(self) + &
                          self%zero_plan%area_xy*(-eps0*self%outer%interface_field)
    call update_interface_diagnostics(self, mesh)
  end subroutine refresh_electrostatic_snapshot

  !> 候補要素電荷から平面平均成分と外部1D状態だけを更新する。
  !!
  !! nonzero_solver は更新しないため、局所 k/=0 場は直前の full refresh の
  !! 状態に固定される。candidate_charge は差分ではなく全要素分の絶対電荷。
  subroutine refresh_electrostatic_snapshot_mean_only( &
    self, mesh, candidate_charge, continuation_stage, continuation_batch &
    )
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: candidate_charge(:)
    character(len=*), intent(in), optional :: continuation_stage
    integer(i32), intent(in), optional :: continuation_batch

    real(dp) :: raw_potential, interface_field

    if (.not. self%use_zero_mode .or. .not. self%use_outer_plasma .or. &
        (.not. self%use_cached_kneq0 .and. .not. self%use_panel_spectral_reference)) then
      error stop 'mean-only refresh requires a split periodic snapshot with an outer plasma.'
    end if
    if (size(candidate_charge) /= mesh%nelem) then
      error stop 'mean-only refresh requires one candidate charge per mesh element.'
    end if
    if (.not. all(ieee_is_finite(candidate_charge))) then
      error stop 'mean-only refresh candidate charges must be finite.'
    end if

    call refresh_periodic_zero_mode_state( &
      self%zero_plan, candidate_charge, zero_mode_bottom_field(self, candidate_charge), &
      self%zero_plan%break_z(1), 0.0_dp, self%zero_state &
      )
    call eval_periodic_zero_mode( &
      self%zero_plan, self%zero_state, self%outer_options%interface_z, zero_mode_trace_plus, &
      raw_potential, interface_field &
      )
    call solve_kinetic_collective( &
      self, interface_field, continuation_stage, continuation_batch &
      )
    call refresh_periodic_zero_mode_state( &
      self%zero_plan, candidate_charge, zero_mode_bottom_field(self, candidate_charge), &
      self%zero_plan%break_z(1), self%outer%interface_potential - raw_potential, self%zero_state &
      )
    self%gauss_residual = upper_flux_charge(self) + &
                          self%zero_plan%area_xy*(-eps0*self%outer%interface_field)
    call update_interface_diagnostics(self, mesh, candidate_charge)
  end subroutine refresh_electrostatic_snapshot_mean_only

  !> root rankで解決済みの外部profileを再solveせず、候補mean chargeへcollectiveに採用する。
  subroutine adopt_electrostatic_snapshot_mean_outer(self, mesh, candidate_charge, resolved_outer)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: candidate_charge(:)
    type(outer_plasma_state_type), intent(inout) :: resolved_outer

    integer(i32) :: status_values(4)
    real(dp) :: scalar_values(19)
    real(dp) :: raw_potential, interface_field, field_tolerance, z_tolerance

    if (.not. self%use_zero_mode .or. .not. self%use_outer_plasma .or. &
        (.not. self%use_cached_kneq0 .and. .not. self%use_panel_spectral_reference)) then
      error stop 'resolved mean outer adoption requires a split periodic outer-plasma snapshot.'
    end if
    if (size(candidate_charge) /= mesh%nelem .or. &
        .not. all(ieee_is_finite(candidate_charge))) then
      error stop 'resolved mean outer adoption received invalid candidate charges.'
    end if

    call refresh_periodic_zero_mode_state( &
      self%zero_plan, candidate_charge, zero_mode_bottom_field(self, candidate_charge), &
      self%zero_plan%break_z(1), 0.0_dp, self%zero_state &
      )
    call eval_periodic_zero_mode( &
      self%zero_plan, self%zero_state, self%outer_options%interface_z, zero_mode_trace_plus, &
      raw_potential, interface_field &
      )

    status_values = 0_i32
    scalar_values = 0.0_dp
    if (mpi_is_root(self%mpi)) then
      status_values(1) = outer_plasma_ok
      if (.not. resolved_outer%ready .or. &
          trim(lower_ascii(resolved_outer%kinetic_closure)) /= 'zhao_charge_driven' .or. &
          index('ABC', resolved_outer%zhao_branch) == 0 .or. &
          resolved_outer%profile_n /= self%kinetic_options%grid_points .or. &
          resolved_outer%photoelectron_source_scale <= 0.0_dp) then
        status_values(1) = 1_i32
      else if (.not. valid_resolved_zhao_outer(resolved_outer, self%kinetic_options)) then
        status_values(1) = 6_i32
      end if
      if (status_values(1) == outer_plasma_ok) then
        field_tolerance = 1.0e-10_dp*max(1.0_dp, abs(interface_field))
        if (abs(resolved_outer%interface_field - interface_field) > field_tolerance) then
          status_values(1) = 4_i32
        end if
      end if
      if (status_values(1) == outer_plasma_ok) then
        z_tolerance = 256.0_dp*epsilon(1.0_dp)*max( &
                      1.0_dp, abs(self%outer_options%interface_z), &
                      abs(resolved_outer%z(resolved_outer%profile_n)) &
                      )
        if (abs(resolved_outer%z(1)) <= z_tolerance) then
          resolved_outer%z = resolved_outer%z + self%outer_options%interface_z
        else if (abs(resolved_outer%z(1) - self%outer_options%interface_z) > z_tolerance) then
          status_values(1) = 5_i32
        end if
        resolved_outer%interface_z = self%outer_options%interface_z
      end if
      status_values(2) = resolved_outer%profile_n
      status_values(3) = resolved_outer%nonlinear_iterations
      status_values(4) = int(iachar(resolved_outer%zhao_branch), i32)
    end if
    call mpi_bcast_i32_array(self%mpi, status_values, 0_i32)
    if (status_values(1) /= outer_plasma_ok) then
      call mpi_world_barrier(self%mpi)
      error stop 'resolved mean Zhao outer state failed collective adoption.'
    end if

    if (.not. mpi_is_root(self%mpi)) then
      resolved_outer = outer_plasma_state_type()
      resolved_outer%profile_n = status_values(2)
      resolved_outer%nonlinear_iterations = status_values(3)
      resolved_outer%zhao_branch = achar(status_values(4))
      allocate (resolved_outer%z(resolved_outer%profile_n))
      allocate (resolved_outer%potential(resolved_outer%profile_n))
      allocate (resolved_outer%field(resolved_outer%profile_n))
      allocate (resolved_outer%charge_density(resolved_outer%profile_n))
    end if
    if (mpi_is_root(self%mpi)) then
      scalar_values = [ &
                      resolved_outer%interface_potential, resolved_outer%infinity_potential, &
                      resolved_outer%debye_length, resolved_outer%interface_field, &
                      resolved_outer%nonlinear_residual, resolved_outer%integrated_charge_per_area, &
                      merge(1.0_dp, 0.0_dp, resolved_outer%ready), &
                      resolved_outer%electron_current_density, resolved_outer%ion_current_density, &
                      resolved_outer%photoelectron_current_density, resolved_outer%total_current_density, &
                      resolved_outer%zhao_phi0, resolved_outer%zhao_phi_minimum, &
                      resolved_outer%zhao_electron_density_infinity, &
                      resolved_outer%photoelectron_population_fraction, &
                      resolved_outer%photoelectron_column_per_area, &
                      resolved_outer%photoelectron_column_target_per_area, &
                      resolved_outer%photoelectron_column_residual_per_area, &
                      resolved_outer%photoelectron_source_scale &
                      ]
    end if
    call mpi_bcast_real_dp_array(self%mpi, scalar_values, 0_i32)
    call mpi_bcast_real_dp_array(self%mpi, resolved_outer%z, 0_i32)
    call mpi_bcast_real_dp_array(self%mpi, resolved_outer%potential, 0_i32)
    call mpi_bcast_real_dp_array(self%mpi, resolved_outer%field, 0_i32)
    call mpi_bcast_real_dp_array(self%mpi, resolved_outer%charge_density, 0_i32)

    resolved_outer%interface_z = self%outer_options%interface_z
    resolved_outer%interface_potential = scalar_values(1)
    resolved_outer%infinity_potential = scalar_values(2)
    resolved_outer%debye_length = scalar_values(3)
    resolved_outer%interface_field = scalar_values(4)
    resolved_outer%nonlinear_residual = scalar_values(5)
    resolved_outer%integrated_charge_per_area = scalar_values(6)
    resolved_outer%ready = scalar_values(7) > 0.5_dp
    resolved_outer%electron_current_density = scalar_values(8)
    resolved_outer%ion_current_density = scalar_values(9)
    resolved_outer%photoelectron_current_density = scalar_values(10)
    resolved_outer%total_current_density = scalar_values(11)
    resolved_outer%zhao_phi0 = scalar_values(12)
    resolved_outer%zhao_phi_minimum = scalar_values(13)
    resolved_outer%zhao_electron_density_infinity = scalar_values(14)
    resolved_outer%photoelectron_population_fraction = scalar_values(15)
    resolved_outer%photoelectron_column_per_area = scalar_values(16)
    resolved_outer%photoelectron_column_target_per_area = scalar_values(17)
    resolved_outer%photoelectron_column_residual_per_area = scalar_values(18)
    resolved_outer%photoelectron_source_scale = scalar_values(19)
    resolved_outer%model = 'kinetic_1d'
    resolved_outer%kinetic_closure = 'zhao_charge_driven'
    resolved_outer%zhao_branch = achar(status_values(4))
    resolved_outer%applicability_status = outer_plasma_ok

    self%outer = resolved_outer
    self%kinetic_options%interface_field = resolved_outer%interface_field
    self%kinetic_options%photoelectron_source_scale = resolved_outer%photoelectron_source_scale
    self%kinetic_options%photoelectron_population_fraction = resolved_outer%photoelectron_population_fraction
    self%kinetic_options%photoelectron_column_closure_enabled = .false.
    self%kinetic_options%zhao_branch = lower_ascii(resolved_outer%zhao_branch)
    call refresh_periodic_zero_mode_state( &
      self%zero_plan, candidate_charge, zero_mode_bottom_field(self, candidate_charge), &
      self%zero_plan%break_z(1), self%outer%interface_potential - raw_potential, self%zero_state &
      )
    self%gauss_residual = upper_flux_charge(self) + &
                          self%zero_plan%area_xy*(-eps0*self%outer%interface_field)
    call update_interface_diagnostics(self, mesh, candidate_charge)
  end subroutine adopt_electrostatic_snapshot_mean_outer

  subroutine get_snapshot_diagnostics(self, diagnostics)
    class(electrostatic_snapshot_type), intent(in) :: self
    type(electrostatic_diagnostics_type), intent(out) :: diagnostics

    diagnostics = self%diagnostics
    diagnostics%outer_nonlinear_iterations = self%outer%nonlinear_iterations
    diagnostics%outer_nonlinear_residual = self%outer%nonlinear_residual
    diagnostics%outer_applicability_status = self%outer%applicability_status
    diagnostics%outer_infinity_potential = self%outer%infinity_potential
    diagnostics%outer_debye_length = self%outer%debye_length
    diagnostics%outer_integrated_charge_per_area = self%outer%integrated_charge_per_area
    diagnostics%outer_electron_current_density = self%outer%electron_current_density
    diagnostics%outer_ion_current_density = self%outer%ion_current_density
    diagnostics%outer_photoelectron_current_density = self%outer%photoelectron_current_density
    diagnostics%outer_total_current_density = self%outer%total_current_density
    diagnostics%outer_kinetic_closure = self%outer%kinetic_closure
    diagnostics%outer_zhao_branch = self%outer%zhao_branch
    diagnostics%outer_zhao_phi0 = self%outer%zhao_phi0
    diagnostics%outer_zhao_phi_minimum = self%outer%zhao_phi_minimum
    diagnostics%outer_zhao_electron_density_infinity = self%outer%zhao_electron_density_infinity
    diagnostics%outer_photoelectron_source_scale = self%outer%photoelectron_source_scale
    diagnostics%outer_photoelectron_population_fraction = self%outer%photoelectron_population_fraction
    diagnostics%outer_photoelectron_column_per_area = self%outer%photoelectron_column_per_area
    diagnostics%outer_photoelectron_column_target_per_area = self%outer%photoelectron_column_target_per_area
    diagnostics%outer_photoelectron_column_residual_per_area = self%outer%photoelectron_column_residual_per_area
    if (allocated(self%outer%z) .and. allocated(self%outer%potential)) then
      diagnostics%outer_profile_z = self%outer%z
      diagnostics%outer_profile_potential = self%outer%potential
      if (allocated(self%outer%field)) diagnostics%outer_profile_field = self%outer%field
      if (allocated(self%outer%charge_density)) then
        diagnostics%outer_profile_charge_density = self%outer%charge_density
      end if
    end if
  end subroutine get_snapshot_diagnostics

  !> Zhao profile の採用・restartで共通に使う副作用なしの内部契約検査。
  pure logical function valid_resolved_zhao_outer(state, options) result(valid)
    type(outer_plasma_state_type), intent(in) :: state
    type(kinetic_outer_plasma_options_type), intent(in) :: options

    real(dp) :: voltage_scale, potential_tolerance, field_tolerance
    real(dp) :: asymptotic_tolerance, shape_tolerance, field_shape_tolerance
    integer(i32) :: minimum_point

    valid = .false.
    if (.not. state%ready .or. state%applicability_status /= outer_plasma_ok .or. &
        trim(lower_ascii(state%model)) /= 'kinetic_1d' .or. &
        trim(lower_ascii(state%kinetic_closure)) /= 'zhao_charge_driven' .or. &
        index('ABC', state%zhao_branch) == 0 .or. state%profile_n /= options%grid_points .or. &
        state%profile_n < 5_i32) return
    if (.not. allocated(state%z) .or. .not. allocated(state%potential) .or. &
        .not. allocated(state%field) .or. .not. allocated(state%charge_density)) return
    if (size(state%z) /= state%profile_n .or. size(state%potential) /= state%profile_n .or. &
        size(state%field) /= state%profile_n .or. size(state%charge_density) /= state%profile_n) return
    if (.not. all(ieee_is_finite([ &
                                 state%interface_z, state%interface_potential, state%infinity_potential, &
                                 state%debye_length, state%interface_field, state%nonlinear_residual, &
                                 state%integrated_charge_per_area, state%electron_current_density, &
                                 state%ion_current_density, state%photoelectron_current_density, &
                                 state%total_current_density, state%zhao_phi0, state%zhao_phi_minimum, &
                                 state%zhao_electron_density_infinity, state%photoelectron_population_fraction, &
                                 state%photoelectron_column_per_area, state%photoelectron_column_target_per_area, &
                                 state%photoelectron_column_residual_per_area, state%photoelectron_source_scale &
                                 ]))) return
    if (.not. all(ieee_is_finite(state%z)) .or. .not. all(ieee_is_finite(state%potential)) .or. &
        .not. all(ieee_is_finite(state%field)) .or. .not. all(ieee_is_finite(state%charge_density))) return
    if (.not. ieee_is_finite(options%residual_tolerance) .or. options%residual_tolerance <= 0.0_dp) return
    if (state%debye_length <= 0.0_dp .or. state%nonlinear_iterations < 0_i32 .or. &
        state%nonlinear_residual < 0.0_dp .or. &
        state%nonlinear_residual > max(options%residual_tolerance, 128.0_dp*epsilon(1.0_dp)) .or. &
        state%zhao_electron_density_infinity <= 0.0_dp .or. &
        state%photoelectron_source_scale < 0.0_dp .or. &
        state%photoelectron_population_fraction < 0.0_dp .or. &
        state%photoelectron_population_fraction > 1.0_dp .or. &
        state%photoelectron_column_per_area < 0.0_dp .or. &
        state%photoelectron_column_target_per_area < 0.0_dp) return
    if (state%photoelectron_source_scale == 0.0_dp .and. &
        state%photoelectron_population_fraction > 64.0_dp*epsilon(1.0_dp)) return
    if (any(state%z(2:) <= state%z(:state%profile_n - 1_i32))) return
    if (.not. ieee_is_finite(options%photoelectron_temperature_j) .or. &
        .not. ieee_is_finite(options%photoelectron_charge) .or. &
        options%photoelectron_temperature_j <= 0.0_dp .or. options%photoelectron_charge >= 0.0_dp) return

    voltage_scale = options%photoelectron_temperature_j/abs(options%photoelectron_charge)
    potential_tolerance = 1.0e-10_dp*max(1.0_dp, abs(state%interface_potential))
    field_tolerance = 1.0e-10_dp*max(1.0_dp, abs(state%interface_field))
    asymptotic_tolerance = max(potential_tolerance, 2.0e-4_dp*voltage_scale)
    shape_tolerance = max(potential_tolerance, 1.0e-9_dp*voltage_scale)
    field_shape_tolerance = max(field_tolerance, 1.0e-9_dp*voltage_scale/state%debye_length)
    if (abs(state%potential(1) - state%interface_potential) > potential_tolerance .or. &
        abs(state%field(1) - state%interface_field) > field_tolerance .or. &
        abs(state%zhao_phi0 - state%interface_potential) > potential_tolerance .or. &
        abs(state%infinity_potential) > potential_tolerance .or. &
        abs(state%potential(state%profile_n) - state%infinity_potential) > asymptotic_tolerance .or. &
        state%z(state%profile_n) <= state%z(1)) return

    select case (state%zhao_branch)
    case ('A')
      minimum_point = minloc(state%potential, dim=1)
      if (minimum_point <= 1_i32 .or. minimum_point >= state%profile_n) return
      if (abs(state%potential(minimum_point) - state%zhao_phi_minimum) > potential_tolerance) return
      if (any(state%potential(2:minimum_point) > &
              state%potential(1:minimum_point - 1_i32) + shape_tolerance)) return
      if (any(state%potential(minimum_point + 1_i32:) < &
              state%potential(minimum_point:state%profile_n - 1_i32) - shape_tolerance)) return
      if (any(state%field(:minimum_point - 1_i32) < -field_shape_tolerance) .or. &
          any(state%field(minimum_point + 1_i32:) > field_shape_tolerance)) return
    case ('B')
      if (abs(state%zhao_phi_minimum - state%zhao_phi0) > potential_tolerance) return
      if (any(state%potential < -shape_tolerance) .or. &
          any(state%potential(2:) > state%potential(:state%profile_n - 1_i32) + shape_tolerance) .or. &
          any(state%field < -field_shape_tolerance)) return
    case ('C')
      if (abs(state%zhao_phi_minimum - state%zhao_phi0) > potential_tolerance) return
      if (any(state%potential > shape_tolerance) .or. &
          any(state%potential(2:) < state%potential(:state%profile_n - 1_i32) - shape_tolerance) .or. &
          any(state%field > field_shape_tolerance)) return
    end select
    valid = .true.
  end function valid_resolved_zhao_outer

  subroutine restore_snapshot_outer_state(self, state, require_charge_consistency)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(electrostatic_restart_state_type), intent(in) :: state
    logical, intent(in), optional :: require_charge_consistency
    real(dp) :: raw_potential, restored_interface_potential, expected_interface_field
    real(dp) :: potential_tolerance, field_tolerance, z_tolerance
    logical :: require_consistent_charge

    require_consistent_charge = .false.
    if (present(require_charge_consistency)) require_consistent_charge = require_charge_consistency
    if (.not. state%outer_ready) return
    if (.not. self%use_outer_plasma) error stop 'checkpoint outer state is incompatible with the active model.'
    if (trim(lower_ascii(self%outer_options%model)) == 'kinetic_1d') then
      if (.not. allocated(state%outer_profile_z) .or. &
          .not. allocated(state%outer_profile_potential)) then
        error stop 'checkpoint kinetic outer profile is missing or has an incompatible grid.'
      end if
      if (size(state%outer_profile_z) /= self%kinetic_options%grid_points .or. &
          size(state%outer_profile_potential) /= self%kinetic_options%grid_points) then
        error stop 'checkpoint kinetic outer profile is missing or has an incompatible grid.'
      end if
      self%outer = outer_plasma_state_type()
      self%outer%model = 'kinetic_1d'
      self%outer%kinetic_closure = self%outer_options%kinetic_closure
      self%outer%ready = state%outer_ready .and. state%outer_profile_complete
      self%outer%applicability_status = state%outer_applicability_status
      self%outer%profile_n = self%kinetic_options%grid_points
      self%outer%interface_z = self%outer_options%interface_z
      self%outer%interface_potential = state%outer_interface_potential
      self%outer%interface_field = state%outer_interface_field
      self%outer%infinity_potential = state%outer_infinity_potential
      self%outer%debye_length = state%outer_debye_length
      self%outer%nonlinear_iterations = state%outer_nonlinear_iterations
      self%outer%nonlinear_residual = state%outer_nonlinear_residual
      self%outer%integrated_charge_per_area = state%outer_integrated_charge_per_area
      self%outer%electron_current_density = state%outer_electron_current_density
      self%outer%ion_current_density = state%outer_ion_current_density
      self%outer%photoelectron_current_density = state%outer_photoelectron_current_density
      self%outer%total_current_density = state%outer_total_current_density
      self%outer%zhao_branch = state%outer_zhao_branch
      self%outer%zhao_phi0 = state%outer_zhao_phi0
      self%outer%zhao_phi_minimum = state%outer_zhao_phi_minimum
      self%outer%zhao_electron_density_infinity = state%outer_zhao_electron_density_infinity
      if (state%outer_zhao_source_scale_complete) then
        self%outer%photoelectron_source_scale = state%outer_photoelectron_source_scale
        self%kinetic_options%photoelectron_source_scale = state%outer_photoelectron_source_scale
        if (.not. self%kinetic_options%photoelectron_column_closure_enabled .and. &
            index('ABC', state%outer_zhao_branch) > 0) then
          self%kinetic_options%zhao_branch = lower_ascii(state%outer_zhao_branch)
        end if
      else
        self%outer%photoelectron_source_scale = self%kinetic_options%photoelectron_source_scale
      end if
      self%outer%photoelectron_population_fraction = state%outer_photoelectron_population_fraction
      self%outer%photoelectron_column_per_area = state%outer_photoelectron_column_per_area
      self%outer%photoelectron_column_target_per_area = state%outer_photoelectron_column_target_per_area
      self%outer%photoelectron_column_residual_per_area = state%outer_photoelectron_column_residual_per_area
      self%kinetic_options%photoelectron_population_fraction = state%outer_photoelectron_population_fraction
      if (self%kinetic_options%photoelectron_column_closure_enabled) then
        self%kinetic_options%photoelectron_column_target_m2 = state%outer_photoelectron_column_target_per_area
      end if
      self%outer%z = state%outer_profile_z
      self%outer%potential = state%outer_profile_potential
      if (state%outer_profile_complete) then
        if (.not. allocated(state%outer_profile_field) .or. &
            .not. allocated(state%outer_profile_charge_density)) then
          error stop 'checkpoint kinetic outer profile is marked complete but E/rho arrays are missing.'
        end if
        if (size(state%outer_profile_field) /= self%kinetic_options%grid_points .or. &
            size(state%outer_profile_charge_density) /= self%kinetic_options%grid_points) then
          error stop 'checkpoint kinetic outer profile is marked complete but E/rho arrays are missing.'
        end if
        self%outer%field = state%outer_profile_field
        self%outer%charge_density = state%outer_profile_charge_density
      end if
      if (trim(lower_ascii(self%outer%kinetic_closure)) == 'zhao_charge_driven') then
        if (.not. valid_resolved_zhao_outer(self%outer, self%kinetic_options)) then
          error stop 'checkpoint Zhao outer profile failed its internal state contract.'
        end if
        z_tolerance = 256.0_dp*epsilon(1.0_dp)*max( &
                      1.0_dp, abs(self%outer_options%interface_z), &
                      abs(self%outer%z(self%outer%profile_n)) &
                      )
        if (abs(self%outer%z(1) - self%outer_options%interface_z) > z_tolerance) then
          error stop 'checkpoint Zhao outer profile does not start at the configured interface.'
        end if
        call eval_periodic_zero_mode( &
          self%zero_plan, self%zero_state, self%outer_options%interface_z, zero_mode_trace_plus, &
          raw_potential, expected_interface_field &
          )
        if (require_consistent_charge) then
          field_tolerance = 1.0e-10_dp*max(1.0_dp, abs(expected_interface_field))
          if (.not. ieee_is_finite(expected_interface_field) .or. &
              abs(self%outer%interface_field - expected_interface_field) > field_tolerance) then
            error stop 'checkpoint implicit Zhao field is inconsistent with the restored surface charge.'
          end if
        end if
        self%zero_state%phi_gauge = self%outer%interface_potential - raw_potential
        call eval_periodic_zero_mode( &
          self%zero_plan, self%zero_state, self%outer_options%interface_z, zero_mode_trace_plus, &
          restored_interface_potential, expected_interface_field &
          )
        potential_tolerance = 1.0e-10_dp*max(1.0_dp, abs(self%outer%interface_potential))
        if (.not. ieee_is_finite(restored_interface_potential) .or. &
            abs(restored_interface_potential - self%outer%interface_potential) > potential_tolerance) then
          error stop 'checkpoint Zhao interface-potential gauge could not be restored.'
        end if
      end if
      return
    end if
    error stop 'checkpoint outer state is incompatible with the active outer-plasma model.'
  end subroutine restore_snapshot_outer_state

  subroutine export_snapshot_restart_state(self, last_outer_update_batch, state)
    class(electrostatic_snapshot_type), intent(in) :: self
    integer(i32), intent(in) :: last_outer_update_batch
    type(electrostatic_restart_state_type), intent(out) :: state

    state%outer_ready = self%outer%ready
    state%outer_profile_complete = self%outer%ready .and. allocated(self%outer%z) .and. &
                                   allocated(self%outer%potential) .and. allocated(self%outer%field) .and. &
                                   allocated(self%outer%charge_density)
    if (self%outer%ready) then
      state%outer_interface_potential = self%outer%interface_potential
      state%outer_interface_field = self%outer%interface_field
      state%outer_applicability_status = self%outer%applicability_status
      state%outer_nonlinear_iterations = self%outer%nonlinear_iterations
      state%outer_nonlinear_residual = self%outer%nonlinear_residual
      state%outer_infinity_potential = self%outer%infinity_potential
      state%outer_debye_length = self%outer%debye_length
      state%outer_integrated_charge_per_area = self%outer%integrated_charge_per_area
      state%outer_electron_current_density = self%outer%electron_current_density
      state%outer_ion_current_density = self%outer%ion_current_density
      state%outer_photoelectron_current_density = self%outer%photoelectron_current_density
      state%outer_total_current_density = self%outer%total_current_density
      state%outer_zhao_branch = self%outer%zhao_branch
      state%outer_zhao_phi0 = self%outer%zhao_phi0
      state%outer_zhao_phi_minimum = self%outer%zhao_phi_minimum
      state%outer_zhao_electron_density_infinity = self%outer%zhao_electron_density_infinity
      state%outer_photoelectron_source_scale = self%outer%photoelectron_source_scale
      state%outer_photoelectron_population_fraction = self%outer%photoelectron_population_fraction
      state%outer_photoelectron_column_per_area = self%outer%photoelectron_column_per_area
      state%outer_photoelectron_column_target_per_area = self%outer%photoelectron_column_target_per_area
      state%outer_photoelectron_column_residual_per_area = self%outer%photoelectron_column_residual_per_area
      state%outer_zhao_state_complete = &
        trim(lower_ascii(self%outer%kinetic_closure)) == 'zhao_charge_driven'
      state%outer_zhao_source_scale_complete = state%outer_zhao_state_complete
      state%outer_zhao_transient_state_complete = state%outer_zhao_state_complete
    end if
    if (state%outer_profile_complete) then
      state%outer_profile_z = self%outer%z
      state%outer_profile_potential = self%outer%potential
      state%outer_profile_field = self%outer%field
      state%outer_profile_charge_density = self%outer%charge_density
    end if
    state%last_outer_update_batch = last_outer_update_batch
  end subroutine export_snapshot_restart_state

  subroutine update_interface_diagnostics(self, mesh, diagnostic_charge)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in), optional :: diagnostic_charge(:)
    type(sim_config) :: evaluation_sim
    integer(i32) :: ix, iy, n
    real(dp) :: position(3), potential, field(3), max_phi, max_field
    real(dp) :: phi_scale, field_scale, gap, surface_scale, surface_total_charge
    real(dp) :: local_charge_estimate

    n = self%periodic_options%interface_sample_n
    evaluation_sim = sim_config()
    max_phi = 0.0_dp
    max_field = 0.0_dp
    do iy = 0_i32, n - 1_i32
      do ix = 0_i32, n - 1_i32
        position = [ &
                   self%periodic_origin(1) + (real(ix, dp) + 0.5_dp)*self%periodic_length(1)/real(n, dp), &
                   self%periodic_origin(2) + (real(iy, dp) + 0.5_dp)*self%periodic_length(2)/real(n, dp), &
                   self%outer_options%interface_z &
                   ]
        if (self%use_cached_kneq0) then
          call self%nonzero_solver%eval_potential(mesh, evaluation_sim, position, potential)
          call self%nonzero_solver%eval_e(mesh, position, field)
        else
          call eval_periodic_nonzero_panel_reference( &
            mesh, position, self%periodic_length(1), self%periodic_length(2), &
            self%reference_mode_layers, self%panel_quadrature_order, potential, field &
            )
        end if
        max_phi = max(max_phi, abs(potential))
        max_field = max(max_field, sqrt(sum(field*field)))
      end do
    end do
    phi_scale = max(abs(self%outer%interface_potential), self%outer_options%thermal_voltage, tiny(1.0_dp))
    field_scale = max(abs(self%outer%interface_field), &
                      self%outer_options%thermal_voltage/minval(self%periodic_length), tiny(1.0_dp))
    gap = self%outer_options%interface_z - self%mesh_top_z
    surface_total_charge = sum(mesh%q_elem)
    if (present(diagnostic_charge)) surface_total_charge = sum(diagnostic_charge)
    surface_scale = max(abs(surface_total_charge), eps0*self%zero_plan%area_xy*field_scale, tiny(1.0_dp))
    local_charge_estimate = abs(eps0*self%zero_plan%area_xy*self%outer%interface_field)* &
                            gap/self%outer_options%debye_length
    self%diagnostics%split_periodic_active = .true.
    self%diagnostics%interface_sample_n = n
    self%diagnostics%interface_potential = self%outer%interface_potential
    self%diagnostics%interface_field = self%outer%interface_field
    self%diagnostics%eta_phi_kneq0 = max_phi/phi_scale
    self%diagnostics%eta_field_kneq0 = max_field/field_scale
    self%diagnostics%eta_gap = gap/self%outer_options%debye_length
    self%diagnostics%eta_local_charge = local_charge_estimate/surface_scale
    self%diagnostics%gauss_residual = self%gauss_residual
    self%diagnostics%outer_integrated_charge = &
      self%zero_plan%area_xy*(-eps0*self%outer%interface_field)
    self%diagnostics%applicable = &
      self%diagnostics%eta_phi_kneq0 <= self%periodic_options%interface_phi_tolerance .and. &
      self%diagnostics%eta_field_kneq0 <= self%periodic_options%interface_field_tolerance .and. &
      self%diagnostics%eta_gap <= self%outer_options%max_gap_ratio .and. &
      self%diagnostics%eta_local_charge <= self%outer_options%max_local_charge_ratio
    if (.not. self%diagnostics%applicable) then
      self%diagnostics%status = 'not_applicable'
      error stop 'split periodic scalar outer model is not applicable at the configured interface.'
    end if
    self%diagnostics%status = 'applicable'
  end subroutine update_interface_diagnostics

  subroutine init_cached_periodic_snapshot( &
    self, mesh, sim, field_config, periodic_config, panel_config, outer_config, kinetic_options &
    )
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(field_physics_config), intent(in) :: field_config
    type(periodic2_physics_config), intent(in) :: periodic_config
    type(panel_kernel_config), intent(in) :: panel_config
    type(outer_plasma_config), intent(in) :: outer_config
    type(kinetic_outer_plasma_options_type), intent(in), optional :: kinetic_options
    integer(i32) :: status
    character(len=256) :: message

    if (trim(lower_ascii(periodic_config%zero_mode_policy)) /= 'exclude_k0' .or. &
        .not. supported_lower_boundary(periodic_config%lower_boundary_model)) then
      error stop 'cached_kneq0 requires exclude_k0 and a supported lower boundary model.'
    end if
    if (trim(lower_ascii(outer_config%model)) /= 'none' .and. &
        trim(lower_ascii(outer_config%model)) /= 'kinetic_1d') then
      error stop 'cached_kneq0 supports none or kinetic_1d.'
    end if
    if (trim(lower_ascii(outer_config%model)) == 'kinetic_1d' .and. .not. present(kinetic_options)) then
      error stop 'cached kinetic_1d requires resolved kinetic options.'
    end if
    if (.not. sim%use_box .or. any(sim%bc_low(1:2) /= bc_periodic) .or. &
        any(sim%bc_high(1:2) /= bc_periodic) .or. sim%bc_low(3) == bc_periodic .or. &
        sim%bc_high(3) == bc_periodic) then
      error stop 'cached_kneq0 snapshot currently requires x/y periodic and z open.'
    end if
    call build_periodic_zero_mode_plan(mesh, product(sim%box_max(1:2) - sim%box_min(1:2)), self%zero_plan, status, message)
    if (status /= periodic_zero_mode_ok) error stop 'cached_kneq0 zero-mode plan: '//trim(message)
    self%periodic_options = periodic_config
    call refresh_periodic_zero_mode_state( &
      self%zero_plan, mesh%q_elem, zero_mode_bottom_field(self, mesh%q_elem), &
      self%zero_plan%break_z(1), 0.0_dp, self%zero_state &
      )
    call self%nonzero_solver%init(mesh, sim, field_config, periodic_config, panel_config)
    self%use_cached_kneq0 = .true.
    self%use_zero_mode = .true.
    self%outer_options = outer_config
    self%periodic_length = sim%box_max(1:2) - sim%box_min(1:2)
    self%periodic_origin = sim%box_min(1:2)
    self%reference_mode_layers = periodic_config%reference_mode_layers
    self%panel_quadrature_order = periodic_config%panel_quadrature_order
    self%mesh_top_z = max(maxval(mesh%v0(3, :)), max(maxval(mesh%v1(3, :)), maxval(mesh%v2(3, :))))
    if (present(kinetic_options)) self%kinetic_options = kinetic_options
    if (trim(lower_ascii(outer_config%model)) == 'kinetic_1d') self%use_outer_plasma = .true.
    self%diagnostics%split_periodic_active = .true.
    self%diagnostics%status = 'cached_kneq0'
    self%diagnostics%periodic_cache_hit = self%nonzero_solver%fmm_core_plan%periodic_cache_hit
    self%diagnostics%periodic_operator_build_count = &
      self%nonzero_solver%fmm_core_plan%periodic_operator_build_count
    self%diagnostics%periodic_cache_fingerprint = &
      self%nonzero_solver%fmm_core_plan%periodic_cache_fingerprint
    self%diagnostics%periodic_cache_path = self%nonzero_solver%fmm_core_plan%periodic_cache_path
  end subroutine init_cached_periodic_snapshot

  subroutine init_split_periodic_snapshot( &
    self, mesh, sim, periodic_config, outer_config, kinetic_options &
    )
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(periodic2_physics_config), intent(in) :: periodic_config
    type(outer_plasma_config), intent(in) :: outer_config
    type(kinetic_outer_plasma_options_type), intent(in), optional :: kinetic_options
    integer(i32) :: status
    character(len=256) :: message

    if (trim(lower_ascii(periodic_config%zero_mode_policy)) /= 'exclude_k0' .or. &
        .not. supported_lower_boundary(periodic_config%lower_boundary_model)) then
      error stop 'panel spectral reference requires exclude_k0 and a supported lower boundary model.'
    end if
    if (trim(lower_ascii(outer_config%model)) /= 'kinetic_1d') then
      error stop 'split periodic snapshot requires kinetic_1d.'
    end if
    if (trim(lower_ascii(outer_config%model)) == 'kinetic_1d' .and. .not. present(kinetic_options)) then
      error stop 'split periodic kinetic_1d requires resolved kinetic options.'
    end if
    if (.not. sim%use_box .or. any(sim%bc_low(1:2) /= bc_periodic) .or. &
        any(sim%bc_high(1:2) /= bc_periodic) .or. sim%bc_low(3) == bc_periodic .or. &
        sim%bc_high(3) == bc_periodic) then
      error stop 'split periodic snapshot currently requires x/y periodic and z open.'
    end if
    if (any(self%prescribed_e /= 0.0_dp)) then
      error stop 'split periodic outer snapshot currently requires sim.e0=0.'
    end if
    self%periodic_length = sim%box_max(1:2) - sim%box_min(1:2)
    self%periodic_origin = sim%box_min(1:2)
    if (any(self%periodic_length <= 0.0_dp)) error stop 'split periodic snapshot requires positive periods.'
    if (outer_config%interface_z <= &
        max(maxval(mesh%v0(3, :)), max(maxval(mesh%v1(3, :)), maxval(mesh%v2(3, :)))) .or. &
        abs(outer_config%interface_z - sim%box_max(3)) > &
        64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(sim%box_max(3)))) then
      error stop 'outer-plasma interface must equal the z-high box face and lie above the mesh.'
    end if
    if (outer_config%debye_length <= 0.0_dp .or. outer_config%thermal_voltage <= 0.0_dp) then
      error stop 'outer plasma requires positive Debye length and thermal voltage.'
    end if
    call build_periodic_zero_mode_plan( &
      mesh, product(self%periodic_length), self%zero_plan, status, message &
      )
    if (status /= periodic_zero_mode_ok) error stop 'zero-mode plan: '//trim(message)
    self%periodic_options = periodic_config
    call refresh_periodic_zero_mode_state( &
      self%zero_plan, mesh%q_elem, zero_mode_bottom_field(self, mesh%q_elem), &
      self%zero_plan%break_z(1), 0.0_dp, self%zero_state &
      )
    self%reference_mode_layers = periodic_config%reference_mode_layers
    self%panel_quadrature_order = periodic_config%panel_quadrature_order
    self%outer_options = outer_config
    if (present(kinetic_options)) self%kinetic_options = kinetic_options
    self%mesh_top_z = max(maxval(mesh%v0(3, :)), max(maxval(mesh%v1(3, :)), maxval(mesh%v2(3, :))))
    self%use_panel_spectral_reference = .true.
    self%use_zero_mode = .true.
    self%use_outer_plasma = .true.
  end subroutine init_split_periodic_snapshot

  subroutine refresh_cached_kinetic_outer( &
    self, mesh, update_outer, continuation_stage, continuation_batch &
    )
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    logical, intent(in), optional :: update_outer
    character(len=*), intent(in), optional :: continuation_stage
    integer(i32), intent(in), optional :: continuation_batch
    real(dp) :: raw_potential, interface_field
    logical :: refresh_outer

    refresh_outer = .true.
    if (present(update_outer)) refresh_outer = update_outer
    if (.not. self%outer%ready) refresh_outer = .true.
    call self%nonzero_solver%refresh(mesh)
    call refresh_periodic_zero_mode_state( &
      self%zero_plan, mesh%q_elem, zero_mode_bottom_field(self, mesh%q_elem), &
      self%zero_plan%break_z(1), 0.0_dp, self%zero_state &
      )
    call eval_periodic_zero_mode( &
      self%zero_plan, self%zero_state, self%outer_options%interface_z, zero_mode_trace_plus, &
      raw_potential, interface_field &
      )
    if (refresh_outer) then
      call solve_kinetic_collective( &
        self, interface_field, continuation_stage, continuation_batch &
        )
    end if
    call refresh_periodic_zero_mode_state( &
      self%zero_plan, mesh%q_elem, zero_mode_bottom_field(self, mesh%q_elem), self%zero_plan%break_z(1), &
      self%outer%interface_potential - raw_potential, self%zero_state &
      )
    self%gauss_residual = upper_flux_charge(self) + &
                          self%zero_plan%area_xy*(-eps0*self%outer%interface_field)
    call update_interface_diagnostics(self, mesh)
  end subroutine refresh_cached_kinetic_outer

  pure real(dp) function zero_mode_bottom_field(self, charge) result(field)
    class(electrostatic_snapshot_type), intent(in) :: self
    real(dp), intent(in) :: charge(:)

    select case (trim(lower_ascii(self%periodic_options%lower_boundary_model)))
    case ('e_bottom_zero')
      field = 0.0_dp
    case ('symmetric_vacuum')
      field = symmetric_vacuum_bottom_field(self%zero_plan, charge)
    case default
      error stop 'unsupported periodic2 lower boundary model.'
    end select
  end function zero_mode_bottom_field

  pure real(dp) function upper_flux_charge(self) result(charge)
    class(electrostatic_snapshot_type), intent(in) :: self

    charge = self%zero_state%total_charge + eps0*self%zero_plan%area_xy*self%zero_state%e_bottom
  end function upper_flux_charge

  pure logical function supported_lower_boundary(model) result(supported)
    character(len=*), intent(in) :: model

    select case (trim(lower_ascii(model)))
    case ('e_bottom_zero', 'symmetric_vacuum')
      supported = .true.
    case default
      supported = .false.
    end select
  end function supported_lower_boundary

  subroutine solve_kinetic_collective( &
    self, interface_field, continuation_stage, continuation_batch &
    )
    class(electrostatic_snapshot_type), intent(inout) :: self
    real(dp), intent(in) :: interface_field
    character(len=*), intent(in), optional :: continuation_stage
    integer(i32), intent(in), optional :: continuation_batch
    type(outer_plasma_state_type) :: solved
    type(zhao_continuation_diagnostics_type) :: zhao_diagnostics
    integer(i32) :: status_values(4), status
    real(dp) :: scalar_values(19)
    character(len=256) :: message

    self%kinetic_options%interface_field = interface_field
    status_values = 0_i32
    if (mpi_is_root(self%mpi)) then
      if (allocated(self%outer%potential)) then
        if (size(self%outer%potential) == self%kinetic_options%grid_points) then
          call solve_outer_plasma_kinetic( &
            self%kinetic_options, solved, status, message, initial_potential=self%outer%potential, &
            initial_state=self%outer, zhao_diagnostics=zhao_diagnostics &
            )
        else
          call solve_outer_plasma_kinetic( &
            self%kinetic_options, solved, status, message, zhao_diagnostics=zhao_diagnostics &
            )
        end if
      else
        call solve_outer_plasma_kinetic( &
          self%kinetic_options, solved, status, message, zhao_diagnostics=zhao_diagnostics &
          )
      end if
      status_values = [status, solved%profile_n, solved%nonlinear_iterations, int(iachar(solved%zhao_branch), i32)]
    end if
    call mpi_bcast_i32_array(self%mpi, status_values, 0_i32)
    status = status_values(1)
    if (status /= outer_plasma_ok) then
      if (mpi_is_root(self%mpi)) then
        if (trim(lower_ascii(self%kinetic_options%kinetic_closure)) == 'zhao_charge_driven' .and. &
            zhao_diagnostics%reason_code /= zhao_continuation_reason_none) then
          call write_zhao_continuation_diagnostics( &
            error_unit, zhao_diagnostics, continuation_stage, continuation_batch &
            )
        end if
        write (error_unit, '(a,i0,a,es24.16,2a)') 'kinetic outer-plasma solve failed without fallback: status=', status, &
          ' interface_field_V_m=', interface_field, ' message=', trim(message)
        flush (error_unit)
      end if
      call mpi_world_barrier(self%mpi)
      error stop 'kinetic outer-plasma solve failed without fallback.'
    end if
    if (.not. mpi_is_root(self%mpi)) then
      solved = outer_plasma_state_type()
      solved%model = 'kinetic_1d'
      solved%kinetic_closure = self%kinetic_options%kinetic_closure
      solved%zhao_branch = achar(status_values(4))
      solved%profile_n = status_values(2)
      solved%nonlinear_iterations = status_values(3)
      allocate (solved%z(solved%profile_n), solved%potential(solved%profile_n), &
                solved%field(solved%profile_n), solved%charge_density(solved%profile_n))
    end if
    if (mpi_is_root(self%mpi)) then
      scalar_values = [ &
                      solved%interface_potential, solved%infinity_potential, solved%debye_length, &
                      solved%interface_field, solved%nonlinear_residual, solved%integrated_charge_per_area, &
                      merge(1.0_dp, 0.0_dp, solved%ready), solved%electron_current_density, &
                      solved%ion_current_density, solved%photoelectron_current_density, solved%total_current_density, &
                      solved%zhao_phi0, solved%zhao_phi_minimum, solved%zhao_electron_density_infinity, &
                      solved%photoelectron_population_fraction, solved%photoelectron_column_per_area, &
                      solved%photoelectron_column_target_per_area, solved%photoelectron_column_residual_per_area, &
                      solved%photoelectron_source_scale &
                      ]
    end if
    call mpi_bcast_real_dp_array(self%mpi, scalar_values, 0_i32)
    call mpi_bcast_real_dp_array(self%mpi, solved%z, 0_i32)
    call mpi_bcast_real_dp_array(self%mpi, solved%potential, 0_i32)
    call mpi_bcast_real_dp_array(self%mpi, solved%field, 0_i32)
    call mpi_bcast_real_dp_array(self%mpi, solved%charge_density, 0_i32)
    solved%interface_z = self%outer_options%interface_z
    solved%z = solved%z + solved%interface_z
    solved%interface_potential = scalar_values(1)
    solved%infinity_potential = scalar_values(2)
    solved%debye_length = scalar_values(3)
    solved%interface_field = scalar_values(4)
    solved%nonlinear_residual = scalar_values(5)
    solved%integrated_charge_per_area = scalar_values(6)
    solved%ready = scalar_values(7) > 0.5_dp
    solved%electron_current_density = scalar_values(8)
    solved%ion_current_density = scalar_values(9)
    solved%photoelectron_current_density = scalar_values(10)
    solved%total_current_density = scalar_values(11)
    solved%zhao_phi0 = scalar_values(12)
    solved%zhao_phi_minimum = scalar_values(13)
    solved%zhao_electron_density_infinity = scalar_values(14)
    solved%photoelectron_population_fraction = scalar_values(15)
    solved%photoelectron_column_per_area = scalar_values(16)
    solved%photoelectron_column_target_per_area = scalar_values(17)
    solved%photoelectron_column_residual_per_area = scalar_values(18)
    solved%photoelectron_source_scale = scalar_values(19)
    solved%model = 'kinetic_1d'
    solved%kinetic_closure = self%kinetic_options%kinetic_closure
    solved%zhao_branch = achar(status_values(4))
    solved%applicability_status = outer_plasma_ok
    self%kinetic_options%photoelectron_source_scale = solved%photoelectron_source_scale
    self%kinetic_options%photoelectron_population_fraction = solved%photoelectron_population_fraction
    self%outer = solved
  end subroutine solve_kinetic_collective

end module bem_electrostatic_snapshot
