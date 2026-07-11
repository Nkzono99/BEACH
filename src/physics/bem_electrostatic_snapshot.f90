module bem_electrostatic_snapshot
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_types, only: mesh_type, sim_config, bc_periodic
  use bem_field_solver, only: field_solver_type
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config, &
                                      outer_plasma_config
  use bem_string_utils, only: lower_ascii
  use bem_periodic_zero_mode_plan, only: periodic_zero_mode_plan_type, periodic_zero_mode_state_type, &
                                         build_periodic_zero_mode_plan, refresh_periodic_zero_mode_state, &
                                         periodic_zero_mode_ok
  use bem_periodic_zero_mode_eval, only: eval_periodic_zero_mode, zero_mode_trace_plus
  use bem_coulomb_fmm_periodic_nonzero_reference, only: eval_periodic_nonzero_panel_reference
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok
  use bem_outer_plasma_linear, only: init_outer_plasma_linear, outer_plasma_integrated_charge_per_area
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
    character(len=64) :: status = 'legacy_or_not_applicable'
    integer(i32) :: last_outer_update_batch = -1_i32
    real(dp) :: max_outer_flight_time = 0.0_dp
    real(dp) :: max_frozen_field_ratio = 0.0_dp
    logical :: periodic_cache_hit = .false.
    integer(i32) :: periodic_operator_build_count = 0_i32
    character(len=128) :: periodic_cache_fingerprint = ''
    character(len=512) :: periodic_cache_path = ''
  end type electrostatic_diagnostics_type

  type, public :: electrostatic_restart_state_type
    logical :: outer_ready = .false.
    real(dp) :: outer_interface_potential = 0.0_dp
    integer(i32) :: last_outer_update_batch = -1_i32
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
    type(periodic2_physics_config) :: periodic_options
    type(electrostatic_diagnostics_type) :: diagnostics
    real(dp) :: mesh_top_z = 0.0_dp
    real(dp) :: gauss_residual = 0.0_dp
  contains
    procedure :: init => init_electrostatic_snapshot
    procedure :: refresh => refresh_electrostatic_snapshot
    procedure :: eval_local_e => eval_snapshot_local_e
    procedure :: eval_local_phi => eval_snapshot_local_phi
    procedure :: compute_mesh_potential => compute_snapshot_mesh_potential
    procedure :: get_diagnostics => get_snapshot_diagnostics
    procedure :: restore_outer_state => restore_snapshot_outer_state
    procedure :: export_restart_state => export_snapshot_restart_state
  end type electrostatic_snapshot_type

contains

  subroutine init_electrostatic_snapshot( &
    self, mesh, sim, field_config, periodic_config, panel_config, outer_config &
    )
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(field_physics_config), intent(in), optional :: field_config
    type(periodic2_physics_config), intent(in), optional :: periodic_config
    type(panel_kernel_config), intent(in), optional :: panel_config
    type(outer_plasma_config), intent(in), optional :: outer_config

    self%prescribed_e = sim%e0
    self%prescribed_phi_origin = 0.0_dp
    self%use_zero_mode = .false.
    self%use_outer_plasma = .false.
    self%use_panel_spectral_reference = .false.
    self%use_cached_kneq0 = .false.
    self%gauss_residual = 0.0_dp
    self%diagnostics = electrostatic_diagnostics_type()
    if (present(periodic_config) .and. present(panel_config) .and. present(outer_config)) then
      if (trim(lower_ascii(periodic_config%nonzero_mode_backend)) == 'panel_spectral_reference') then
        call init_split_periodic_snapshot(self, mesh, sim, periodic_config, panel_config, outer_config)
        return
      else if (trim(lower_ascii(periodic_config%nonzero_mode_backend)) == 'cached_kneq0') then
        if (.not. present(field_config)) error stop 'cached_kneq0 requires typed field configuration.'
        call init_cached_periodic_snapshot(self, mesh, sim, field_config, periodic_config, panel_config, outer_config)
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

  subroutine refresh_electrostatic_snapshot(self, mesh, update_outer)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    logical, intent(in), optional :: update_outer

    real(dp) :: interface_potential, interface_field, raw_potential, linearity_ratio
    integer(i32) :: status
    character(len=256) :: message
    logical :: refresh_outer

    if (self%use_cached_kneq0) then
      call self%nonzero_solver%refresh(mesh)
      call refresh_periodic_zero_mode_state( &
        self%zero_plan, mesh%q_elem, 0.0_dp, self%zero_plan%break_z(1), 0.0_dp, self%zero_state &
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
      self%zero_plan, mesh%q_elem, 0.0_dp, self%zero_plan%break_z(1), 0.0_dp, self%zero_state &
      )
    call eval_periodic_zero_mode( &
      self%zero_plan, self%zero_state, self%outer_options%interface_z, zero_mode_trace_plus, &
      raw_potential, interface_field &
      )
    if (refresh_outer) then
      interface_potential = self%outer_options%infinity_potential + self%outer_options%debye_length*interface_field
    else
      interface_potential = self%outer%interface_potential
    end if
    call refresh_periodic_zero_mode_state( &
      self%zero_plan, mesh%q_elem, 0.0_dp, self%zero_plan%break_z(1), &
      interface_potential - raw_potential, self%zero_state &
      )
    if (refresh_outer) then
      linearity_ratio = abs(interface_potential - self%outer_options%infinity_potential)/ &
                        self%outer_options%thermal_voltage
      call init_outer_plasma_linear( &
        self%outer_options%interface_z, interface_potential, self%outer_options%infinity_potential, &
        self%outer_options%debye_length, linearity_ratio, self%outer_options%max_linearity_ratio, &
        self%outer, status, message &
        )
      if (status /= outer_plasma_ok) error stop 'outer plasma refresh: '//trim(message)
    end if
    self%gauss_residual = sum(mesh%q_elem) + self%zero_plan%area_xy*outer_plasma_integrated_charge_per_area(self%outer)
    call update_interface_diagnostics(self, mesh)
  end subroutine refresh_electrostatic_snapshot

  subroutine eval_snapshot_local_e(self, mesh, position, electric_field)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: position(3)
    real(dp), intent(out) :: electric_field(3)

    real(dp) :: zero_potential, zero_field

    if (self%use_cached_kneq0) then
      call self%nonzero_solver%eval_e(mesh, position, electric_field)
      call eval_periodic_zero_mode( &
        self%zero_plan, self%zero_state, position(3), zero_mode_trace_plus, zero_potential, zero_field &
        )
      electric_field(3) = electric_field(3) + zero_field
    else if (self%use_panel_spectral_reference) then
      call eval_periodic_nonzero_panel_reference( &
        mesh, position, self%periodic_length(1), self%periodic_length(2), &
        self%reference_mode_layers, self%panel_quadrature_order, zero_potential, electric_field &
        )
      call eval_periodic_zero_mode( &
        self%zero_plan, self%zero_state, position(3), zero_mode_trace_plus, zero_potential, zero_field &
        )
      electric_field(3) = electric_field(3) + zero_field
    else
      call self%nonzero_solver%eval_e(mesh, position, electric_field)
    end if
    electric_field = electric_field + self%prescribed_e
  end subroutine eval_snapshot_local_e

  subroutine eval_snapshot_local_phi(self, mesh, sim, position, potential)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    real(dp), intent(in) :: position(3)
    real(dp), intent(out) :: potential

    real(dp) :: zero_potential, zero_field, electric_field_dummy(3)

    if (self%use_cached_kneq0) then
      call self%nonzero_solver%eval_potential(mesh, sim, position, potential)
      call eval_periodic_zero_mode( &
        self%zero_plan, self%zero_state, position(3), zero_mode_trace_plus, zero_potential, zero_field &
        )
      potential = potential + zero_potential
    else if (self%use_panel_spectral_reference) then
      call eval_periodic_nonzero_panel_reference( &
        mesh, position, self%periodic_length(1), self%periodic_length(2), &
        self%reference_mode_layers, self%panel_quadrature_order, potential, electric_field_dummy &
        )
      call eval_periodic_zero_mode( &
        self%zero_plan, self%zero_state, position(3), zero_mode_trace_plus, zero_potential, zero_field &
        )
      potential = potential + zero_potential
    else
      call self%nonzero_solver%eval_potential(mesh, sim, position, potential)
    end if
    potential = potential - dot_product(self%prescribed_e, position - self%prescribed_phi_origin)
  end subroutine eval_snapshot_local_phi

  subroutine compute_snapshot_mesh_potential(self, mesh, sim, potential)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    real(dp), intent(out) :: potential(:)
    integer(i32) :: element

    if (size(potential) /= mesh%nelem) error stop 'snapshot mesh-potential size mismatch.'
    if (self%use_panel_spectral_reference .or. self%use_cached_kneq0) then
      do element = 1, mesh%nelem
        call self%eval_local_phi(mesh, sim, mesh%centers(:, element), potential(element))
      end do
    else
      call self%nonzero_solver%compute_mesh_potential(mesh, sim, potential)
      do element = 1, mesh%nelem
        potential(element) = potential(element) - &
                             dot_product(self%prescribed_e, mesh%centers(:, element) - self%prescribed_phi_origin)
      end do
    end if
  end subroutine compute_snapshot_mesh_potential

  subroutine get_snapshot_diagnostics(self, diagnostics)
    class(electrostatic_snapshot_type), intent(in) :: self
    type(electrostatic_diagnostics_type), intent(out) :: diagnostics

    diagnostics = self%diagnostics
  end subroutine get_snapshot_diagnostics

  subroutine restore_snapshot_outer_state(self, state)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(electrostatic_restart_state_type), intent(in) :: state
    integer(i32) :: status
    character(len=256) :: message
    real(dp) :: linearity_ratio

    if (.not. state%outer_ready) return
    if (.not. self%use_outer_plasma) error stop 'checkpoint outer state is incompatible with the active model.'
    linearity_ratio = abs(state%outer_interface_potential - self%outer_options%infinity_potential)/ &
                      self%outer_options%thermal_voltage
    call init_outer_plasma_linear( &
      self%outer_options%interface_z, state%outer_interface_potential, self%outer_options%infinity_potential, &
      self%outer_options%debye_length, linearity_ratio, self%outer_options%max_linearity_ratio, &
      self%outer, status, message &
      )
    if (status /= outer_plasma_ok) error stop 'checkpoint outer state: '//trim(message)
  end subroutine restore_snapshot_outer_state

  subroutine export_snapshot_restart_state(self, last_outer_update_batch, state)
    class(electrostatic_snapshot_type), intent(in) :: self
    integer(i32), intent(in) :: last_outer_update_batch
    type(electrostatic_restart_state_type), intent(out) :: state

    state%outer_ready = self%outer%ready
    if (self%outer%ready) state%outer_interface_potential = self%outer%interface_potential
    state%last_outer_update_batch = last_outer_update_batch
  end subroutine export_snapshot_restart_state

  subroutine update_interface_diagnostics(self, mesh)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    integer(i32) :: ix, iy, n
    real(dp) :: position(3), potential, field(3), max_phi, max_field
    real(dp) :: phi_scale, field_scale, gap, surface_scale, local_charge_estimate

    n = self%periodic_options%interface_sample_n
    max_phi = 0.0_dp
    max_field = 0.0_dp
    do iy = 0_i32, n - 1_i32
      do ix = 0_i32, n - 1_i32
        position = [ &
                   self%periodic_origin(1) + (real(ix, dp) + 0.5_dp)*self%periodic_length(1)/real(n, dp), &
                   self%periodic_origin(2) + (real(iy, dp) + 0.5_dp)*self%periodic_length(2)/real(n, dp), &
                   self%outer_options%interface_z &
                   ]
        call eval_periodic_nonzero_panel_reference( &
          mesh, position, self%periodic_length(1), self%periodic_length(2), &
          self%reference_mode_layers, self%panel_quadrature_order, potential, field &
          )
        max_phi = max(max_phi, abs(potential))
        max_field = max(max_field, sqrt(sum(field*field)))
      end do
    end do
    phi_scale = max(abs(self%outer%interface_potential), self%outer_options%thermal_voltage, tiny(1.0_dp))
    field_scale = max(abs(self%outer%interface_field), &
                      self%outer_options%thermal_voltage/minval(self%periodic_length), tiny(1.0_dp))
    gap = self%outer_options%interface_z - self%mesh_top_z
    surface_scale = max(abs(sum(mesh%q_elem)), eps0*self%zero_plan%area_xy*field_scale, tiny(1.0_dp))
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
      self%zero_plan%area_xy*outer_plasma_integrated_charge_per_area(self%outer)
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

  subroutine init_cached_periodic_snapshot(self, mesh, sim, field_config, periodic_config, panel_config, outer_config)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(field_physics_config), intent(in) :: field_config
    type(periodic2_physics_config), intent(in) :: periodic_config
    type(panel_kernel_config), intent(in) :: panel_config
    type(outer_plasma_config), intent(in) :: outer_config
    integer(i32) :: status
    character(len=256) :: message

    if (trim(lower_ascii(periodic_config%zero_mode_policy)) /= 'exclude_k0' .or. &
        trim(lower_ascii(periodic_config%lower_boundary_model)) /= 'e_bottom_zero') then
      error stop 'cached_kneq0 requires exclude_k0 and e_bottom_zero.'
    end if
    if (trim(lower_ascii(outer_config%model)) /= 'none') then
      error stop 'cached_kneq0 currently requires outer_plasma.model=none.'
    end if
    if (.not. sim%use_box .or. any(sim%bc_low(1:2) /= bc_periodic) .or. &
        any(sim%bc_high(1:2) /= bc_periodic) .or. sim%bc_low(3) == bc_periodic .or. &
        sim%bc_high(3) == bc_periodic) then
      error stop 'cached_kneq0 snapshot currently requires x/y periodic and z open.'
    end if
    call build_periodic_zero_mode_plan(mesh, product(sim%box_max(1:2) - sim%box_min(1:2)), self%zero_plan, status, message)
    if (status /= periodic_zero_mode_ok) error stop 'cached_kneq0 zero-mode plan: '//trim(message)
    call refresh_periodic_zero_mode_state( &
      self%zero_plan, mesh%q_elem, 0.0_dp, self%zero_plan%break_z(1), 0.0_dp, self%zero_state &
      )
    call self%nonzero_solver%init(mesh, sim, field_config, periodic_config, panel_config)
    self%use_cached_kneq0 = .true.
    self%use_zero_mode = .true.
    self%periodic_options = periodic_config
    self%diagnostics%split_periodic_active = .true.
    self%diagnostics%status = 'cached_kneq0'
    self%diagnostics%periodic_cache_hit = self%nonzero_solver%fmm_core_plan%periodic_cache_hit
    self%diagnostics%periodic_operator_build_count = &
      self%nonzero_solver%fmm_core_plan%periodic_operator_build_count
    self%diagnostics%periodic_cache_fingerprint = &
      self%nonzero_solver%fmm_core_plan%periodic_cache_fingerprint
    self%diagnostics%periodic_cache_path = self%nonzero_solver%fmm_core_plan%periodic_cache_path
  end subroutine init_cached_periodic_snapshot

  subroutine init_split_periodic_snapshot(self, mesh, sim, periodic_config, panel_config, outer_config)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(periodic2_physics_config), intent(in) :: periodic_config
    type(panel_kernel_config), intent(in) :: panel_config
    type(outer_plasma_config), intent(in) :: outer_config
    integer(i32) :: status
    character(len=256) :: message

    if (trim(lower_ascii(periodic_config%zero_mode_policy)) /= 'exclude_k0' .or. &
        trim(lower_ascii(periodic_config%lower_boundary_model)) /= 'e_bottom_zero') then
      error stop 'panel spectral reference requires exclude_k0 and e_bottom_zero.'
    end if
    if (trim(lower_ascii(panel_config%source_model)) /= 'triangle_p0') then
      error stop 'panel spectral reference requires triangle_p0 sources.'
    end if
    if (trim(lower_ascii(outer_config%model)) /= 'linear_debye') then
      error stop 'split periodic snapshot requires outer_plasma.model=linear_debye.'
    end if
    if (.not. sim%use_box .or. any(sim%bc_low(1:2) /= bc_periodic) .or. &
        any(sim%bc_high(1:2) /= bc_periodic) .or. sim%bc_low(3) == bc_periodic .or. &
        sim%bc_high(3) == bc_periodic) then
      error stop 'split periodic snapshot currently requires x/y periodic and z open.'
    end if
    if (any(self%prescribed_e /= 0.0_dp)) then
      error stop 'split periodic linear outer snapshot currently requires sim.e0=0.'
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
      error stop 'linear outer plasma requires positive Debye length and thermal voltage.'
    end if
    call build_periodic_zero_mode_plan( &
      mesh, product(self%periodic_length), self%zero_plan, status, message &
      )
    if (status /= periodic_zero_mode_ok) error stop 'zero-mode plan: '//trim(message)
    self%reference_mode_layers = periodic_config%reference_mode_layers
    self%panel_quadrature_order = periodic_config%panel_quadrature_order
    self%periodic_options = periodic_config
    self%outer_options = outer_config
    self%mesh_top_z = max(maxval(mesh%v0(3, :)), max(maxval(mesh%v1(3, :)), maxval(mesh%v2(3, :))))
    self%use_panel_spectral_reference = .true.
    self%use_zero_mode = .true.
    self%use_outer_plasma = .true.
  end subroutine init_split_periodic_snapshot

end module bem_electrostatic_snapshot
