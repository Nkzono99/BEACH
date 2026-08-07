!> `bem_field_solver` の初期化・設定補助手続きを実装する submodule。
submodule(bem_field_solver) bem_field_solver_config
  use bem_coulomb_fmm_core, only: build_panel_plan, update_state, destroy_plan, destroy_state
  use bem_physics_config_types, only: validate_phase1_panel_config, physics_config_ok
  implicit none
contains

  !> 設定値から direct/treecode/fmm の実行モードを確定し、必要なら木構造を初期化する。
  module procedure init_field_solver
  character(len=16) :: requested_mode, field_bc_mode
  integer(i32) :: axis, n_periodic
  real(dp) :: span
  real(dp), allocatable :: panel_v0(:, :), panel_v1(:, :), panel_v2(:, :)
  integer(i32) :: panel_status
  character(len=256) :: panel_message

  self%fmm_core_state = fmm_state_type()
  call destroy_plan(self%fmm_core_plan)
  call destroy_state(self%fmm_core_state)
  self%fmm_use_core = .false.
  self%fmm_core_ready = .false.
  call reset_tree_storage(self)

  self%field_normalization = lower_ascii(trim(sim%field_normalization))
  self%field_length_scale = resolve_field_length_scale(mesh, sim)
  self%field_origin = resolve_field_origin(mesh, sim, self%field_normalization)
  self%field_inv_length_scale = 1.0d0/self%field_length_scale
  self%field_output_scale = k_coulomb*self%field_inv_length_scale*self%field_inv_length_scale
  self%potential_output_scale = k_coulomb*self%field_inv_length_scale
  self%min_nelem = sim%tree_min_nelem
  self%nelem = mesh%nelem
  self%use_periodic2 = .false.
  self%periodic_axes = 0_i32
  self%periodic_len = 0.0d0
  self%periodic_image_layers = max(0_i32, sim%field_periodic_image_layers)
  self%periodic_far_correction = lower_ascii(trim(sim%field_periodic_far_correction))
  self%periodic_ewald_alpha = sim%field_periodic_ewald_alpha
  self%periodic_ewald_layers = max(0_i32, sim%field_periodic_ewald_layers)
  self%periodic_cache_dir = trim(sim%field_periodic_cache_dir)
  self%periodic_generation_tolerance = sim%field_periodic_generation_tolerance
  self%target_box_min = 0.0d0
  self%target_box_max = 0.0d0

  requested_mode = lower_ascii(trim(sim%field_solver))
  field_bc_mode = lower_ascii(trim(sim%field_bc_mode))
  if (present(periodic_config)) continue
  if (present(panel_config)) then
    call validate_phase1_panel_config(sim, panel_config, panel_status, panel_message)
    if (panel_status /= physics_config_ok) error stop 'field solver panel config: '//trim(panel_message)
  end if
  if (present(field_config)) then
    if (trim(lower_ascii(field_config%backend)) /= trim(requested_mode) .or. &
        trim(lower_ascii(field_config%normalization)) /= trim(self%field_normalization)) then
      error stop 'typed field config does not match sim field config.'
    end if
  end if
  if (mesh%nelem > 0_i32) then
    if (any(mesh%panel_area <= 0.0_dp)) then
      error stop 'triangle_p0 requires finite, non-degenerate triangles.'
    end if
    if (any(abs(mesh%elem_vacuum_sign) /= 1_i32)) then
      error stop 'triangle_p0 requires resolved element vacuum sides.'
    end if
  end if
  self%field_bc_mode = field_bc_mode
  select case (trim(field_bc_mode))
  case ('free')
    self%periodic_far_correction = 'none'
    continue
  case ('periodic2')
    if (trim(requested_mode) /= 'fmm') then
      error stop 'sim.field_bc_mode must be "free" unless sim.field_solver="fmm".'
    end if
    if (.not. sim%use_box) then
      error stop 'sim.field_bc_mode="periodic2" requires sim.use_box=true.'
    end if
    n_periodic = 0_i32
    do axis = 1_i32, 3_i32
      if ((sim%bc_low(axis) == bc_periodic) .neqv. (sim%bc_high(axis) == bc_periodic)) then
        error stop 'periodic2 requires bc_low(axis)=bc_high(axis)=periodic for periodic axes.'
      end if
      if (sim%bc_low(axis) == bc_periodic) then
        n_periodic = n_periodic + 1_i32
        if (n_periodic <= 2_i32) self%periodic_axes(n_periodic) = axis
      end if
    end do
    if (n_periodic /= 2_i32) then
      error stop 'sim.field_bc_mode="periodic2" requires exactly two periodic axes.'
    end if
    do axis = 1_i32, 2_i32
      span = sim%box_max(self%periodic_axes(axis)) - sim%box_min(self%periodic_axes(axis))
      if (span <= 0.0d0) error stop 'periodic2 requires positive box length on periodic axes.'
      self%periodic_len(axis) = span
    end do
    self%target_box_min = sim%box_min
    self%target_box_max = sim%box_max
    select case (trim(self%periodic_far_correction))
    case ('auto')
      self%periodic_far_correction = 'none'
    case ('none')
      continue
    case ('m2l_root_oracle')
      error stop 'periodic2 far correction "m2l_root_oracle" was removed; use "cached_kneq0".'
    case ('cached_kneq0')
      continue
    case default
      error stop 'periodic2 far correction supports "auto", "none", or "cached_kneq0" only.'
    end select
    self%use_periodic2 = .true.
  case default
    error stop 'Unknown sim.field_bc_mode in field solver init.'
  end select

  self%mode = resolve_field_solver_mode(mesh%nelem, sim)

  if (trim(self%mode) == 'fmm' .and. sim%use_box) then
    if (any(sim%box_max <= sim%box_min)) then
      error stop 'sim.use_box=true requires positive box extents for fmm dual-target.'
    end if
    self%target_box_min = sim%box_min
    self%target_box_max = sim%box_max
  end if

  call resolve_field_solver_tree_params(mesh%nelem, sim, self%theta, self%leaf_max)

  self%fmm_core_options = fmm_options_type()
  self%fmm_core_options%theta = self%theta
  self%fmm_core_options%leaf_max = self%leaf_max
  self%fmm_core_options%order = field_solver_fmm_expansion_order
  self%fmm_core_options%softening = 0.0_dp
  self%fmm_core_options%use_periodic2 = self%use_periodic2
  self%fmm_core_options%periodic_far_correction = self%periodic_far_correction
  self%fmm_core_options%periodic_axes = self%periodic_axes
  self%fmm_core_options%periodic_len = self%periodic_len*self%field_inv_length_scale
  self%fmm_core_options%periodic_image_layers = self%periodic_image_layers
  self%fmm_core_options%periodic_ewald_alpha = self%periodic_ewald_alpha*self%field_length_scale
  self%fmm_core_options%periodic_ewald_layers = self%periodic_ewald_layers
  self%fmm_core_options%periodic_cache_dir = self%periodic_cache_dir
  self%fmm_core_options%periodic_generation_tolerance = self%periodic_generation_tolerance
  self%fmm_core_options%target_box_min = (self%target_box_min - self%field_origin)*self%field_inv_length_scale
  self%fmm_core_options%target_box_max = (self%target_box_max - self%field_origin)*self%field_inv_length_scale

  if (trim(self%mode) == 'fmm') then
    select case (trim(self%periodic_far_correction))
    case ('auto', 'none', 'cached_kneq0')
      continue
    case default
      error stop 'FMM core received an unsupported periodic far correction.'
    end select
    self%fmm_use_core = .true.
    if (mesh%nelem > 0_i32) then
      allocate (panel_v0(3, mesh%nelem), panel_v1(3, mesh%nelem), panel_v2(3, mesh%nelem))
      panel_v0 = (mesh%v0 - spread(self%field_origin, 2, mesh%nelem))*self%field_inv_length_scale
      panel_v1 = (mesh%v1 - spread(self%field_origin, 2, mesh%nelem))*self%field_inv_length_scale
      panel_v2 = (mesh%v2 - spread(self%field_origin, 2, mesh%nelem))*self%field_inv_length_scale
      call build_panel_plan(self%fmm_core_plan, panel_v0, panel_v1, panel_v2, self%fmm_core_options)
      deallocate (panel_v0, panel_v1, panel_v2)
      call update_state(self%fmm_core_plan, self%fmm_core_state, mesh%q_elem)
      self%fmm_core_ready = self%fmm_core_plan%built .and. self%fmm_core_state%ready
      call sync_core_plan_view(self)
    end if
    return
  end if

  if (trim(self%mode) == 'treecode' .and. mesh%nelem > 0_i32) then
    call build_tree_topology(self, mesh)
    call refresh_field_solver(self, mesh)
  else
    call reset_tree_storage(self)
    self%tree_ready = .false.
  end if
  end procedure init_field_solver

  !> FMM core 内部で使う長さ正規化スケール `L0` [m] を解決する。
  real(dp) function resolve_field_length_scale(mesh, sim) result(length_scale)
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    real(dp) :: bbox_min(3), bbox_max(3), span(3)
    character(len=16) :: mode

    mode = lower_ascii(trim(sim%field_normalization))
    select case (trim(mode))
    case ('si')
      length_scale = 1.0d0
    case ('length')
      length_scale = sim%field_length_scale
    case ('box')
      if (.not. sim%use_box) error stop 'sim.field_normalization="box" requires sim.use_box=true.'
      span = sim%box_max - sim%box_min
      if (any(span <= 0.0d0)) error stop 'sim.field_normalization="box" requires positive sim.box extents.'
      length_scale = maxval(span)
    case ('mesh')
      if (mesh%nelem <= 0_i32) then
        length_scale = sim%field_length_scale
      else
        bbox_min = min(minval(mesh%v0, dim=2), min(minval(mesh%v1, dim=2), minval(mesh%v2, dim=2)))
        bbox_max = max(maxval(mesh%v0, dim=2), max(maxval(mesh%v1, dim=2), maxval(mesh%v2, dim=2)))
        span = bbox_max - bbox_min
        length_scale = maxval(span)
      end if
    case default
      error stop 'Unknown sim.field_normalization in field solver init.'
    end select

    if (length_scale <= 0.0d0) error stop 'Resolved field normalization length must be > 0.'
  end function resolve_field_length_scale

  !> FMM core 内部で使う座標原点 [m] を解決する。
  function resolve_field_origin(mesh, sim, normalization) result(origin)
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    character(len=*), intent(in) :: normalization
    real(dp) :: origin(3)
    character(len=16) :: mode

    mode = lower_ascii(trim(normalization))
    select case (trim(mode))
    case ('box')
      origin = sim%box_min
    case ('mesh')
      if (mesh%nelem <= 0_i32) then
        origin = 0.0d0
      else
        origin = min(minval(mesh%v0, dim=2), min(minval(mesh%v1, dim=2), minval(mesh%v2, dim=2)))
      end if
    case default
      origin = 0.0d0
    end select
  end function resolve_field_origin

  !> 要素数レンジに応じた treecode/FMM の `theta` と `leaf_max` 推奨値を返す。
  module procedure estimate_auto_tree_params
  if (nelem < 1500_i32) then
    theta = 0.40d0
    leaf_max = 12_i32
  else if (nelem < 10000_i32) then
    theta = 0.50d0
    leaf_max = 16_i32
  else if (nelem < 50000_i32) then
    theta = 0.58d0
    leaf_max = 20_i32
  else
    theta = 0.65d0
    leaf_max = 24_i32
  end if
  end procedure estimate_auto_tree_params

  module procedure resolve_field_solver_tree_params
  character(len=16) :: requested_mode
  logical :: uses_tree

  requested_mode = lower_ascii(trim(sim%field_solver))
  uses_tree = trim(requested_mode) == 'treecode' .or. trim(requested_mode) == 'fmm' .or. &
              (trim(requested_mode) == 'auto' .and. nelem >= sim%tree_min_nelem)
  theta = sim%tree_theta
  leaf_max = sim%tree_leaf_max
  if (.not. uses_tree) return

  call estimate_auto_tree_params(nelem, theta, leaf_max)
  if (sim%has_tree_theta) theta = sim%tree_theta
  if (sim%has_tree_leaf_max) leaf_max = sim%tree_leaf_max
  end procedure resolve_field_solver_tree_params

  module procedure resolve_field_solver_mode
  character(len=16) :: requested_mode

  requested_mode = lower_ascii(trim(sim%field_solver))
  select case (trim(requested_mode))
  case ('direct', 'treecode', 'fmm')
    mode = requested_mode
  case ('auto')
    if (nelem >= sim%tree_min_nelem) then
      mode = 'fmm'
    else
      mode = 'direct'
    end if
  case default
    error stop 'Unknown sim.field_solver in field solver mode resolution.'
  end select
  end procedure resolve_field_solver_mode

end submodule bem_field_solver_config
