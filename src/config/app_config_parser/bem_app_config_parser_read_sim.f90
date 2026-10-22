!> 領域別の TOML 読み取り。意味検証と派生値の確定は preflight が担当する。
submodule(bem_app_config_parser) bem_app_config_parser_read_sim
  implicit none
contains

  module procedure apply_periodic2_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey, stat
  character(len=:), allocatable :: k, switch_text

  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    ! adaptive batch の上限も同じ内部 table に入るが split 設定ではない。それだけでは split backend を選ばない。
    if (trim(k) /= 'max_nonzero_mode_potential_step') authoring%periodic2%present = .true.
    select case (trim(k))
    case ('nonzero_mode_backend')
      call get_toml_string( &
        table, keys(ikey), authoring%periodic2%nonzero_mode_backend, 'fields.periodic.backend' &
        )
      authoring%periodic2%nonzero_mode_backend = lower_ascii(trim(authoring%periodic2%nonzero_mode_backend))
      authoring%periodic2%has_nonzero_mode_backend = .true.
    case ('zero_mode_policy')
      call get_toml_string(table, keys(ikey), authoring%periodic2%zero_mode_policy, 'fields.periodic.backend')
      authoring%periodic2%zero_mode_policy = lower_ascii(trim(authoring%periodic2%zero_mode_policy))
      authoring%periodic2%has_zero_mode_policy = .true.
    case ('lower_boundary_model')
      call get_toml_string( &
        table, keys(ikey), authoring%periodic2%lower_boundary_model, 'fields.periodic.lower_boundary_model' &
        )
      authoring%periodic2%lower_boundary_model = lower_ascii(trim(authoring%periodic2%lower_boundary_model))
      authoring%periodic2%has_lower_boundary_model = .true.
    case ('upper_fourier_switch_height')
      if (allocated(switch_text)) deallocate (switch_text)
      call get_value(table, keys(ikey), switch_text, stat=stat)
      if (allocated(switch_text)) then
        select case (trim(switch_text))
        case ('off', 'auto')
          authoring%periodic2%upper_fourier_switch = switch_text
        case default
          call stop_config_error('fields.periodic.upper_fourier.switch_height must be off, auto, or a height in metres.')
        end select
      else
        call get_toml_real(table, keys(ikey), authoring%periodic2%upper_fourier_switch_height, &
                           'fields.periodic.upper_fourier.switch_height')
        authoring%periodic2%upper_fourier_switch = 'height'
      end if
    case ('upper_fourier_mode_layers')
      call get_toml_int(table, keys(ikey), authoring%periodic2%upper_fourier_mode_layers, &
                        'fields.periodic.upper_fourier.mode_layers')
    case ('reference_mode_layers')
      call get_toml_int( &
        table, keys(ikey), authoring%periodic2%reference_mode_layers, 'fields.periodic.reference.mode_layers' &
        )
    case ('panel_quadrature_order')
      call get_toml_int( &
        table, keys(ikey), authoring%periodic2%panel_quadrature_order, 'fields.periodic.reference.panel_quadrature_order' &
        )
    case ('max_nonzero_mode_potential_step')
      call get_toml_real( &
        table, keys(ikey), authoring%periodic2%max_nonzero_mode_potential_step, &
        'run.batch.adaptive.max_nonzero_mode_potential_step_v' &
        )
    case default
      error stop 'Unknown key in [fields.periodic]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_periodic2_toml_table

  module procedure apply_domain_toml_table
  type(toml_key), allocatable :: keys(:)
  type(config_toml_array), pointer :: array
  integer :: ikey, iaxis, stat
  character(len=:), allocatable :: axis_name, k

  domain%present = .true.
  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('box_origin')
      call get_toml_real_array(table, keys(ikey), domain%box_origin, 'domain.box_origin')
      domain%has_box_origin = .true.
    case ('box_size')
      call get_toml_real_array(table, keys(ikey), domain%box_size, 'domain.box_size')
      domain%has_box_size = .true.
    case ('box_min')
      call get_toml_real_array(table, keys(ikey), domain%box_min, 'domain.box_min')
      domain%has_box_min = .true.
    case ('box_max')
      call get_toml_real_array(table, keys(ikey), domain%box_max, 'domain.box_max')
      domain%has_box_max = .true.
    case ('periodic_axes')
      nullify (array)
      call get_value(table, keys(ikey), array, stat=stat)
      call require_toml_success(stat, 'domain.periodic_axes')
      if (.not. associated(array)) error stop 'domain.periodic_axes must be an array of axis names.'
      if (toml_len(array) > 3) error stop 'domain.periodic_axes may contain at most x, y, and z.'
      domain%periodic_axis = .false.
      do iaxis = 1, toml_len(array)
        if (allocated(axis_name)) deallocate (axis_name)
        call get_value(array, iaxis, axis_name, stat=stat)
        call require_toml_success(stat, 'domain.periodic_axes')
        select case (trim(lower_ascii(axis_name)))
        case ('x')
          if (domain%periodic_axis(1)) error stop 'domain.periodic_axes contains duplicate "x".'
          domain%periodic_axis(1) = .true.
        case ('y')
          if (domain%periodic_axis(2)) error stop 'domain.periodic_axes contains duplicate "y".'
          domain%periodic_axis(2) = .true.
        case ('z')
          if (domain%periodic_axis(3)) error stop 'domain.periodic_axes contains duplicate "z".'
          domain%periodic_axis(3) = .true.
        case default
          error stop 'domain.periodic_axes entries must be "x", "y", or "z".'
        end select
      end do
    case default
      error stop 'Unknown key in [domain]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_domain_toml_table

  module procedure apply_field_boundary_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  field%present = .true.
  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('mode')
      call get_toml_string(table, keys(ikey), field%mode, 'fields.boundary')
      field%mode = lower_ascii(trim(field%mode))
    case default
      error stop 'Unknown key in [fields]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_field_boundary_toml_table

  module procedure apply_particle_boundary_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  particles%present = .true.
  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('x_low')
      call get_toml_particle_boundary_mode(table, keys(ikey), particles%low(1), 'particles.boundary.x_low', .false.)
    case ('x_high')
      call get_toml_particle_boundary_mode(table, keys(ikey), particles%high(1), 'particles.boundary.x_high', .false.)
    case ('y_low')
      call get_toml_particle_boundary_mode(table, keys(ikey), particles%low(2), 'particles.boundary.y_low', .false.)
    case ('y_high')
      call get_toml_particle_boundary_mode(table, keys(ikey), particles%high(2), 'particles.boundary.y_high', .false.)
    case ('z_low')
      call get_toml_particle_boundary_mode(table, keys(ikey), particles%low(3), 'particles.boundary.z_low', .false.)
    case ('z_high')
      call get_toml_particle_boundary_mode(table, keys(ikey), particles%high(3), 'particles.boundary.z_high', .false.)
    case default
      error stop 'Unknown key in [particles.boundary]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_particle_boundary_toml_table

  module procedure apply_reservoir_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  reservoir%present = .true.
  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('inflow_model')
      call get_toml_string(table, keys(ikey), reservoir%inflow_model, 'particles.reservoir.inflow_model')
      reservoir%inflow_model = lower_ascii(trim(reservoir%inflow_model))
    case ('phi_infty')
      call get_toml_real(table, keys(ikey), reservoir%phi_infty, 'particles.reservoir.phi_infty_v')
    case ('face_potential_grid_n')
      call get_toml_int( &
        table, keys(ikey), reservoir%face_potential_grid_n, 'particles.reservoir.face_potential_grid_n' &
        )
    case default
      error stop 'Unknown key in [particles.reservoir]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_reservoir_toml_table

  module procedure apply_sim_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('dt')
      call get_toml_real(table, keys(ikey), cfg%sim%dt, 'particles.tracking.dt_s')
    case ('rng_seed')
      call get_toml_int(table, keys(ikey), cfg%sim%rng_seed, 'run.rng_seed')
    case ('batch_count')
      call get_toml_int(table, keys(ikey), cfg%sim%batch_count, 'run.batch_count')
    case ('batch_duration')
      call get_toml_real(table, keys(ikey), cfg%sim%batch_duration, 'run.batch.duration_s')
      cfg%sim%has_batch_duration = .true.
    case ('batch_duration_step')
      call get_toml_real(table, keys(ikey), cfg%sim%batch_duration_step, 'run.batch.duration_steps')
      cfg%sim%has_batch_duration_step = .true.
    case ('max_step')
      call get_toml_int(table, keys(ikey), cfg%sim%max_step, 'particles.tracking.max_steps_per_particle')
    case ('substep_field_change_max')
      call get_toml_real(table, keys(ikey), cfg%sim%substep_field_change_max, &
                         'particles.tracking.substep_field_change_max')
    case ('substep_max_halvings')
      call get_toml_int(table, keys(ikey), cfg%sim%substep_max_halvings, 'particles.tracking.substep_max_halvings')
    case ('tol_rel')
      call get_toml_real(table, keys(ikey), cfg%sim%tol_rel, 'output.diagnostics.charge_rel_change_threshold')
    case ('q_floor')
      call get_toml_real(table, keys(ikey), cfg%sim%q_floor, 'output.diagnostics.charge_floor_c')
    case ('field_solver')
      call get_toml_string(table, keys(ikey), cfg%sim%field_solver, 'fields.solver.method')
      cfg%sim%field_solver = lower_ascii(trim(cfg%sim%field_solver))
    case ('field_normalization')
      call get_toml_string(table, keys(ikey), cfg%sim%field_normalization, 'fields.solver.normalization')
      cfg%sim%field_normalization = lower_ascii(trim(cfg%sim%field_normalization))
    case ('field_length_scale')
      call get_toml_real(table, keys(ikey), cfg%sim%field_length_scale, 'fields.solver.length_scale_m')
    case ('field_periodic_image_layers')
      call get_toml_int(table, keys(ikey), cfg%sim%field_periodic_image_layers, 'fields.periodic.image_layers')
    case ('field_periodic_far_correction')
      call get_toml_string( &
        table, keys(ikey), cfg%sim%field_periodic_far_correction, 'fields.periodic.backend' &
        )
      cfg%sim%field_periodic_far_correction = lower_ascii(trim(cfg%sim%field_periodic_far_correction))
    case ('field_periodic_ewald_alpha')
      call get_toml_real(table, keys(ikey), cfg%sim%field_periodic_ewald_alpha, 'fields.periodic.ewald_alpha')
    case ('field_periodic_ewald_layers')
      call get_toml_int(table, keys(ikey), cfg%sim%field_periodic_ewald_layers, 'fields.periodic.ewald_layers')
    case ('field_periodic_cache_dir')
      call get_toml_string(table, keys(ikey), cfg%sim%field_periodic_cache_dir, 'fields.periodic.cache_dir')
    case ('field_periodic_generation_tolerance')
      call get_toml_real( &
        table, keys(ikey), cfg%sim%field_periodic_generation_tolerance, &
        'fields.periodic.generation_tolerance' &
        )
    case ('tree_theta')
      call get_toml_real(table, keys(ikey), cfg%sim%tree_theta, 'fields.solver.tree.theta')
      cfg%sim%has_tree_theta = .true.
    case ('tree_leaf_max')
      call get_toml_int(table, keys(ikey), cfg%sim%tree_leaf_max, 'fields.solver.tree.leaf_max')
      cfg%sim%has_tree_leaf_max = .true.
    case ('tree_min_nelem')
      call get_toml_int(table, keys(ikey), cfg%sim%tree_min_nelem, 'fields.solver.auto_fmm_min_elements')
    case ('e0')
      call get_toml_real_array(table, keys(ikey), cfg%sim%e0, 'fields.external.electric_v_m')
      cfg%sim%has_e0_vector = .true.
    case ('e0_abs')
      call get_toml_real(table, keys(ikey), cfg%sim%e0_abs, 'fields.external.electric_magnitude_v_m')
      cfg%sim%has_e0_abs = .true.
    case ('e0_phi_xy_deg')
      call get_toml_real(table, keys(ikey), cfg%sim%e0_phi_xy_deg, 'fields.external.electric_azimuth_deg')
      cfg%sim%has_e0_phi_xy_deg = .true.
    case ('e0_phi_z_deg')
      call get_toml_real(table, keys(ikey), cfg%sim%e0_phi_z_deg, 'fields.external.electric_elevation_deg')
      cfg%sim%has_e0_phi_z_deg = .true.
    case ('b0')
      call get_toml_real_array(table, keys(ikey), cfg%sim%b0, 'fields.external.magnetic_t')
    case ('multiple_box_events_policy')
      call get_toml_string( &
        table, keys(ikey), cfg%sim%multiple_box_events_policy, 'particles.tracking.events.policy' &
        )
      cfg%sim%multiple_box_events_policy = lower_ascii(trim(cfg%sim%multiple_box_events_policy))
    case ('multiple_box_events_retry_backend')
      call get_toml_string( &
        table, keys(ikey), cfg%sim%multiple_box_events_retry_backend, &
        'particles.tracking.events.retry_backend' &
        )
      cfg%sim%multiple_box_events_retry_backend = lower_ascii(trim(cfg%sim%multiple_box_events_retry_backend))
    case ('multiple_box_events_soft_discard_count_grace')
      call get_toml_int( &
        table, keys(ikey), cfg%sim%multiple_box_events_soft_discard_count_grace, &
        'particles.tracking.events.discard_count_grace' &
        )
    case ('multiple_box_events_soft_discard_fraction_limit')
      call get_toml_real( &
        table, keys(ikey), cfg%sim%multiple_box_events_soft_discard_fraction_limit, &
        'particles.tracking.events.discard_fraction_limit' &
        )
    case ('multiple_box_events_soft_discard_abs_charge_limit')
      call get_toml_real( &
        table, keys(ikey), cfg%sim%multiple_box_events_soft_discard_abs_charge_limit, &
        'particles.tracking.events.discard_charge_warning_c' &
        )
    case ('raycast_max_bounce')
      call get_toml_int(table, keys(ikey), cfg%sim%raycast_max_bounce, 'particles.raycast.max_bounces')
    case default
      error stop 'Unknown key in [run], [particles.tracking], or [fields]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_sim_toml_table

  module procedure apply_output_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('write_files')
      call get_toml_logical(table, keys(ikey), cfg%write_output, 'output.enabled')
    case ('write_mesh_potential')
      call get_toml_logical(table, keys(ikey), cfg%write_mesh_potential, 'output.final.mesh_potential')
    case ('write_potential_history')
      call get_toml_logical(table, keys(ikey), cfg%write_potential_history, 'output.history.potential')
    case ('dir')
      call get_toml_string(table, keys(ikey), cfg%output_dir, 'output.dir')
    case ('history_stride')
      call get_toml_int(table, keys(ikey), cfg%history_stride, 'output.history.stride_batches')
    case ('checkpoint_stride')
      call get_toml_int(table, keys(ikey), cfg%checkpoint_stride, 'output.checkpoint.stride_batches')
    case ('resume')
      call get_toml_logical(table, keys(ikey), cfg%resume_output, 'run.restart')
    case ('restart_from')
      call get_toml_string(table, keys(ikey), cfg%output_restart_from, 'run.restart.from')
    case default
      error stop 'Unknown key in [output]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_output_toml_table

end submodule bem_app_config_parser_read_sim
