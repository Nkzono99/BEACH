!> 領域別の TOML 読み取り。意味検証と派生値の確定は preflight が担当する。
submodule(bem_app_config_parser) bem_app_config_parser_read_sim
  implicit none
contains

  module procedure apply_periodic2_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  authoring%periodic2%present = .true.
  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('nonzero_mode_backend')
      call get_toml_string( &
        table, keys(ikey), authoring%periodic2%nonzero_mode_backend, 'periodic2.nonzero_mode_backend' &
        )
      authoring%periodic2%nonzero_mode_backend = lower_ascii(trim(authoring%periodic2%nonzero_mode_backend))
      authoring%periodic2%has_nonzero_mode_backend = .true.
    case ('zero_mode_policy')
      call get_toml_string(table, keys(ikey), authoring%periodic2%zero_mode_policy, 'periodic2.zero_mode_policy')
      authoring%periodic2%zero_mode_policy = lower_ascii(trim(authoring%periodic2%zero_mode_policy))
      authoring%periodic2%has_zero_mode_policy = .true.
    case ('lower_boundary_model')
      call get_toml_string( &
        table, keys(ikey), authoring%periodic2%lower_boundary_model, 'periodic2.lower_boundary_model' &
        )
      authoring%periodic2%lower_boundary_model = lower_ascii(trim(authoring%periodic2%lower_boundary_model))
      authoring%periodic2%has_lower_boundary_model = .true.
    case ('reference_mode_layers')
      call get_toml_int( &
        table, keys(ikey), authoring%periodic2%reference_mode_layers, 'periodic2.reference_mode_layers' &
        )
    case ('panel_quadrature_order')
      call get_toml_int( &
        table, keys(ikey), authoring%periodic2%panel_quadrature_order, 'periodic2.panel_quadrature_order' &
        )
    case ('max_nonzero_mode_potential_step')
      call get_toml_real( &
        table, keys(ikey), authoring%periodic2%max_nonzero_mode_potential_step, &
        'periodic2.max_nonzero_mode_potential_step' &
        )
    case default
      error stop 'Unknown key in [periodic2]: '//trim(keys(ikey)%key)
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
      call get_toml_string(table, keys(ikey), field%mode, 'field_boundary.mode')
      field%mode = lower_ascii(trim(field%mode))
    case default
      error stop 'Unknown key in [field_boundary]: '//trim(keys(ikey)%key)
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
      call get_toml_particle_boundary_mode(table, keys(ikey), particles%low(1), 'particle_boundary.x_low', .false.)
    case ('x_high')
      call get_toml_particle_boundary_mode(table, keys(ikey), particles%high(1), 'particle_boundary.x_high', .false.)
    case ('y_low')
      call get_toml_particle_boundary_mode(table, keys(ikey), particles%low(2), 'particle_boundary.y_low', .false.)
    case ('y_high')
      call get_toml_particle_boundary_mode(table, keys(ikey), particles%high(2), 'particle_boundary.y_high', .false.)
    case ('z_low')
      call get_toml_particle_boundary_mode(table, keys(ikey), particles%low(3), 'particle_boundary.z_low', .false.)
    case ('z_high')
      call get_toml_particle_boundary_mode(table, keys(ikey), particles%high(3), 'particle_boundary.z_high', .false.)
    case ('ordinary_open_model')
      call get_toml_string( &
        table, keys(ikey), particles%ordinary_open_model, 'particle_boundary.ordinary_open_model' &
        )
      particles%ordinary_open_model = lower_ascii(trim(particles%ordinary_open_model))
    case default
      error stop 'Unknown key in [particle_boundary]: '//trim(keys(ikey)%key)
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
      call get_toml_string(table, keys(ikey), reservoir%inflow_model, 'reservoir.inflow_model')
      reservoir%inflow_model = lower_ascii(trim(reservoir%inflow_model))
    case ('phi_infty')
      call get_toml_real(table, keys(ikey), reservoir%phi_infty, 'reservoir.phi_infty')
    case ('face_potential_grid_n')
      call get_toml_int( &
        table, keys(ikey), reservoir%face_potential_grid_n, 'reservoir.face_potential_grid_n' &
        )
    case default
      error stop 'Unknown key in [reservoir]: '//trim(keys(ikey)%key)
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
      call get_toml_real(table, keys(ikey), cfg%sim%dt, 'sim.dt')
    case ('rng_seed')
      call get_toml_int(table, keys(ikey), cfg%sim%rng_seed, 'sim.rng_seed')
    case ('batch_count')
      call get_toml_int(table, keys(ikey), cfg%sim%batch_count, 'sim.batch_count')
    case ('batch_duration')
      call get_toml_real(table, keys(ikey), cfg%sim%batch_duration, 'sim.batch_duration')
      cfg%sim%has_batch_duration = .true.
    case ('batch_duration_step')
      call get_toml_real(table, keys(ikey), cfg%sim%batch_duration_step, 'sim.batch_duration_step')
      cfg%sim%has_batch_duration_step = .true.
    case ('max_step')
      call get_toml_int(table, keys(ikey), cfg%sim%max_step, 'sim.max_step')
    case ('tol_rel')
      call get_toml_real(table, keys(ikey), cfg%sim%tol_rel, 'sim.tol_rel')
    case ('q_floor')
      call get_toml_real(table, keys(ikey), cfg%sim%q_floor, 'sim.q_floor')
    case ('field_solver')
      call get_toml_string(table, keys(ikey), cfg%sim%field_solver, 'sim.field_solver')
      cfg%sim%field_solver = lower_ascii(trim(cfg%sim%field_solver))
    case ('field_normalization')
      call get_toml_string(table, keys(ikey), cfg%sim%field_normalization, 'sim.field_normalization')
      cfg%sim%field_normalization = lower_ascii(trim(cfg%sim%field_normalization))
    case ('field_length_scale')
      call get_toml_real(table, keys(ikey), cfg%sim%field_length_scale, 'sim.field_length_scale')
    case ('field_periodic_image_layers')
      call get_toml_int(table, keys(ikey), cfg%sim%field_periodic_image_layers, 'sim.field_periodic_image_layers')
    case ('field_periodic_far_correction')
      call get_toml_string( &
        table, keys(ikey), cfg%sim%field_periodic_far_correction, 'sim.field_periodic_far_correction' &
        )
      cfg%sim%field_periodic_far_correction = lower_ascii(trim(cfg%sim%field_periodic_far_correction))
    case ('field_periodic_ewald_alpha')
      call get_toml_real(table, keys(ikey), cfg%sim%field_periodic_ewald_alpha, 'sim.field_periodic_ewald_alpha')
    case ('field_periodic_ewald_layers')
      call get_toml_int(table, keys(ikey), cfg%sim%field_periodic_ewald_layers, 'sim.field_periodic_ewald_layers')
    case ('field_periodic_cache_dir')
      call get_toml_string(table, keys(ikey), cfg%sim%field_periodic_cache_dir, 'sim.field_periodic_cache_dir')
    case ('field_periodic_generation_tolerance')
      call get_toml_real( &
        table, keys(ikey), cfg%sim%field_periodic_generation_tolerance, &
        'sim.field_periodic_generation_tolerance' &
        )
    case ('tree_theta')
      call get_toml_real(table, keys(ikey), cfg%sim%tree_theta, 'sim.tree_theta')
      cfg%sim%has_tree_theta = .true.
    case ('tree_leaf_max')
      call get_toml_int(table, keys(ikey), cfg%sim%tree_leaf_max, 'sim.tree_leaf_max')
      cfg%sim%has_tree_leaf_max = .true.
    case ('tree_min_nelem')
      call get_toml_int(table, keys(ikey), cfg%sim%tree_min_nelem, 'sim.tree_min_nelem')
    case ('e0')
      call get_toml_real_array(table, keys(ikey), cfg%sim%e0, 'sim.e0')
      cfg%sim%has_e0_vector = .true.
    case ('e0_abs')
      call get_toml_real(table, keys(ikey), cfg%sim%e0_abs, 'sim.e0_abs')
      cfg%sim%has_e0_abs = .true.
    case ('e0_phi_xy_deg')
      call get_toml_real(table, keys(ikey), cfg%sim%e0_phi_xy_deg, 'sim.e0_phi_xy_deg')
      cfg%sim%has_e0_phi_xy_deg = .true.
    case ('e0_phi_z_deg')
      call get_toml_real(table, keys(ikey), cfg%sim%e0_phi_z_deg, 'sim.e0_phi_z_deg')
      cfg%sim%has_e0_phi_z_deg = .true.
    case ('b0')
      call get_toml_real_array(table, keys(ikey), cfg%sim%b0, 'sim.b0')
    case ('multiple_box_events_policy')
      call get_toml_string( &
        table, keys(ikey), cfg%sim%multiple_box_events_policy, 'sim.multiple_box_events_policy' &
        )
      cfg%sim%multiple_box_events_policy = lower_ascii(trim(cfg%sim%multiple_box_events_policy))
    case ('multiple_box_events_retry_backend')
      call get_toml_string( &
        table, keys(ikey), cfg%sim%multiple_box_events_retry_backend, &
        'sim.multiple_box_events_retry_backend' &
        )
      cfg%sim%multiple_box_events_retry_backend = lower_ascii(trim(cfg%sim%multiple_box_events_retry_backend))
    case ('multiple_box_events_soft_discard_count_grace')
      call get_toml_int( &
        table, keys(ikey), cfg%sim%multiple_box_events_soft_discard_count_grace, &
        'sim.multiple_box_events_soft_discard_count_grace' &
        )
    case ('multiple_box_events_soft_discard_fraction_limit')
      call get_toml_real( &
        table, keys(ikey), cfg%sim%multiple_box_events_soft_discard_fraction_limit, &
        'sim.multiple_box_events_soft_discard_fraction_limit' &
        )
    case ('multiple_box_events_soft_discard_abs_charge_limit')
      call get_toml_real( &
        table, keys(ikey), cfg%sim%multiple_box_events_soft_discard_abs_charge_limit, &
        'sim.multiple_box_events_soft_discard_abs_charge_limit' &
        )
    case ('raycast_max_bounce')
      call get_toml_int(table, keys(ikey), cfg%sim%raycast_max_bounce, 'sim.raycast_max_bounce')
    case default
      error stop 'Unknown key in [sim]: '//trim(keys(ikey)%key)
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
      call get_toml_logical(table, keys(ikey), cfg%write_output, 'output.write_files')
    case ('write_mesh_potential')
      call get_toml_logical(table, keys(ikey), cfg%write_mesh_potential, 'output.write_mesh_potential')
    case ('write_potential_history')
      call get_toml_logical(table, keys(ikey), cfg%write_potential_history, 'output.write_potential_history')
    case ('dir')
      call get_toml_string(table, keys(ikey), cfg%output_dir, 'output.dir')
    case ('history_stride')
      call get_toml_int(table, keys(ikey), cfg%history_stride, 'output.history_stride')
    case ('checkpoint_stride')
      call get_toml_int(table, keys(ikey), cfg%checkpoint_stride, 'output.checkpoint_stride')
    case ('resume')
      call get_toml_logical(table, keys(ikey), cfg%resume_output, 'output.resume')
    case ('restart_from')
      call get_toml_string(table, keys(ikey), cfg%output_restart_from, 'output.restart_from')
    case default
      error stop 'Unknown key in [output]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_output_toml_table

end submodule bem_app_config_parser_read_sim
