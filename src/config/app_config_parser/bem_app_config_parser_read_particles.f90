!> 領域別の TOML 読み取り。意味検証と派生値の確定は preflight が担当する。
submodule(bem_app_config_parser) bem_app_config_parser_read_particles
  implicit none
contains

  module procedure ensure_particle_species_capacity
  type(particle_species_spec), allocatable :: grown(:)
  integer :: old_capacity, new_capacity

  if (required_size <= 0) return
  if (allocated(cfg%particle_species)) then
    old_capacity = size(cfg%particle_species)
  else
    old_capacity = 0
  end if
  if (old_capacity >= required_size) return

  new_capacity = max(required_size, max(max_particle_species, max(1, 2*old_capacity)))
  allocate (grown(new_capacity))
  grown = particle_species_spec()
  if (old_capacity > 0) grown(1:old_capacity) = cfg%particle_species(1:old_capacity)
  call move_alloc(grown, cfg%particle_species)
  end procedure ensure_particle_species_capacity

  module procedure apply_particles_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('species')
      call read_particle_species_array(cfg, table, keys(ikey), authoring)
    case default
      error stop 'Unknown key in [particles]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_particles_toml_table

  module procedure read_particle_species_array
  type(config_toml_array), pointer :: array
  type(config_toml_table), pointer :: child
  integer :: ispec, n, stat

  nullify (array)
  call get_value(table, key, array, stat=stat)
  call require_toml_success(stat, 'particles.species')
  if (.not. associated(array)) error stop 'particles.species must be an array of tables.'

  n = toml_len(array)
  call ensure_particle_species_capacity(cfg, n)
  call ensure_authoring_particle_capacity(authoring, n)
  if (n > 0) cfg%particle_species(1:n) = particle_species_spec()
  if (n > 0) authoring%particle_species(1:n) = particle_authoring_spec()
  do ispec = 1, n
    nullify (child)
    call get_value(array, ispec, child, stat=stat)
    call require_toml_success(stat, 'particles.species entry')
    if (.not. associated(child)) error stop 'particles.species entries must be tables.'
    cfg%particle_species(ispec) = species_from_defaults()
    cfg%particle_species(ispec)%enabled = .true.
    call apply_particles_species_toml_table(cfg%particle_species(ispec), child, authoring%particle_species(ispec))
  end do
  cfg%n_particle_species = int(n, i32)
  end procedure read_particle_species_array

  module procedure apply_particles_species_toml_table
  type(toml_key), allocatable :: keys(:)
  type(config_toml_table), pointer :: child
  integer :: ikey, stat
  character(len=:), allocatable :: k

  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('species_key')
      call get_toml_string(table, keys(ikey), spec%species_key, 'particles.species.species_key')
    case ('enabled')
      call get_toml_logical(table, keys(ikey), spec%enabled, 'particles.species.enabled')
    case ('npcls_per_step')
      call get_toml_int(table, keys(ikey), spec%npcls_per_step, 'particles.species.npcls_per_step')
    case ('source_mode')
      call get_toml_string(table, keys(ikey), spec%source_mode, 'particles.species.source_mode')
      spec%source_mode = lower_ascii(trim(spec%source_mode))
    case ('number_density_cm3')
      call get_toml_real(table, keys(ikey), spec%number_density_cm3, 'particles.species.number_density_cm3')
      spec%has_number_density_cm3 = .true.
    case ('number_density_m3')
      call get_toml_real(table, keys(ikey), spec%number_density_m3, 'particles.species.number_density_m3')
    case ('q_particle')
      call get_toml_real(table, keys(ikey), spec%q_particle, 'particles.species.q_particle')
    case ('m_particle')
      call get_toml_real(table, keys(ikey), spec%m_particle, 'particles.species.m_particle')
    case ('w_particle')
      call get_toml_real(table, keys(ikey), spec%w_particle, 'particles.species.w_particle')
      spec%has_w_particle = .true.
    case ('target_macro_particles_per_batch')
      call get_toml_int( &
        table, keys(ikey), spec%target_macro_particles_per_batch, &
        'particles.species.target_macro_particles_per_batch' &
        )
      spec%has_target_macro_particles_per_batch = .true.
    case ('pos_low')
      call get_toml_real_array(table, keys(ikey), spec%pos_low, 'particles.species.pos_low')
      auth%has_pos_low = .true.
    case ('pos_high')
      call get_toml_real_array(table, keys(ikey), spec%pos_high, 'particles.species.pos_high')
      auth%has_pos_high = .true.
    case ('velocity_distribution')
      call get_toml_string(table, keys(ikey), spec%velocity_distribution, 'particles.species.velocity_distribution')
      spec%velocity_distribution = lower_ascii(trim(spec%velocity_distribution))
    case ('velocity_grid_path')
      call get_toml_string(table, keys(ikey), spec%velocity_grid_path, 'particles.species.velocity_grid_path')
    case ('velocity_grid_pdf_kind')
      call get_toml_string(table, keys(ikey), spec%velocity_grid_pdf_kind, 'particles.species.velocity_grid_pdf_kind')
      spec%velocity_grid_pdf_kind = lower_ascii(trim(spec%velocity_grid_pdf_kind))
    case ('velocity_grid_sampling')
      call get_toml_string(table, keys(ikey), spec%velocity_grid_sampling, 'particles.species.velocity_grid_sampling')
      spec%velocity_grid_sampling = lower_ascii(trim(spec%velocity_grid_sampling))
    case ('particle_flux_m2_s')
      call get_toml_real(table, keys(ikey), spec%particle_flux_m2_s, 'particles.species.particle_flux_m2_s')
    case ('current_density_a_m2')
      call get_toml_real(table, keys(ikey), spec%current_density_a_m2, 'particles.species.current_density_a_m2')
      spec%has_current_density_a_m2 = .true.
    case ('drift_velocity')
      call get_toml_real_array(table, keys(ikey), spec%drift_velocity, 'particles.species.drift_velocity')
    case ('temperature_k')
      call get_toml_real(table, keys(ikey), spec%temperature_k, 'particles.species.temperature_k')
    case ('temperature_ev')
      call get_toml_real(table, keys(ikey), spec%temperature_ev, 'particles.species.temperature_ev')
      spec%has_temperature_ev = .true.
    case ('emit_current_density_a_m2')
      call get_toml_real(table, keys(ikey), spec%emit_current_density_a_m2, 'particles.species.emit_current_density_a_m2')
    case ('rays_per_batch')
      call get_toml_int(table, keys(ikey), spec%rays_per_batch, 'particles.species.rays_per_batch')
    case ('deposit_opposite_charge_on_emit')
      call get_toml_logical( &
        table, keys(ikey), spec%deposit_opposite_charge_on_emit, &
        'particles.species.deposit_opposite_charge_on_emit' &
        )
    case ('normal_drift_speed')
      call get_toml_real(table, keys(ikey), spec%normal_drift_speed, 'particles.species.normal_drift_speed')
    case ('ray_direction')
      call get_toml_real_array(table, keys(ikey), spec%ray_direction, 'particles.species.ray_direction')
      spec%has_ray_direction = .true.
    case ('source_normal')
      call get_toml_real_array(table, keys(ikey), spec%source_normal, 'particles.species.source_normal')
    case ('inject_face')
      call get_toml_string(table, keys(ikey), spec%inject_face, 'particles.species.inject_face')
      spec%inject_face = lower_ascii(trim(spec%inject_face))
    case ('boundary')
      nullify (child)
      call get_value(table, keys(ikey), child, requested=.false., stat=stat)
      call require_toml_success(stat, 'particles.species.boundary')
      if (.not. associated(child)) error stop 'particles.species.boundary must be a table.'
      call apply_species_boundary_toml_table(spec, child)
    case ('boundary_inflow')
      nullify (child)
      call get_value(table, keys(ikey), child, requested=.false., stat=stat)
      call require_toml_success(stat, 'particles.species.boundary_inflow')
      if (.not. associated(child)) error stop 'particles.species.boundary_inflow must be a table.'
      call apply_species_boundary_inflow_toml_table(spec, child)
    case ('surface_charge_closure')
      call get_toml_string( &
        table, keys(ikey), spec%surface_charge_closure, 'particles.species.surface_charge_closure' &
        )
      spec%surface_charge_closure = lower_ascii(trim(spec%surface_charge_closure))
    case ('target_absorbed_current_a')
      call get_toml_real( &
        table, keys(ikey), spec%target_absorbed_current_a, 'particles.species.target_absorbed_current_a' &
        )
      spec%has_target_absorbed_current_a = .true.
    case ('target_emission_current_a')
      call get_toml_real( &
        table, keys(ikey), spec%target_emission_current_a, 'particles.species.target_emission_current_a' &
        )
      spec%has_target_emission_current_a = .true.
    case ('inject_region_mode')
      call get_toml_string(table, keys(ikey), auth%inject_region_mode, 'particles.species.inject_region_mode')
      auth%inject_region_mode = lower_ascii(trim(auth%inject_region_mode))
      auth%has_inject_region_mode = .true.
    case ('uv_low')
      call get_toml_real_array(table, keys(ikey), auth%uv_low, 'particles.species.uv_low')
      auth%has_uv_low = .true.
    case ('uv_high')
      call get_toml_real_array(table, keys(ikey), auth%uv_high, 'particles.species.uv_high')
      auth%has_uv_high = .true.
    case default
      error stop 'Unknown key in [[particles.species]]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_particles_species_toml_table

  module procedure apply_species_boundary_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('x_low')
      call get_toml_particle_boundary_mode( &
        table, keys(ikey), spec%boundary_low(1), 'particles.species.boundary.x_low', .true. &
        )
    case ('x_high')
      call get_toml_particle_boundary_mode( &
        table, keys(ikey), spec%boundary_high(1), 'particles.species.boundary.x_high', .true. &
        )
    case ('y_low')
      call get_toml_particle_boundary_mode( &
        table, keys(ikey), spec%boundary_low(2), 'particles.species.boundary.y_low', .true. &
        )
    case ('y_high')
      call get_toml_particle_boundary_mode( &
        table, keys(ikey), spec%boundary_high(2), 'particles.species.boundary.y_high', .true. &
        )
    case ('z_low')
      call get_toml_particle_boundary_mode( &
        table, keys(ikey), spec%boundary_low(3), 'particles.species.boundary.z_low', .true. &
        )
    case ('z_high')
      call get_toml_particle_boundary_mode( &
        table, keys(ikey), spec%boundary_high(3), 'particles.species.boundary.z_high', .true. &
        )
    case default
      error stop 'Unknown key in [particles.species.boundary]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_species_boundary_toml_table

  module procedure apply_species_boundary_inflow_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('x_low')
      call get_toml_boundary_inflow_mode( &
        table, keys(ikey), spec%boundary_inflow_low(1), 'particles.species.boundary_inflow.x_low' &
        )
    case ('x_high')
      call get_toml_boundary_inflow_mode( &
        table, keys(ikey), spec%boundary_inflow_high(1), 'particles.species.boundary_inflow.x_high' &
        )
    case ('y_low')
      call get_toml_boundary_inflow_mode( &
        table, keys(ikey), spec%boundary_inflow_low(2), 'particles.species.boundary_inflow.y_low' &
        )
    case ('y_high')
      call get_toml_boundary_inflow_mode( &
        table, keys(ikey), spec%boundary_inflow_high(2), 'particles.species.boundary_inflow.y_high' &
        )
    case ('z_low')
      call get_toml_boundary_inflow_mode( &
        table, keys(ikey), spec%boundary_inflow_low(3), 'particles.species.boundary_inflow.z_low' &
        )
    case ('z_high')
      call get_toml_boundary_inflow_mode( &
        table, keys(ikey), spec%boundary_inflow_high(3), 'particles.species.boundary_inflow.z_high' &
        )
    case default
      error stop 'Unknown key in [particles.species.boundary_inflow]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_species_boundary_inflow_toml_table

end submodule bem_app_config_parser_read_particles
