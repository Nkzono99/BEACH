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
    case ('secondary_emission_max_generations')
      call get_toml_int(table, keys(ikey), cfg%secondary_emission_max_generations, &
                        'particles.secondary_emission.max_generations')
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
    call require_toml_success(stat, species_config_path(ispec))
    if (.not. associated(child)) error stop species_config_path(ispec)//' must be a table.'
    cfg%particle_species(ispec) = species_from_defaults()
    cfg%particle_species(ispec)%enabled = .true.
    call apply_particles_species_toml_table(cfg%particle_species(ispec), child, authoring%particle_species(ispec), &
                                            species_config_path(ispec))
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
    case ('npcls_per_step', 'w_particle', 'target_macro_particles_per_batch', &
          'rays_per_batch', 'emit_current_density_a_m2', 'boundary_inflow')
      spec%has_secondary_forbidden_source_option = .true.
    end select
    select case (trim(k))
    case ('secondary_emission_enabled')
      call get_toml_logical(table, keys(ikey), spec%secondary_emission_enabled, &
                            species_path//'.surface_interaction.secondary_emission.enabled')
    case ('secondary_emission_model')
      call get_toml_string(table, keys(ikey), spec%secondary_emission_model, &
                           species_path//'.surface_interaction.secondary_emission.model')
      spec%secondary_emission_model = lower_ascii(trim(spec%secondary_emission_model))
    case ('secondary_emission_yield_max')
      call get_toml_real(table, keys(ikey), spec%secondary_emission_yield_max, &
                         species_path//'.surface_interaction.secondary_emission.yield_max')
    case ('secondary_emission_energy_max_ev')
      call get_toml_real(table, keys(ikey), spec%secondary_emission_energy_max_ev, &
                         species_path//'.surface_interaction.secondary_emission.energy_max_ev')
    case ('secondary_emission_const_yield')
      call get_toml_real(table, keys(ikey), spec%secondary_emission_const_yield, &
                         species_path//'.surface_interaction.secondary_emission.const_yield')
    case ('secondary_emission_target_species')
      call get_toml_string(table, keys(ikey), spec%secondary_emission_target_species, &
                           species_path//'.surface_interaction.secondary_emission.target_species')
    case ('species_key')
      call get_toml_string(table, keys(ikey), spec%species_key, species_path//'.species_key')
    case ('enabled')
      call get_toml_logical(table, keys(ikey), spec%enabled, species_path//'.enabled')
    case ('npcls_per_step')
      call get_toml_int(table, keys(ikey), spec%npcls_per_step, species_path//'.sampling.volume_macro_particles_per_batch')
    case ('source_mode')
      call get_toml_string(table, keys(ikey), spec%source_mode, species_path//'.source.mode')
      spec%source_mode = lower_ascii(trim(spec%source_mode))
    case ('number_density_cm3')
      call get_toml_real(table, keys(ikey), spec%number_density_cm3, species_path//'.distribution.number_density_cm3')
      spec%has_number_density_cm3 = .true.
    case ('number_density_m3')
      call get_toml_real(table, keys(ikey), spec%number_density_m3, species_path//'.distribution.number_density_m3')
      spec%has_number_density_m3 = .true.
    case ('q_particle')
      call get_toml_real(table, keys(ikey), spec%q_particle, species_path//'.charge_c')
    case ('m_particle')
      call get_toml_real(table, keys(ikey), spec%m_particle, species_path//'.mass_kg')
    case ('w_particle')
      call get_toml_real(table, keys(ikey), spec%w_particle, species_path//'.sampling.weight')
      spec%has_w_particle = .true.
    case ('target_macro_particles_per_batch')
      call get_toml_int( &
        table, keys(ikey), spec%target_macro_particles_per_batch, &
        species_path//'.sampling.target_macro_particles_per_batch' &
        )
      spec%has_target_macro_particles_per_batch = .true.
    case ('pos_low')
      call get_toml_real_array(table, keys(ikey), spec%pos_low, species_path//'.source.pos_low')
      auth%has_pos_low = .true.
    case ('pos_high')
      call get_toml_real_array(table, keys(ikey), spec%pos_high, species_path//'.source.pos_high')
      auth%has_pos_high = .true.
    case ('drift_velocity')
      call get_toml_real_array(table, keys(ikey), spec%drift_velocity, species_path//'.distribution.drift_velocity_m_s')
    case ('temperature_k')
      call get_toml_real(table, keys(ikey), spec%temperature_k, species_path//'.distribution.temperature_k')
    case ('temperature_ev')
      call get_toml_real(table, keys(ikey), spec%temperature_ev, species_path//'.distribution.temperature_ev')
      spec%has_temperature_ev = .true.
    case ('emit_current_density_a_m2')
      call get_toml_real(table, keys(ikey), spec%emit_current_density_a_m2, species_path//'.source.emit_current_density_a_m2')
    case ('rays_per_batch')
      call get_toml_int(table, keys(ikey), spec%rays_per_batch, species_path//'.sampling.rays_per_batch')
    case ('deposit_opposite_charge_on_emit')
      call get_toml_logical( &
        table, keys(ikey), spec%deposit_opposite_charge_on_emit, &
        species_path//'.source.deposit_opposite_charge_on_emit' &
        )
    case ('normal_drift_speed')
      call get_toml_real(table, keys(ikey), spec%normal_drift_speed, species_path//'.source.normal_drift_speed_m_s')
    case ('ray_direction')
      call get_toml_real_array(table, keys(ikey), spec%ray_direction, species_path//'.source.ray_direction')
      spec%has_ray_direction = .true.
    case ('inject_face')
      call get_toml_string(table, keys(ikey), spec%inject_face, species_path//'.source.inject_face')
      spec%inject_face = lower_ascii(trim(spec%inject_face))
    case ('boundary')
      nullify (child)
      call get_value(table, keys(ikey), child, requested=.false., stat=stat)
      call require_toml_success(stat, species_path//'.boundary')
      if (.not. associated(child)) error stop species_path//'.boundary must be a table.'
      call apply_species_boundary_toml_table(spec, child, species_path)
    case ('boundary_inflow')
      nullify (child)
      call get_value(table, keys(ikey), child, requested=.false., stat=stat)
      call require_toml_success(stat, species_path//'.inflow')
      if (.not. associated(child)) error stop species_path//'.inflow must be a table.'
      call apply_species_boundary_inflow_toml_table(spec, child, species_path)
    case ('surface_charge_closure')
      spec%has_surface_charge_closure = .true.
      call get_toml_string( &
        table, keys(ikey), spec%surface_charge_closure, species_path//'.charging.closure' &
        )
      spec%surface_charge_closure = lower_ascii(trim(spec%surface_charge_closure))
    case ('target_absorbed_current_a')
      call get_toml_real( &
        table, keys(ikey), spec%target_absorbed_current_a, species_path//'.charging.target_absorbed_current_a' &
        )
      spec%has_target_absorbed_current_a = .true.
    case ('target_emission_current_a')
      call get_toml_real( &
        table, keys(ikey), spec%target_emission_current_a, species_path//'.charging.target_emission_current_a' &
        )
      spec%has_target_emission_current_a = .true.
    case ('inject_region_mode')
      call get_toml_string(table, keys(ikey), auth%inject_region_mode, species_path//'.source.region_mode')
      auth%inject_region_mode = lower_ascii(trim(auth%inject_region_mode))
      auth%has_inject_region_mode = .true.
    case ('uv_low')
      call get_toml_real_array(table, keys(ikey), auth%uv_low, species_path//'.source.uv_low')
      auth%has_uv_low = .true.
    case ('uv_high')
      call get_toml_real_array(table, keys(ikey), auth%uv_high, species_path//'.source.uv_high')
      auth%has_uv_high = .true.
    case default
      error stop 'Unknown key in '//species_path//': '//trim(keys(ikey)%key)
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
        table, keys(ikey), spec%boundary_low(1), species_path//'.boundary.x_low', .true. &
        )
    case ('x_high')
      call get_toml_particle_boundary_mode( &
        table, keys(ikey), spec%boundary_high(1), species_path//'.boundary.x_high', .true. &
        )
    case ('y_low')
      call get_toml_particle_boundary_mode( &
        table, keys(ikey), spec%boundary_low(2), species_path//'.boundary.y_low', .true. &
        )
    case ('y_high')
      call get_toml_particle_boundary_mode( &
        table, keys(ikey), spec%boundary_high(2), species_path//'.boundary.y_high', .true. &
        )
    case ('z_low')
      call get_toml_particle_boundary_mode( &
        table, keys(ikey), spec%boundary_low(3), species_path//'.boundary.z_low', .true. &
        )
    case ('z_high')
      call get_toml_particle_boundary_mode( &
        table, keys(ikey), spec%boundary_high(3), species_path//'.boundary.z_high', .true. &
        )
    case default
      error stop 'Unknown key in ['//species_path//'.boundary]: '//trim(keys(ikey)%key)
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
        table, keys(ikey), spec%boundary_inflow_low(1), species_path//'.inflow.x_low' &
        )
    case ('x_high')
      call get_toml_boundary_inflow_mode( &
        table, keys(ikey), spec%boundary_inflow_high(1), species_path//'.inflow.x_high' &
        )
    case ('y_low')
      call get_toml_boundary_inflow_mode( &
        table, keys(ikey), spec%boundary_inflow_low(2), species_path//'.inflow.y_low' &
        )
    case ('y_high')
      call get_toml_boundary_inflow_mode( &
        table, keys(ikey), spec%boundary_inflow_high(2), species_path//'.inflow.y_high' &
        )
    case ('z_low')
      call get_toml_boundary_inflow_mode( &
        table, keys(ikey), spec%boundary_inflow_low(3), species_path//'.inflow.z_low' &
        )
    case ('z_high')
      call get_toml_boundary_inflow_mode( &
        table, keys(ikey), spec%boundary_inflow_high(3), species_path//'.inflow.z_high' &
        )
    case default
      error stop 'Unknown key in ['//species_path//'.inflow]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_species_boundary_inflow_toml_table

end submodule bem_app_config_parser_read_particles
