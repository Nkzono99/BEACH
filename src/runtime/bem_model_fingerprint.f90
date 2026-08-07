!> Restart compatibility fingerprints for the ordered physical model contract.
module bem_model_fingerprint
  use bem_kinds, only: dp, i32, i64
  use bem_types, only: mesh_type, bc_redistributed_reflect
  use bem_app_config_types, only: app_config, particle_species_spec, particle_bc_inherit
  implicit none
  private

  integer(i64), parameter :: hash_modulus = 2147483647_i64
  integer(i64), parameter :: hash_multiplier_a = 65599_i64
  integer(i64), parameter :: hash_multiplier_b = 131071_i64

  type :: hash_state
    integer(i64) :: a = 146959810_i64
    integer(i64) :: b = 109951162_i64
  end type hash_state

  public :: model_fingerprint, mesh_fingerprint, species_fingerprint

contains

  function model_fingerprint(cfg) result(fingerprint)
    type(app_config), intent(in) :: cfg
    character(len=16) :: fingerprint
    type(hash_state) :: hash

    call feed_string(hash, 'same_time_midpoint_boris')
    call feed_string(hash, cfg%field%backend)
    call feed_string(hash, cfg%field%normalization)
    call feed_string(hash, cfg%periodic2%nonzero_mode_backend)
    call feed_string(hash, cfg%periodic2%zero_mode_policy)
    call feed_string(hash, cfg%periodic2%lower_boundary_model)
    call feed_integer(hash, cfg%periodic2%reference_mode_layers)
    call feed_integer(hash, cfg%periodic2%panel_quadrature_order)
    if (cfg%periodic2%max_nonzero_mode_potential_step > 0.0_dp) then
      call feed_string(hash, 'adaptive_nonzero_mode_v1')
      call feed_real(hash, cfg%periodic2%max_nonzero_mode_potential_step)
    end if
    ! Preserve the former source-model slot at the only supported model so
    ! existing triangle_p0 checkpoints keep the same ordered fingerprint.
    call feed_string(hash, 'triangle_p0')
    call feed_string(hash, cfg%panel%kernel_id)
    call feed_string(hash, cfg%panel%surface_side_policy)
    call feed_real(hash, cfg%sim%dt)
    call feed_integer(hash, cfg%sim%max_step)
    call feed_real(hash, cfg%sim%batch_duration)
    call feed_logical(hash, cfg%sim%has_batch_duration)
    call feed_real(hash, cfg%sim%batch_duration_step)
    call feed_logical(hash, cfg%sim%has_batch_duration_step)
    ! Preserve the retired softening slot at the triangle_p0 contract value.
    call feed_real(hash, 0.0_dp)
    call feed_real(hash, cfg%sim%field_length_scale)
    call feed_string(hash, cfg%sim%field_bc_mode)
    call feed_integer(hash, cfg%sim%field_periodic_image_layers)
    call feed_real(hash, cfg%sim%field_periodic_ewald_alpha)
    call feed_integer(hash, cfg%sim%field_periodic_ewald_layers)
    call feed_real(hash, cfg%sim%field_periodic_generation_tolerance)
    call feed_real(hash, cfg%sim%tree_theta)
    call feed_integer(hash, cfg%sim%tree_leaf_max)
    call feed_integer(hash, cfg%sim%tree_min_nelem)
    call feed_real_vector(hash, cfg%sim%e0)
    call feed_real_vector(hash, cfg%sim%b0)
    call feed_string(hash, cfg%sim%reservoir_potential_model)
    call feed_real(hash, cfg%sim%phi_infty)
    call feed_string(hash, cfg%sim%open_boundary_model)
    if (trim(cfg%sim%multiple_box_events_policy) /= 'abort') then
      call feed_string(hash, cfg%sim%multiple_box_events_policy)
      call feed_integer(hash, cfg%sim%multiple_box_events_soft_discard_count_limit)
      call feed_real(hash, cfg%sim%multiple_box_events_soft_discard_abs_charge_limit)
    end if
    call feed_integer(hash, cfg%sim%injection_face_phi_grid_n)
    call feed_integer(hash, cfg%sim%raycast_max_bounce)
    call feed_logical(hash, cfg%sim%use_box)
    call feed_real_vector(hash, cfg%sim%box_min)
    call feed_real_vector(hash, cfg%sim%box_max)
    call feed_integer_vector(hash, cfg%sim%bc_low)
    call feed_integer_vector(hash, cfg%sim%bc_high)
    if (any(cfg%particle_boundary_low /= particle_bc_inherit) .or. &
        any(cfg%particle_boundary_high /= particle_bc_inherit)) then
      call feed_string(hash, 'particle_boundary_v1')
      call feed_integer_vector(hash, cfg%particle_boundary_low)
      call feed_integer_vector(hash, cfg%particle_boundary_high)
    end if
    if (uses_redistributed_reflect(cfg)) then
      call feed_string(hash, 'redistributed_reflect_rng_v1')
      call feed_integer(hash, cfg%sim%rng_seed)
    end if
    if (trim(cfg%surface_current%model) /= 'none') then
      call feed_string(hash, 'surface_current_model_v2')
      call feed_string(hash, cfg%surface_current%model)
      call feed_string(hash, cfg%surface_current%zhao_branch)
      call feed_string(hash, cfg%surface_current%electron_species)
      call feed_string(hash, cfg%surface_current%ion_species)
      call feed_string(hash, cfg%surface_current%photoelectron_species)
      call feed_real(hash, cfg%surface_current%solar_elevation_deg)
      call feed_real(hash, cfg%surface_current%photoelectron_ref_density_m3)
      call feed_real(hash, cfg%surface_current%photoelectron_source_scale)
      call feed_real(hash, cfg%surface_current%reference_area_m2)
      call feed_logical(hash, cfg%surface_current%has_reference_area_m2)
    end if
    fingerprint = finish_hash(hash)
  end function model_fingerprint

  !> 面内再配置反射のcounter RNGがmodel fingerprintへseedを要求するかを返す。
  pure logical function uses_redistributed_reflect(cfg) result(uses)
    type(app_config), intent(in) :: cfg
    integer(i32) :: species

    uses = any(cfg%sim%bc_low == bc_redistributed_reflect) .or. &
           any(cfg%sim%bc_high == bc_redistributed_reflect) .or. &
           any(cfg%particle_boundary_low == bc_redistributed_reflect) .or. &
           any(cfg%particle_boundary_high == bc_redistributed_reflect)
    if (uses) return
    do species = 1_i32, cfg%n_particle_species
      uses = any(cfg%particle_species(species)%boundary_low == bc_redistributed_reflect) .or. &
             any(cfg%particle_species(species)%boundary_high == bc_redistributed_reflect)
      if (uses) return
    end do
  end function uses_redistributed_reflect

  function mesh_fingerprint(mesh) result(fingerprint)
    type(mesh_type), intent(in) :: mesh
    character(len=16) :: fingerprint
    type(hash_state) :: hash
    integer(i32) :: elem

    call feed_integer(hash, mesh%nelem)
    do elem = 1, mesh%nelem
      call feed_real_vector(hash, mesh%v0(:, elem))
      call feed_real_vector(hash, mesh%v1(:, elem))
      call feed_real_vector(hash, mesh%v2(:, elem))
      if (allocated(mesh%elem_mesh_id)) call feed_integer(hash, mesh%elem_mesh_id(elem))
      if (allocated(mesh%elem_surface_model)) call feed_integer(hash, mesh%elem_surface_model(elem))
      if (allocated(mesh%elem_epsilon_r)) call feed_real(hash, mesh%elem_epsilon_r(elem))
      if (allocated(mesh%elem_vacuum_sign)) call feed_integer(hash, mesh%elem_vacuum_sign(elem))
    end do
    fingerprint = finish_hash(hash)
  end function mesh_fingerprint

  function species_fingerprint(cfg) result(fingerprint)
    type(app_config), intent(in) :: cfg
    character(len=16) :: fingerprint
    type(hash_state) :: hash
    integer(i32) :: species

    call feed_integer(hash, cfg%n_particle_species)
    do species = 1, cfg%n_particle_species
      call feed_species(hash, cfg%particle_species(species))
    end do
    fingerprint = finish_hash(hash)
  end function species_fingerprint

  subroutine feed_species(hash, spec)
    type(hash_state), intent(inout) :: hash
    type(particle_species_spec), intent(in) :: spec

    call feed_string(hash, spec%species_key)
    call feed_logical(hash, spec%enabled)
    call feed_integer(hash, spec%npcls_per_step)
    call feed_logical(hash, spec%has_npcls_per_step)
    call feed_string(hash, spec%source_mode)
    call feed_real(hash, spec%number_density_cm3)
    call feed_logical(hash, spec%has_number_density_cm3)
    call feed_real(hash, spec%number_density_m3)
    call feed_logical(hash, spec%has_number_density_m3)
    call feed_real(hash, spec%q_particle)
    call feed_real(hash, spec%m_particle)
    call feed_real(hash, spec%w_particle)
    call feed_logical(hash, spec%has_w_particle)
    call feed_integer(hash, spec%target_macro_particles_per_batch)
    call feed_logical(hash, spec%has_target_macro_particles_per_batch)
    call feed_real_vector(hash, spec%pos_low)
    call feed_real_vector(hash, spec%pos_high)
    call feed_string(hash, spec%velocity_distribution)
    call feed_string(hash, spec%velocity_grid_path)
    call feed_string(hash, spec%velocity_grid_pdf_kind)
    call feed_string(hash, spec%velocity_grid_sampling)
    call feed_real(hash, spec%particle_flux_m2_s)
    call feed_logical(hash, spec%has_particle_flux_m2_s)
    call feed_real(hash, spec%current_density_a_m2)
    call feed_logical(hash, spec%has_current_density_a_m2)
    call feed_real_vector(hash, spec%drift_velocity)
    call feed_real(hash, spec%temperature_k)
    call feed_logical(hash, spec%has_temperature_k)
    call feed_real(hash, spec%temperature_ev)
    call feed_logical(hash, spec%has_temperature_ev)
    call feed_real(hash, spec%emit_current_density_a_m2)
    call feed_integer(hash, spec%rays_per_batch)
    call feed_logical(hash, spec%deposit_opposite_charge_on_emit)
    call feed_logical(hash, spec%has_deposit_opposite_charge_on_emit)
    call feed_real(hash, spec%normal_drift_speed)
    call feed_real_vector(hash, spec%ray_direction)
    call feed_logical(hash, spec%has_ray_direction)
    call feed_string(hash, spec%inject_face)
    if (trim(spec%source_mode) == 'plane_source') then
      call feed_string(hash, 'plane_source_v1')
      call feed_real_vector(hash, spec%source_normal)
      call feed_logical(hash, spec%has_source_normal)
    end if
    if (any(spec%boundary_inflow_low /= 0_i32) .or. any(spec%boundary_inflow_high /= 0_i32)) then
      call feed_string(hash, 'species_boundary_inflow_v1')
      call feed_integer_vector(hash, spec%boundary_inflow_low)
      call feed_integer_vector(hash, spec%boundary_inflow_high)
    end if
    if (any(spec%boundary_low /= particle_bc_inherit) .or. &
        any(spec%boundary_high /= particle_bc_inherit)) then
      call feed_string(hash, 'species_particle_boundary_v1')
      call feed_integer_vector(hash, spec%boundary_low)
      call feed_integer_vector(hash, spec%boundary_high)
    end if
    if (trim(spec%surface_charge_closure) /= 'explicit') then
      call feed_string(hash, 'surface_charge_closure_v1')
      call feed_string(hash, spec%surface_charge_closure)
    end if
    if (spec%has_target_absorbed_current_a .or. spec%has_target_emission_current_a) then
      call feed_string(hash, 'fixed_surface_current_v1')
      call feed_real(hash, spec%target_absorbed_current_a)
      call feed_logical(hash, spec%has_target_absorbed_current_a)
      call feed_real(hash, spec%target_emission_current_a)
      call feed_logical(hash, spec%has_target_emission_current_a)
    end if
  end subroutine feed_species

  subroutine feed_string(hash, value)
    type(hash_state), intent(inout) :: hash
    character(len=*), intent(in) :: value
    integer :: index

    call feed_byte(hash, len_trim(value))
    do index = 1, len_trim(value)
      call feed_byte(hash, iachar(value(index:index)))
    end do
  end subroutine feed_string

  subroutine feed_integer(hash, value)
    type(hash_state), intent(inout) :: hash
    integer(i32), intent(in) :: value
    character(len=32) :: encoded

    write (encoded, '(i0)') value
    call feed_string(hash, trim(encoded))
  end subroutine feed_integer

  subroutine feed_integer_vector(hash, values)
    type(hash_state), intent(inout) :: hash
    integer(i32), intent(in) :: values(:)
    integer :: index

    call feed_integer(hash, int(size(values), i32))
    do index = 1, size(values)
      call feed_integer(hash, values(index))
    end do
  end subroutine feed_integer_vector

  subroutine feed_real(hash, value)
    type(hash_state), intent(inout) :: hash
    real(dp), intent(in) :: value
    character(len=32) :: encoded

    write (encoded, '(es24.16e3)') value
    call feed_string(hash, trim(adjustl(encoded)))
  end subroutine feed_real

  subroutine feed_real_vector(hash, values)
    type(hash_state), intent(inout) :: hash
    real(dp), intent(in) :: values(:)
    integer :: index

    call feed_integer(hash, int(size(values), i32))
    do index = 1, size(values)
      call feed_real(hash, values(index))
    end do
  end subroutine feed_real_vector

  subroutine feed_logical(hash, value)
    type(hash_state), intent(inout) :: hash
    logical, intent(in) :: value

    if (value) then
      call feed_string(hash, 'true')
    else
      call feed_string(hash, 'false')
    end if
  end subroutine feed_logical

  subroutine feed_byte(hash, value)
    type(hash_state), intent(inout) :: hash
    integer, intent(in) :: value

    hash%a = modulo(hash%a*hash_multiplier_a + int(value, i64) + 1_i64, hash_modulus)
    hash%b = modulo(hash%b*hash_multiplier_b + int(value, i64) + 1_i64, hash_modulus)
  end subroutine feed_byte

  function finish_hash(hash) result(fingerprint)
    type(hash_state), intent(in) :: hash
    character(len=16) :: fingerprint

    write (fingerprint, '(z8.8,z8.8)') hash%a, hash%b
  end function finish_hash

end module bem_model_fingerprint
