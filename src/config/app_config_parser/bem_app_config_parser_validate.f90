!> lint 済みの設定から電場・流束・粒子重みを導出し、計算結果を検査する。
submodule(bem_app_config_parser) bem_app_config_parser_validate
  use bem_config_helpers, only: resolve_inward_normal, species_number_density_m3, species_temperature_k
  implicit none
contains

  !> 時間刻み数から batch_duration を導出し、積の overflow/underflow を検査する。
  module procedure resolve_batch_duration
  real(dp) :: batch_duration

  if (cfg%sim%has_batch_duration .and. cfg%sim%has_batch_duration_step) then
    error stop 'run.batch.duration_s and run.batch.duration_steps cannot be specified together.'
  end if
  if (cfg%sim%has_batch_duration_step) then
    batch_duration = cfg%sim%dt*cfg%sim%batch_duration_step
    if (.not. ieee_is_finite(batch_duration) .or. batch_duration <= 0.0d0) then
      error stop 'run.batch.duration_steps produced invalid run.batch.duration_s.'
    end if
    cfg%sim%batch_duration = batch_duration
    cfg%sim%has_batch_duration = .true.
  end if
  end procedure resolve_batch_duration

  !> 一様外部電場の大きさ・角度を `sim%e0` へ変換する。
  module procedure resolve_external_e_field
  real(dp), parameter :: deg2rad = acos(-1.0d0)/180.0d0
  real(dp) :: phi_xy, phi_z

  if (cfg%sim%has_e0_vector .and. (cfg%sim%has_e0_abs .or. cfg%sim%has_e0_phi_xy_deg .or. cfg%sim%has_e0_phi_z_deg)) then
    error stop 'fields.external.electric_v_m cannot be combined with '// &
      'fields.external.electric_magnitude_v_m/fields.external.electric_azimuth_deg/fields.external.electric_elevation_deg.'
  end if
  if ((cfg%sim%has_e0_phi_xy_deg .or. cfg%sim%has_e0_phi_z_deg) .and. .not. cfg%sim%has_e0_abs) then
    error stop 'fields.external.electric_azimuth_deg/fields.external.electric_elevation_deg require '// &
      'fields.external.electric_magnitude_v_m.'
  end if
  if (cfg%sim%has_e0_vector) return
  if (.not. cfg%sim%has_e0_abs) then
    cfg%sim%e0 = [0.0d0, 0.0d0, 0.0d0]
    return
  end if

  phi_xy = cfg%sim%e0_phi_xy_deg*deg2rad
  phi_z = cfg%sim%e0_phi_z_deg*deg2rad
  cfg%sim%e0(1) = cfg%sim%e0_abs*cos(phi_z)*cos(phi_xy)
  cfg%sim%e0(2) = cfg%sim%e0_abs*cos(phi_z)*sin(phi_xy)
  cfg%sim%e0(3) = cfg%sim%e0_abs*sin(phi_z)
  end procedure resolve_external_e_field

  !> 有効な box 面からの総流入率を求め、粒子重みを導出する。
  module procedure validate_boundary_inflow_species
  type(particle_species_spec) :: spec
  real(dp) :: physical_rate, area, inward_normal(3), gamma_in
  integer :: face
  character(len=6) :: face_name

  spec = cfg%particle_species(species_idx)
  call validate_flux_driven_parameters(cfg, species_idx, spec)
  physical_rate = 0.0_dp
  do face = 1, 6
    if (.not. boundary_inflow_face_enabled(spec, face)) cycle
    call boundary_face_name(face, face_name)
    call resolve_inward_normal(face_name, inward_normal)
    area = compute_face_area_from_bounds(face_name, cfg%sim%box_min, cfg%sim%box_max)
    gamma_in = compute_inflow_flux_from_drifting_maxwellian( &
               species_number_density_m3(spec), species_temperature_k(spec), spec%m_particle, &
               spec%drift_velocity, inward_normal &
               )
    physical_rate = physical_rate + gamma_in*area
  end do
  call resolve_flux_driven_weight(cfg, spec, physical_rate, species_idx)
  cfg%particle_species(species_idx) = spec
  end procedure validate_boundary_inflow_species

  !> photo_raycast の面積と正規化した入射方向を確定する。
  module procedure validate_photo_raycast_species
  real(dp) :: area, direction_norm, inward_dot
  real(dp) :: inward_normal(3)
  type(particle_species_spec) :: spec

  spec = cfg%particle_species(species_idx)

  area = compute_face_area_from_bounds(spec%inject_face, spec%pos_low, spec%pos_high)
  if (.not. ieee_is_finite(area) .or. area <= 0.0d0) then
    error stop 'photo_raycast opening area must be positive.'
  end if

  call resolve_inward_normal(spec%inject_face, inward_normal)
  if (spec%has_ray_direction) then
    direction_norm = sqrt(sum(spec%ray_direction*spec%ray_direction))
    if (.not. ieee_is_finite(direction_norm) .or. direction_norm <= 0.0d0) then
      error stop species_config_path(species_idx)//'.source.ray_direction norm must be > 0.'
    end if
    spec%ray_direction = spec%ray_direction/direction_norm
  else
    spec%ray_direction = inward_normal
  end if
  inward_dot = dot_product(spec%ray_direction, inward_normal)
  if (.not. ieee_is_finite(inward_dot) .or. inward_dot <= 0.0d0) then
    error stop species_config_path(species_idx)//'.source.ray_direction must point inward from '// &
      species_config_path(species_idx)//'.source.inject_face.'
  end if

  cfg%particle_species(species_idx) = spec
  end procedure validate_photo_raycast_species

  subroutine validate_flux_driven_parameters(cfg, species_idx, spec)
    type(app_config), intent(in) :: cfg
    integer, intent(in) :: species_idx
    type(particle_species_spec), intent(inout) :: spec

    if (spec%has_target_macro_particles_per_batch .and. spec%target_macro_particles_per_batch /= -1_i32) then
      if (spec%target_macro_particles_per_batch <= 0_i32) then
        error stop species_config_path(species_idx)// &
          '.sampling.target_macro_particles_per_batch must be positive for weight division.'
      end if
    end if
    if (spec%has_target_macro_particles_per_batch .and. spec%target_macro_particles_per_batch == -1_i32) then
      if (species_idx == 1) error stop 'particles.species[1].sampling.target_macro_particles_per_batch cannot be -1.'
      if (.not. cfg%particle_species(1)%enabled .or. .not. cfg%particle_species(1)%has_w_particle .or. &
          .not. ieee_is_finite(cfg%particle_species(1)%w_particle) .or. cfg%particle_species(1)%w_particle <= 0.0_dp) then
        error stop species_config_path(species_idx)// &
          '.sampling.target_macro_particles_per_batch=-1 requires a positive particles.species[1].sampling.weight.'
      end if
    end if

    if (spec%has_number_density_cm3 .and. spec%has_number_density_m3) then
      error stop 'Specify either '//species_config_path(species_idx)//'.distribution.number_density_cm3 or '// &
        species_config_path(species_idx)//'.distribution.number_density_m3, not both.'
    end if
    if (.not. ieee_is_finite(species_number_density_m3(spec)) .or. species_number_density_m3(spec) <= 0.0_dp) then
      error stop species_config_path(species_idx)// &
        '.distribution must resolve to a finite number density > 0.'
    end if
    if (.not. ieee_is_finite(species_temperature_k(spec)) .or. species_temperature_k(spec) < 0.0_dp) then
      error stop species_config_path(species_idx)// &
        '.distribution must resolve to a finite temperature >= 0.'
    end if
  end subroutine validate_flux_driven_parameters

  subroutine resolve_flux_driven_weight(cfg, spec, physical_rate, species_idx)
    type(app_config), intent(in) :: cfg
    type(particle_species_spec), intent(inout) :: spec
    real(dp), intent(in) :: physical_rate
    integer, intent(in) :: species_idx

    if (.not. spec%has_target_macro_particles_per_batch) return
    if (spec%target_macro_particles_per_batch == -1_i32) then
      spec%w_particle = cfg%particle_species(1)%w_particle
    else
      spec%w_particle = physical_rate*cfg%sim%batch_duration/real(spec%target_macro_particles_per_batch, dp)
    end if
    if (.not. ieee_is_finite(spec%w_particle) .or. spec%w_particle <= 0.0_dp) then
      call stop_config_error(species_config_path(species_idx)// &
                             '.sampling.target_macro_particles_per_batch produced invalid '// &
                             species_config_path(species_idx)//'.sampling.weight.')
    end if
    spec%has_w_particle = .true.
  end subroutine resolve_flux_driven_weight

  pure logical function boundary_inflow_face_enabled(spec, face) result(enabled)
    type(particle_species_spec), intent(in) :: spec
    integer, intent(in) :: face

    if (mod(face, 2) == 1) then
      enabled = spec%boundary_inflow_low((face + 1)/2) == particle_inflow_reservoir
    else
      enabled = spec%boundary_inflow_high(face/2) == particle_inflow_reservoir
    end if
  end function boundary_inflow_face_enabled

  pure subroutine boundary_face_name(face, name)
    integer, intent(in) :: face
    character(len=*), intent(out) :: name

    select case (face)
    case (1)
      name = 'x_low'
    case (2)
      name = 'x_high'
    case (3)
      name = 'y_low'
    case (4)
      name = 'y_high'
    case (5)
      name = 'z_low'
    case (6)
      name = 'z_high'
    case default
      error stop 'invalid boundary face index.'
    end select
  end subroutine boundary_face_name

end submodule bem_app_config_parser_validate
