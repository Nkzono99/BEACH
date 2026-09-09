!> lint 済みの粒子種識別子と有効な粒子源の派生量を確定する。
submodule(bem_app_config_parser) bem_app_config_parser_preflight_particles
  use bem_app_config_types, only: particle_inflow_none
  implicit none
contains

  module procedure validate_particle_species_config
  integer :: i
  character(len=64) :: generated_species_key
  logical :: has_boundary_inflow

  do i = 1, cfg%n_particle_species
    if (len_trim(cfg%particle_species(i)%species_key) == 0) then
      write (generated_species_key, '(a,i0)') 'species_', i
      cfg%particle_species(i)%species_key = trim(generated_species_key)
    end if
    if (.not. cfg%particle_species(i)%enabled) cycle

    cfg%particle_species(i)%source_mode = lower_ascii(trim(cfg%particle_species(i)%source_mode))
    cfg%particle_species(i)%velocity_distribution = lower_ascii(trim(cfg%particle_species(i)%velocity_distribution))
    cfg%particle_species(i)%velocity_grid_pdf_kind = lower_ascii(trim(cfg%particle_species(i)%velocity_grid_pdf_kind))
    cfg%particle_species(i)%velocity_grid_sampling = lower_ascii(trim(cfg%particle_species(i)%velocity_grid_sampling))
    cfg%particle_species(i)%surface_charge_closure = &
      lower_ascii(trim(cfg%particle_species(i)%surface_charge_closure))
    has_boundary_inflow = any(cfg%particle_species(i)%boundary_inflow_low /= particle_inflow_none) .or. &
                          any(cfg%particle_species(i)%boundary_inflow_high /= particle_inflow_none)

    ! 位置は box/face 相対指定の展開結果なので、入力値の有限性だけでは保証できない。
    if (.not. all(ieee_is_finite(cfg%particle_species(i)%pos_low)) .or. &
        .not. all(ieee_is_finite(cfg%particle_species(i)%pos_high))) then
      error stop 'particles.species.pos_low/pos_high must contain finite values.'
    end if
    select case (trim(cfg%particle_species(i)%source_mode))
    case ('volume_seed')
      if (has_boundary_inflow) call validate_boundary_inflow_species(cfg, i)
    case ('reservoir_face')
      call validate_reservoir_species(cfg, i)
    case ('plane_source')
      call validate_plane_source_species(cfg, i)
    case ('photo_raycast')
      call validate_photo_raycast_species(cfg, i)
    end select
  end do
  end procedure validate_particle_species_config

end submodule bem_app_config_parser_preflight_particles
