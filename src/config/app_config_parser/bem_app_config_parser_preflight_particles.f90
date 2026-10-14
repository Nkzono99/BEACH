!> lint 済みの粒子種識別子と有効な粒子源の派生量を確定する。
submodule(bem_app_config_parser) bem_app_config_parser_preflight_particles
  use bem_app_config_types, only: particle_inflow_none, particles_per_batch_from_config
  use bem_config_helpers, only: resolve_particle_boundaries, particle_boundary_action_for_face
  implicit none
contains

  module procedure validate_particle_species_config
  integer :: i, j, axis
  integer(i32) :: batch_n, effective_low(3), effective_high(3)
  character(len=64) :: generated_species_key
  logical :: has_boundary_inflow

  do i = 1, cfg%n_particle_species
    if (len_trim(cfg%particle_species(i)%species_key) == 0) then
      write (generated_species_key, '(a,i0)') 'species_', i
      cfg%particle_species(i)%species_key = trim(generated_species_key)
    end if
    do j = 1, i - 1
      if (trim(cfg%particle_species(j)%species_key) == trim(cfg%particle_species(i)%species_key)) then
        call stop_config_error('particles.species.species_key values must be unique: '// &
                               trim(cfg%particle_species(i)%species_key))
      end if
    end do
    if (.not. cfg%particle_species(i)%enabled) cycle

    cfg%particle_species(i)%source_mode = lower_ascii(trim(cfg%particle_species(i)%source_mode))
    cfg%particle_species(i)%surface_charge_closure = &
      lower_ascii(trim(cfg%particle_species(i)%surface_charge_closure))
    has_boundary_inflow = any(cfg%particle_species(i)%boundary_inflow_low /= particle_inflow_none) .or. &
                          any(cfg%particle_species(i)%boundary_inflow_high /= particle_inflow_none)

    if (.not. ieee_is_finite(cfg%particle_species(i)%q_particle) .or. cfg%particle_species(i)%q_particle == 0.0_dp) then
      error stop 'particles.species.q_particle must be finite and non-zero for an enabled species.'
    end if
    if (.not. cfg%sim%use_box .and. &
        (any(cfg%particle_species(i)%boundary_low /= particle_bc_inherit) .or. &
         any(cfg%particle_species(i)%boundary_high /= particle_bc_inherit))) then
      error stop 'particles.species.boundary requires a finite [domain].'
    end if
    do axis = 1, 3
      if (cfg%sim%bc_low(axis) /= bc_periodic) cycle
      if (cfg%particle_species(i)%boundary_low(axis) /= particle_bc_inherit .or. &
          cfg%particle_species(i)%boundary_high(axis) /= particle_bc_inherit) then
        error stop 'particles.species.boundary cannot override a periodic domain face.'
      end if
      if (cfg%particle_species(i)%boundary_inflow_low(axis) /= particle_inflow_none .or. &
          cfg%particle_species(i)%boundary_inflow_high(axis) /= particle_inflow_none) then
        error stop 'particles.species.boundary_inflow cannot inject through a periodic domain face.'
      end if
    end do
    call resolve_particle_boundaries( &
      cfg%sim, cfg%particle_boundary_low, cfg%particle_boundary_high, cfg%particle_species(i), &
      effective_low, effective_high &
      )
    do axis = 1, 3
      if ((cfg%particle_species(i)%boundary_inflow_low(axis) /= particle_inflow_none .and. &
           effective_low(axis) /= bc_open) .or. &
          (cfg%particle_species(i)%boundary_inflow_high(axis) /= particle_inflow_none .and. &
           effective_high(axis) /= bc_open)) then
        error stop 'particles.species.boundary_inflow requires an open particle boundary action on the inflow face.'
      end if
    end do
    if (has_boundary_inflow .and. trim(cfg%particle_species(i)%source_mode) /= 'volume_seed') then
      call stop_config_error('particles.species.boundary_inflow cannot be combined with source_mode="'// &
                             trim(cfg%particle_species(i)%source_mode)//'".')
    end if

    ! 位置は box/face 相対指定の展開結果なので、入力値の有限性だけでは保証できない。
    if (.not. all(ieee_is_finite(cfg%particle_species(i)%pos_low)) .or. &
        .not. all(ieee_is_finite(cfg%particle_species(i)%pos_high))) then
      error stop 'particles.species.pos_low/pos_high must contain finite values.'
    end if
    select case (trim(cfg%particle_species(i)%source_mode))
    case ('volume_seed')
      if (has_boundary_inflow) call validate_boundary_inflow_species(cfg, i)
    case ('photo_raycast')
      call validate_photo_raycast_species(cfg, i)
    end select
    call validate_surface_charge_closure(cfg, i, effective_low, effective_high)
  end do
  batch_n = particles_per_batch_from_config(cfg)
  end procedure validate_particle_species_config

  !> 表面 source 電荷 closure が、その closure を適用できる粒子源と目標電流を持つことを確認する。
  subroutine validate_surface_charge_closure(cfg, species_idx, effective_low, effective_high)
    type(app_config), intent(in) :: cfg
    integer, intent(in) :: species_idx
    integer(i32), intent(in) :: effective_low(3), effective_high(3)
    integer(i32) :: inject_action
    logical :: automatic_current_species

    associate (spec => cfg%particle_species(species_idx))
      select case (trim(spec%surface_charge_closure))
      case ('explicit')
        if (spec%has_target_absorbed_current_a .or. spec%has_target_emission_current_a) then
          error stop 'target surface currents require surface_charge_closure="fixed_current".'
        end if
      case ('fixed_current')
        if (.not. ieee_is_finite(cfg%sim%batch_duration) .or. cfg%sim%batch_duration <= 0.0_dp) then
          error stop 'sim.batch_duration must be > 0 for fixed_current.'
        end if
        automatic_current_species = trim(lower_ascii(cfg%surface_current%model)) == 'zhao_stationary' .and. &
                                    (trim(spec%species_key) == trim(cfg%surface_current%electron_species) .or. &
                                     trim(spec%species_key) == trim(cfg%surface_current%ion_species) .or. &
                                     trim(spec%species_key) == trim(cfg%surface_current%photoelectron_species))
        if (.not. (spec%has_target_absorbed_current_a .or. spec%has_target_emission_current_a .or. &
                   automatic_current_species)) then
          error stop 'fixed_current requires at least one target current.'
        end if
        if (spec%has_target_absorbed_current_a) then
          if (.not. ieee_is_finite(spec%target_absorbed_current_a) .or. &
              spec%target_absorbed_current_a*spec%q_particle < 0.0_dp) then
            error stop 'target_absorbed_current_a must be finite and have the same sign as q_particle.'
          end if
        end if
        if (spec%has_target_emission_current_a) then
          if (trim(spec%source_mode) /= 'photo_raycast' .or. .not. spec%deposit_opposite_charge_on_emit .or. &
              .not. ieee_is_finite(spec%target_emission_current_a) .or. &
              spec%target_emission_current_a*spec%q_particle > 0.0_dp) then
            error stop 'target_emission_current_a requires photo_raycast, an opposite-charge emission deposit, '// &
              'and current sign opposite to q_particle.'
          end if
        end if
      case ('neutral_return')
        if (.not. cfg%sim%use_box) error stop 'neutral_return requires a finite [domain].'
        inject_action = particle_boundary_action_for_face(effective_low, effective_high, spec%inject_face)
        if (trim(spec%source_mode) /= 'photo_raycast' .or. spec%q_particle >= 0.0_dp .or. &
            .not. spec%deposit_opposite_charge_on_emit .or. &
            (inject_action /= bc_reflect .and. inject_action /= bc_redistributed_reflect)) then
          error stop 'neutral_return requires a negative photo_raycast species, deposit_opposite_charge_on_emit=true, '// &
            'and a reflecting action on inject_face.'
        end if
      case default
        error stop 'particles.species.surface_charge_closure must be "explicit", "fixed_current", or "neutral_return".'
      end select
    end associate
  end subroutine validate_surface_charge_closure

end submodule bem_app_config_parser_preflight_particles
