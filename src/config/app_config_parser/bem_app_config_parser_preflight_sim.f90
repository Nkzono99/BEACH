!> 実行時間・領域・出力と境界条件の preflight。
submodule(bem_app_config_parser) bem_app_config_parser_preflight_sim
  implicit none
contains

  module procedure validate_simulation_config
  integer :: axis
  if (cfg%sim%batch_count <= 0_i32) error stop 'sim.batch_count must be > 0.'
  if (.not. ieee_is_finite(cfg%sim%dt) .or. cfg%sim%dt <= 0.0d0) then
    error stop 'sim.dt must be finite and > 0.'
  end if
  if (cfg%sim%max_step <= 0_i32) error stop 'sim.max_step must be > 0.'
  if (.not. ieee_is_finite(cfg%sim%tol_rel) .or. cfg%sim%tol_rel < 0.0d0) then
    error stop 'sim.tol_rel must be finite and >= 0.'
  end if
  if (.not. ieee_is_finite(cfg%sim%q_floor) .or. cfg%sim%q_floor <= 0.0d0) then
    error stop 'sim.q_floor must be finite and > 0.'
  end if
  if (.not. all(ieee_is_finite(cfg%sim%b0))) error stop 'sim.b0 must contain finite values.'
  if (cfg%sim%use_box) then
    if (.not. all(ieee_is_finite(cfg%sim%box_min)) .or. .not. all(ieee_is_finite(cfg%sim%box_max))) then
      error stop 'domain.box_min/box_max must contain finite values.'
    end if
    if (any(cfg%sim%box_max <= cfg%sim%box_min)) then
      error stop 'domain.box_max must be greater than domain.box_min on all axes in [domain].'
    end if
  end if
  do axis = 1, 3
    call validate_particle_boundary_override( &
      cfg%particle_boundary_low(axis), cfg%sim%bc_low(axis), 'particle_boundary low face' &
      )
    call validate_particle_boundary_override( &
      cfg%particle_boundary_high(axis), cfg%sim%bc_high(axis), 'particle_boundary high face' &
      )
  end do
  if (.not. cfg%sim%use_box .and. &
      (any(cfg%particle_boundary_low /= particle_bc_inherit) .or. &
       any(cfg%particle_boundary_high /= particle_bc_inherit))) then
    error stop '[particle_boundary] requires a finite [domain].'
  end if
  if (cfg%n_particle_species <= 0_i32) error stop 'At least one [[particles.species]] entry is required.'
  if (cfg%resume_output .and. .not. cfg%write_output) then
    error stop 'output.resume requires output.write_files = true.'
  end if
  if (len_trim(cfg%output_restart_from) > 0 .and. .not. cfg%resume_output) then
    error stop 'output.restart_from requires output.resume = true.'
  end if
  if (cfg%checkpoint_stride < 0_i32) then
    error stop 'output.checkpoint_stride must be >= 0.'
  end if
  if (cfg%checkpoint_stride > 0_i32 .and. .not. cfg%write_output) then
    error stop 'output.checkpoint_stride > 0 requires output.write_files = true.'
  end if
  call resolve_external_e_field(cfg)
  cfg%sim%reservoir_potential_model = lower_ascii(trim(cfg%sim%reservoir_potential_model))
  select case (trim(cfg%sim%reservoir_potential_model))
  case ('none', 'infinity_barrier')
    continue
  case default
    error stop 'reservoir.inflow_model must be "none" or "infinity_barrier".'
  end select
  cfg%sim%open_boundary_model = lower_ascii(trim(cfg%sim%open_boundary_model))
  select case (trim(cfg%sim%open_boundary_model))
  case ('escape', 'potential_barrier')
    continue
  case default
    error stop 'particle_boundary.ordinary_open_model must be "escape" or "potential_barrier".'
  end select
  cfg%sim%multiple_box_events_policy = lower_ascii(trim(cfg%sim%multiple_box_events_policy))
  select case (trim(cfg%sim%multiple_box_events_policy))
  case ('abort', 'soft_discard')
    continue
  case default
    error stop 'sim.multiple_box_events_policy must be "abort" or "soft_discard".'
  end select
  cfg%sim%multiple_box_events_retry_backend = lower_ascii(trim(cfg%sim%multiple_box_events_retry_backend))
  select case (trim(cfg%sim%multiple_box_events_retry_backend))
  case ('none', 'upper_panel_fourier')
    continue
  case default
    error stop 'sim.multiple_box_events_retry_backend must be "none" or "upper_panel_fourier".'
  end select
  if (trim(cfg%sim%multiple_box_events_policy) == 'soft_discard') then
    if (cfg%sim%multiple_box_events_soft_discard_count_grace < 0_i32) then
      error stop 'sim.multiple_box_events_soft_discard_count_grace must be >= 0 for soft_discard.'
    end if
    if (.not. ieee_is_finite(cfg%sim%multiple_box_events_soft_discard_fraction_limit) .or. &
        cfg%sim%multiple_box_events_soft_discard_fraction_limit <= 0.0_dp .or. &
        cfg%sim%multiple_box_events_soft_discard_fraction_limit > 1.0_dp) then
      error stop 'sim.multiple_box_events_soft_discard_fraction_limit must be finite and in (0, 1] for soft_discard.'
    end if
    if (.not. ieee_is_finite(cfg%sim%multiple_box_events_soft_discard_abs_charge_limit) .or. &
        cfg%sim%multiple_box_events_soft_discard_abs_charge_limit <= 0.0_dp) then
      error stop 'sim.multiple_box_events_soft_discard_abs_charge_limit must be finite and > 0 for soft_discard.'
    end if
  end if
  if (cfg%sim%injection_face_phi_grid_n < 1_i32) then
    error stop 'reservoir.face_potential_grid_n must be >= 1.'
  end if
  if (.not. ieee_is_finite(cfg%sim%phi_infty)) then
    error stop 'reservoir.phi_infty must be finite.'
  end if
  call resolve_batch_duration(cfg)
  end procedure validate_simulation_config

end submodule bem_app_config_parser_preflight_sim
