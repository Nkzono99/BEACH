!> lint 済み設定の基本的な計算条件を確認し、時間と外部場を確定する。
submodule(bem_app_config_parser) bem_app_config_parser_preflight_sim
  implicit none
contains

  module procedure validate_simulation_config
  integer :: axis
  character(len=6), parameter :: low_faces(3) = ['x_low ', 'y_low ', 'z_low ']
  character(len=6), parameter :: high_faces(3) = ['x_high', 'y_high', 'z_high']

  if (cfg%sim%dt <= 0.0_dp) error stop 'particles.tracking.dt_s must be > 0.'
  if (cfg%sim%q_floor <= 0.0_dp) error stop 'output.diagnostics.charge_floor_c must be > 0.'
  if (cfg%sim%use_box) then
    if (.not. all(ieee_is_finite(cfg%sim%box_min)) .or. .not. all(ieee_is_finite(cfg%sim%box_max)) .or. &
        any(cfg%sim%box_max <= cfg%sim%box_min)) then
      error stop '[domain] must resolve to finite bounds with positive dimensions.'
    end if
  else if (any(cfg%particle_boundary_low /= particle_bc_inherit) .or. &
           any(cfg%particle_boundary_high /= particle_bc_inherit)) then
    error stop '[particles.boundary] requires a finite [domain].'
  end if
  ! 周期軸の面は domain.periodic_axes が決める。粒子境界の上書きは面の意味を変えるので拒否する。
  do axis = 1, 3
    if (cfg%sim%bc_low(axis) /= bc_periodic) cycle
    if (cfg%particle_boundary_low(axis) /= particle_bc_inherit) then
      error stop 'particles.boundary.'//trim(low_faces(axis))//' cannot override a periodic domain face.'
    end if
    if (cfg%particle_boundary_high(axis) /= particle_bc_inherit) then
      error stop 'particles.boundary.'//trim(high_faces(axis))//' cannot override a periodic domain face.'
    end if
  end do
  if (cfg%n_particle_species <= 0_i32) error stop 'At least one [[particles.species]] entry is required.'
  if (cfg%sim%injection_face_phi_grid_n < 1_i32) error stop 'particles.reservoir.face_potential_grid_n must be >= 1.'
  call resolve_external_e_field(cfg)
  call resolve_batch_duration(cfg)
  end procedure validate_simulation_config

end submodule bem_app_config_parser_preflight_sim
