!> lint 済み設定の基本的な計算条件を確認し、時間と外部場を確定する。
submodule(bem_app_config_parser) bem_app_config_parser_preflight_sim
  implicit none
contains

  module procedure validate_simulation_config
  if (cfg%sim%dt <= 0.0_dp) error stop 'sim.dt must be > 0.'
  if (cfg%sim%q_floor <= 0.0_dp) error stop 'sim.q_floor must be > 0.'
  if (cfg%sim%use_box) then
    if (.not. all(ieee_is_finite(cfg%sim%box_min)) .or. .not. all(ieee_is_finite(cfg%sim%box_max)) .or. &
        any(cfg%sim%box_max <= cfg%sim%box_min)) then
      error stop '[domain] must resolve to finite bounds with positive dimensions.'
    end if
  end if
  if (cfg%n_particle_species <= 0_i32) error stop 'At least one [[particles.species]] entry is required.'
  if (cfg%sim%injection_face_phi_grid_n < 1_i32) error stop 'reservoir.face_potential_grid_n must be >= 1.'
  call resolve_external_e_field(cfg)
  call resolve_batch_duration(cfg)
  end procedure validate_simulation_config

end submodule bem_app_config_parser_preflight_sim
