!> lint 済みの粒子種識別子と有効な粒子源の派生量を確定する。
submodule(bem_app_config_parser) bem_app_config_parser_preflight_particles
  use bem_app_config_types, only: particle_inflow_none, particles_per_batch_from_config
  use bem_config_helpers, only: resolve_particle_boundaries, particle_boundary_action_for_face, species_temperature_k
  implicit none
contains

  !> 粒子種ごとに key を確定し、境界・粒子源・電荷 closure を検査する。二次電子の target は全種の key が
  !! 確定してから解決する（2 パス）。最後に batch 粒子数の上限を検査する。
  module procedure validate_particle_species_config
  integer :: i
  integer(i32) :: batch_n, effective_low(3), effective_high(3)
  logical :: has_boundary_inflow

  if (cfg%secondary_emission_max_generations < 0_i32) then
    call stop_config_error("particles.secondary_emission.max_generations must be >= 0.")
  end if
  do i = 1, cfg%n_particle_species
    call resolve_species_key(cfg, i)
    if (.not. cfg%particle_species(i)%enabled) cycle

    cfg%particle_species(i)%source_mode = lower_ascii(trim(cfg%particle_species(i)%source_mode))
    cfg%particle_species(i)%surface_charge_closure = &
      lower_ascii(trim(cfg%particle_species(i)%surface_charge_closure))
    has_boundary_inflow = any(cfg%particle_species(i)%boundary_inflow_low /= particle_inflow_none) .or. &
                          any(cfg%particle_species(i)%boundary_inflow_high /= particle_inflow_none)
    call validate_species_boundaries(cfg, i, effective_low, effective_high)
    call validate_species_source(cfg, i, has_boundary_inflow)
    call validate_surface_charge_closure(cfg, i, effective_low, effective_high)
  end do
  ! Resolve targets after every species has received its unique key.
  do i = 1, cfg%n_particle_species
    if (.not. cfg%particle_species(i)%enabled .or. .not. cfg%particle_species(i)%secondary_emission_enabled) cycle
    call resolve_secondary_emission_target(cfg, i)
  end do
  batch_n = particles_per_batch_from_config(cfg)
  end procedure validate_particle_species_config

  !> 省略された species key を `species_<i>` で補い、それより前の種と重複しないことを確かめる。
  subroutine resolve_species_key(cfg, species_idx)
    type(app_config), intent(inout) :: cfg
    integer, intent(in) :: species_idx
    character(len=64) :: generated_species_key
    integer :: j

    if (len_trim(cfg%particle_species(species_idx)%species_key) == 0) then
      write (generated_species_key, '(a,i0)') 'species_', species_idx
      cfg%particle_species(species_idx)%species_key = trim(generated_species_key)
    end if
    do j = 1, species_idx - 1
      if (trim(cfg%particle_species(j)%species_key) == trim(cfg%particle_species(species_idx)%species_key)) then
        call stop_config_error(species_config_path(species_idx)//'.species_key values must be unique: '// &
                               trim(cfg%particle_species(species_idx)%species_key))
      end if
    end do
  end subroutine resolve_species_key

  !> 電荷と粒子境界の組合せ（有限 box、周期面、流入面の open 動作）を検査し、実効境界動作を返す。
  subroutine validate_species_boundaries(cfg, species_idx, effective_low, effective_high)
    type(app_config), intent(in) :: cfg
    integer, intent(in) :: species_idx
    integer(i32), intent(out) :: effective_low(3), effective_high(3)
    integer :: axis

    associate (spec => cfg%particle_species(species_idx))
      if (.not. ieee_is_finite(spec%q_particle) .or. spec%q_particle == 0.0_dp) then
        error stop species_config_path(species_idx)//'.charge_c must be finite and non-zero for an enabled species.'
      end if
      if (.not. cfg%sim%use_box .and. &
          (any(spec%boundary_low /= particle_bc_inherit) .or. any(spec%boundary_high /= particle_bc_inherit))) then
        error stop species_config_path(species_idx)//'.boundary requires a finite [domain].'
      end if
      do axis = 1, 3
        if (cfg%sim%bc_low(axis) /= bc_periodic) cycle
        if (spec%boundary_low(axis) /= particle_bc_inherit .or. spec%boundary_high(axis) /= particle_bc_inherit) then
          error stop species_config_path(species_idx)//'.boundary cannot override a periodic domain face.'
        end if
        if (spec%boundary_inflow_low(axis) /= particle_inflow_none .or. &
            spec%boundary_inflow_high(axis) /= particle_inflow_none) then
          error stop species_config_path(species_idx)//'.inflow cannot inject through a periodic domain face.'
        end if
      end do
      call resolve_particle_boundaries( &
        cfg%sim, cfg%particle_boundary_low, cfg%particle_boundary_high, spec, effective_low, effective_high &
        )
      do axis = 1, 3
        if ((spec%boundary_inflow_low(axis) /= particle_inflow_none .and. effective_low(axis) /= bc_open) .or. &
            (spec%boundary_inflow_high(axis) /= particle_inflow_none .and. effective_high(axis) /= bc_open)) then
          error stop species_config_path(species_idx)//'.inflow requires an open particle boundary action on the inflow face.'
        end if
      end do
    end associate
  end subroutine validate_species_boundaries

  !> source mode ごとの条件を検査し、流入重みと photo ray 方向を確定する。
  subroutine validate_species_source(cfg, species_idx, has_boundary_inflow)
    type(app_config), intent(inout) :: cfg
    integer, intent(in) :: species_idx
    logical, intent(in) :: has_boundary_inflow

    if (has_boundary_inflow .and. trim(cfg%particle_species(species_idx)%source_mode) /= 'volume_seed') then
      call stop_config_error(species_config_path(species_idx)//'.inflow cannot be combined with '// &
                             species_config_path(species_idx)//'.source.mode="'// &
                             trim(cfg%particle_species(species_idx)%source_mode)//'".')
    end if
    ! 位置は box/face 相対指定の展開結果なので、入力値の有限性だけでは保証できない。
    if (.not. all(ieee_is_finite(cfg%particle_species(species_idx)%pos_low)) .or. &
        .not. all(ieee_is_finite(cfg%particle_species(species_idx)%pos_high))) then
      error stop species_config_path(species_idx)//'.source.pos_low/'//species_config_path(species_idx)// &
        '.source.pos_high must contain finite values.'
    end if
    select case (trim(cfg%particle_species(species_idx)%source_mode))
    case ('volume_seed')
      call validate_volume_seed_species(cfg%particle_species(species_idx), species_idx, has_boundary_inflow)
      if (has_boundary_inflow) call resolve_boundary_inflow_weight(cfg, species_idx)
    case ('photo_raycast')
      call resolve_photo_raycast_direction(cfg, species_idx)
    case ('secondary_emission')
      associate (spec => cfg%particle_species(species_idx))
        if (spec%has_secondary_forbidden_source_option) then
          call stop_config_error('secondary_emission species must omit sampling, inflow, and emission source options.')
        end if
        if (spec%has_surface_charge_closure .or. trim(spec%surface_charge_closure) /= 'explicit') then
          call stop_config_error('secondary_emission species must omit '//species_config_path(species_idx)// &
                                 '.charging.closure and remain fully tracked.')
        end if
        if (spec%q_particle >= 0.0_dp .or. spec%m_particle <= 0.0_dp .or. species_temperature_k(spec) < 0.0_dp) then
          call stop_config_error('secondary_emission requires negative charge, positive mass, and nonnegative temperature.')
        end if
      end associate
    case default
      call stop_config_error('Unknown '//species_config_path(species_idx)//'.source.mode.')
    end select
  end subroutine validate_species_source

  !> 二次電子の収率モデルを検査し、target species を有効な secondary_emission 種の index へ解決する。
  subroutine resolve_secondary_emission_target(cfg, species_idx)
    type(app_config), intent(inout) :: cfg
    integer, intent(in) :: species_idx
    integer :: j

    associate (spec => cfg%particle_species(species_idx))
      select case (trim(spec%secondary_emission_model))
      case ('whipple')
        if (spec%secondary_emission_yield_max <= 0.0_dp .or. spec%secondary_emission_energy_max_ev <= 0.0_dp) then
          call stop_config_error('whipple SEE requires '//species_config_path(species_idx)// &
                                 '.surface_interaction.secondary_emission.yield_max > 0 and '// &
                                 species_config_path(species_idx)//'.surface_interaction.secondary_emission.energy_max_ev > 0.')
        end if
      case ('constant')
        if (spec%secondary_emission_const_yield < 0.0_dp) then
          call stop_config_error('constant SEE requires '//species_config_path(species_idx)// &
                                 '.surface_interaction.secondary_emission.const_yield >= 0.')
        end if
      case default
        call stop_config_error(species_config_path(species_idx)// &
                               '.surface_interaction.secondary_emission.model must be "whipple" or "constant".')
      end select
      spec%secondary_emission_target_idx = 0_i32
      do j = 1, cfg%n_particle_species
        if (len_trim(spec%secondary_emission_target_species) == 0) exit
        if (.not. cfg%particle_species(j)%enabled) cycle
        if (trim(spec%secondary_emission_target_species) /= trim(cfg%particle_species(j)%species_key)) cycle
        if (trim(cfg%particle_species(j)%source_mode) /= 'secondary_emission') exit
        spec%secondary_emission_target_idx = int(j, i32)
        exit
      end do
      if (spec%secondary_emission_target_idx == 0_i32) then
        call stop_config_error(species_config_path(species_idx)// &
                               '.surface_interaction.secondary_emission.target_species must name an enabled '// &
                               'secondary_emission species.')
      end if
    end associate
  end subroutine resolve_secondary_emission_target

  !> volume_seed の体積生成と reservoir 流入が共有する粒子属性と生成領域を検査する。
  !! slot 生成はこれらを再検査せずに使うので、native 入口ではここが唯一の所有者になる。
  !! 点状・平面状の生成領域（`pos_low == pos_high` の成分）と温度 0 は有効な入力として受け入れる。
  subroutine validate_volume_seed_species(spec, species_idx, has_boundary_inflow)
    type(particle_species_spec), intent(in) :: spec
    integer, intent(in) :: species_idx
    logical, intent(in) :: has_boundary_inflow
    real(dp) :: temperature_k

    if (.not. ieee_is_finite(spec%m_particle) .or. spec%m_particle <= 0.0_dp) then
      call stop_config_error(species_config_path(species_idx)//'.mass_kg must be finite and > 0.')
    end if
    temperature_k = species_temperature_k(spec)
    if (.not. ieee_is_finite(temperature_k) .or. temperature_k < 0.0_dp) then
      call stop_config_error(species_config_path(species_idx)//'.distribution must resolve to a finite temperature >= 0.')
    end if
    if (any(spec%pos_high < spec%pos_low)) then
      call stop_config_error(species_config_path(species_idx)//'.source.pos_high must be >= '// &
                             species_config_path(species_idx)//'.source.pos_low in every component.')
    end if
    ! target_macro_particles_per_batch は流入 flux から重みを導くので、reservoir 流入面がないと意味を持たない。
    if (spec%has_target_macro_particles_per_batch .and. .not. has_boundary_inflow) then
      call stop_config_error(species_config_path(species_idx)//'.sampling.target_macro_particles_per_batch requires '// &
                             species_config_path(species_idx)//'.inflow; volume_seed uses '// &
                             species_config_path(species_idx)//'.sampling.weight.')
    end if
    ! 流入重みは、解決する resolve_boundary_inflow_weight が検査する。
    if (has_boundary_inflow .and. spec%has_target_macro_particles_per_batch) return
    if (.not. ieee_is_finite(spec%w_particle) .or. spec%w_particle <= 0.0_dp) then
      call stop_config_error(species_config_path(species_idx)//'.sampling.weight must be finite and > 0.')
    end if
  end subroutine validate_volume_seed_species

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
          error stop 'target surface currents require '//species_config_path(species_idx)//'.charging.closure="fixed_current".'
        end if
      case ('fixed_current')
        if (.not. ieee_is_finite(cfg%sim%batch_duration) .or. cfg%sim%batch_duration <= 0.0_dp) then
          error stop 'run.batch.duration_s must be > 0 for fixed_current.'
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
            error stop species_config_path(species_idx)// &
              '.charging.target_absorbed_current_a must be finite and have the same sign as '// &
              species_config_path(species_idx)//'.charge_c.'
          end if
        end if
        if (spec%has_target_emission_current_a) then
          if (trim(spec%source_mode) /= 'photo_raycast' .or. .not. spec%deposit_opposite_charge_on_emit .or. &
              .not. ieee_is_finite(spec%target_emission_current_a) .or. &
              spec%target_emission_current_a*spec%q_particle > 0.0_dp) then
            error stop species_config_path(species_idx)// &
              '.charging.target_emission_current_a requires photo_raycast, an opposite-charge emission deposit, '// &
              'and current sign opposite to '//species_config_path(species_idx)//'.charge_c.'
          end if
        end if
      case ('neutral_return')
        if (.not. cfg%sim%use_box) error stop 'neutral_return requires a finite [domain].'
        inject_action = particle_boundary_action_for_face(effective_low, effective_high, spec%inject_face)
        if (trim(spec%source_mode) /= 'photo_raycast' .or. spec%q_particle >= 0.0_dp .or. &
            .not. spec%deposit_opposite_charge_on_emit .or. &
            (inject_action /= bc_reflect .and. inject_action /= bc_redistributed_reflect)) then
          error stop 'neutral_return requires a negative photo_raycast species, '// &
            species_config_path(species_idx)//'.source.deposit_opposite_charge_on_emit=true, '// &
            'and a reflecting action on '//species_config_path(species_idx)//'.source.inject_face.'
        end if
      case default
        error stop species_config_path(species_idx)//'.charging.closure must be "explicit", "fixed_current", or "neutral_return".'
      end select
    end associate
  end subroutine validate_surface_charge_closure

end submodule bem_app_config_parser_preflight_particles
