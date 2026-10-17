!> 外部モデルから species 別の固定表面電流を解決する。
module bem_surface_current_model
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann, pi, qe
  use bem_app_config_types, only: app_config
  use bem_irradiation_cell_contract, only: &
    irradiation_cell_contract_type, resolve_irradiation_cell_contract, irradiation_cell_ok
  use bem_surface_closure_contract, only: surface_closure_contract_type
  use bem_config_helpers, only: species_number_density_m3, species_temperature_k
  use sheath_model, only: plasma_input, sheath_solver, sheath_matching_conditions, sheath_matching_iteration, &
                          emission_source, maxwellian_emission, observed_outflow_source, SHEATH_OK
  use bem_outflow_spectrum, only: outflow_spectrum_bins, outflow_spectrum_edges_ev
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  type, extends(surface_closure_contract_type), public :: surface_current_model_result_type
    character(len=32) :: model = 'none'
    character(len=1) :: zhao_branch = ' '
    integer(i32) :: electron_species_idx = 0_i32
    integer(i32) :: ion_species_idx = 0_i32
    integer(i32) :: photoelectron_species_idx = 0_i32
    logical :: photoelectron_active = .false.
    integer(i32) :: secondary_electron_species_idx = 0_i32
    real(dp) :: secondary_share = 0.0_dp
    real(dp) :: reference_area_m2 = 0.0_dp
    real(dp) :: phi0_v = 0.0_dp
    real(dp) :: phi_m_v = 0.0_dp
    !> 上流に接する E^2<0 の電位幅 [V] と、その内側の場の不足 sqrt(-min E^2) [V/m]。
    !! 電子drift時に許容幅または場の下限で受理したA/C根だけが非ゼロ。
    real(dp) :: upstream_negative_band_v = 0.0_dp
    real(dp) :: upstream_field_deficit_v_m = 0.0_dp
    real(dp) :: ambient_electron_density_m3 = 0.0_dp
    real(dp) :: electron_current_density_a_m2 = 0.0_dp
    real(dp) :: ion_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_emission_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_escape_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_return_current_density_a_m2 = 0.0_dp
    real(dp) :: net_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_budget_residual_current_density_a_m2 = 0.0_dp
    real(dp) :: surface_budget_residual_current_density_a_m2 = 0.0_dp
    !> 外部Zhaoシースが壁放出として見るPE束 [1/(m2 s)] と平均法線エネルギー [eV]。
    real(dp) :: outer_photoelectron_flux_m2_s = 0.0_dp
    real(dp) :: outer_photoelectron_mean_energy_ev = 0.0_dp
    !> 外部源に使ったH外向きスペクトル [1/(m2 s)]、bin別。全て0は表面放出の半Maxwell源。
    real(dp) :: outer_photoelectron_spectrum_m2_s(outflow_spectrum_bins) = 0.0_dp
    character(len=32) :: kinetic_contract = 'none'
  end type surface_current_model_result_type

  public :: evaluate_surface_current_model
  public :: evaluate_surface_closure
  public :: solve_zhao_outflow_closure, refresh_zhao_outflow_closure, split_emission_escape_target

contains

  !> モデル固有の診断値をシミュレータへ漏らさず、境界契約だけを返す。
  subroutine evaluate_surface_closure(app, contract, cell)
    type(app_config), intent(in) :: app
    type(surface_closure_contract_type), intent(out) :: contract
    type(irradiation_cell_contract_type), intent(in), optional :: cell
    type(surface_current_model_result_type) :: detailed_result

    call evaluate_surface_current_model(app, detailed_result, cell)
    contract = detailed_result%surface_closure_contract_type
  end subroutine evaluate_surface_closure

  !> 設定されたmodelをdispatchし、固定電流closure用のtarget配列を返す。
  !! zhao_stationary で outer を渡すと、後の更新に使う外部解を返す。
  subroutine evaluate_surface_current_model(app, result, cell, outer)
    type(app_config), intent(in) :: app
    type(surface_current_model_result_type), intent(out) :: result
    type(irradiation_cell_contract_type), intent(in), optional :: cell
    type(sheath_matching_iteration), intent(out), optional :: outer
    type(sheath_matching_iteration) :: local_outer
    type(irradiation_cell_contract_type) :: active_cell
    integer(i32) :: status
    character(len=256) :: message

    if (present(cell)) then
      active_cell = cell
    else
      call resolve_irradiation_cell_contract(app, active_cell, status, message)
      if (status /= irradiation_cell_ok) error stop trim(message)
    end if
    call allocate_closure_channels(app, result)
    result%model = trim(lower_ascii(app%surface_current%model))

    select case (trim(result%model))
    case ('none')
      return
    case ('zhao_stationary')
      call evaluate_zhao_stationary_current(app, result, active_cell, local_outer)
      if (present(outer)) outer = local_outer
    case default
      error stop 'Unknown surface current model dispatch.'
    end select
  end subroutine evaluate_surface_current_model

  !> species別channelを未使用状態で確保する。
  subroutine allocate_closure_channels(app, result)
    type(app_config), intent(in) :: app
    type(surface_current_model_result_type), intent(inout) :: result

    allocate ( &
      result%has_absorbed_target(app%n_particle_species), &
      result%has_emission_target(app%n_particle_species), &
      result%has_escape_target(app%n_particle_species), &
      result%has_inflow_kinetic_map(app%n_particle_species), &
      result%has_outflow_kinetic_barrier(app%n_particle_species), &
      result%has_inflow_number_flux(app%n_particle_species), &
      result%absorbed_current_a(app%n_particle_species), &
      result%emission_current_a(app%n_particle_species), &
      result%escaped_particle_current_a(app%n_particle_species), &
      result%inflow_reservoir_potential_v(app%n_particle_species), &
      result%inflow_access_potential_v(app%n_particle_species), &
      result%inflow_reservoir_density_m3(app%n_particle_species), &
      result%inflow_kinetic_face(app%n_particle_species), &
      result%outflow_barrier_potential_v(app%n_particle_species), &
      result%outflow_barrier_face(app%n_particle_species), &
      result%inflow_number_flux_m2_s(app%n_particle_species) &
      )
    result%has_absorbed_target = .false.
    result%has_emission_target = .false.
    result%has_escape_target = .false.
    result%has_inflow_kinetic_map = .false.
    result%has_outflow_kinetic_barrier = .false.
    result%has_inflow_number_flux = .false.
    result%absorbed_current_a = 0.0_dp
    result%emission_current_a = 0.0_dp
    result%escaped_particle_current_a = 0.0_dp
    result%inflow_reservoir_potential_v = 0.0_dp
    result%inflow_access_potential_v = 0.0_dp
    result%inflow_reservoir_density_m3 = 0.0_dp
    result%inflow_kinetic_face = 0_i32
    result%outflow_barrier_potential_v = 0.0_dp
    result%outflow_barrier_face = 0_i32
    result%inflow_number_flux_m2_s = 0.0_dp
  end subroutine allocate_closure_channels

  !> 設定の表面放出を外部源として零電流根を解く。PEなしではType Cだけが成立する。
  subroutine evaluate_zhao_stationary_current(app, result, cell, outer)
    type(app_config), intent(in) :: app
    type(surface_current_model_result_type), intent(inout) :: result
    type(irradiation_cell_contract_type), intent(in) :: cell
    type(sheath_matching_iteration), intent(out) :: outer
    character(len=3) :: order
    character(len=256) :: message
    real(dp) :: source_density_m3, source_temperature_ev
    integer(i32) :: status
    logical :: success

    call surface_photoelectron_source(app, source_density_m3, source_temperature_ev)
    order = elevation_order(app%surface_current%solar_elevation_deg)
    if (app%surface_current%photoelectron_source_scale <= 0.0_dp) then
      ! 外部根の更新ではSEE流出で他のbranchへ移れるよう、autoのままCを先に試す。
      select case (trim(lower_ascii(app%surface_current%zhao_branch)))
      case ('auto', 'c')
        order = 'CAB'
      case default
        error stop 'sheath.photoelectrons.source_scale=0 requires sheath.zhao.branch="auto" or "c".'
      end select
    end if
    call start_outer(app, cell, maxwellian_emission(source_density_m3, source_temperature_ev), order, outer, status, message)
    if (status /= SHEATH_OK) error stop 'Zhao stationary surface-current root solve failed: '//trim(message)
    call finish_zhao_closure(app, outer%conditions(), result, success, message, cell=cell)
    if (.not. success) error stop 'Zhao stationary '//trim(message)
    call record_outer_source(app, result, source_density_m3, source_temperature_ev)
  end subroutine evaluate_zhao_stationary_current

  !> 表面放出はそのままに、matching planeで観測したPE流出を外部Zhaoシースの放出源として零電流根を解き直す。
  !!
  !! 周期セルは外部シースから見て一様な壁なので、外部の1-D解は壁から出るPEとしてH通過流出だけを見る。
  !! 合算した流出スペクトルを使い、省略時だけ half-Maxwellian とする。前回の外部状態は持たず、
  !! preferred_branch を先に試す。outer を渡すと、後の更新に使う外部解を返す。
  subroutine solve_zhao_outflow_closure( &
    app, outer_flux_m2_s, outer_mean_energy_ev, preferred_branch, result, success, message, &
    outer_spectrum_m2_s, secondary_share, cell, outer &
    )
    type(app_config), intent(in) :: app
    real(dp), intent(in) :: outer_flux_m2_s, outer_mean_energy_ev
    character(len=1), intent(in) :: preferred_branch
    type(surface_current_model_result_type), intent(out) :: result
    logical, intent(out) :: success
    character(len=*), intent(out) :: message
    !> H外向きスペクトル [1/(m2 s)]、bin別。省略時は束と平均法線エネルギーの半Maxwell源とする。
    real(dp), intent(in), optional :: outer_spectrum_m2_s(outflow_spectrum_bins)
    real(dp), intent(in), optional :: secondary_share
    type(irradiation_cell_contract_type), intent(in), optional :: cell
    type(sheath_matching_iteration), intent(out), optional :: outer
    type(sheath_matching_iteration) :: local_outer
    type(irradiation_cell_contract_type) :: active_cell
    type(emission_source) :: source
    integer(i32) :: status
    character(len=256) :: root_message

    success = .false.
    message = ''
    if (present(cell)) then
      active_cell = cell
    else
      call resolve_irradiation_cell_contract(app, active_cell, status, message)
      if (status /= irradiation_cell_ok) return
    end if
    call outer_emission(app, outer_flux_m2_s, outer_mean_energy_ev, source, success, message, outer_spectrum_m2_s)
    if (.not. success) return
    ! 外部根の解き直しでは表面の太陽高度ではなく、前回のbranchを優先して連続性を保つ。
    call start_outer(app, active_cell, source, preferred_order(preferred_branch), local_outer, status, root_message)
    success = status == SHEATH_OK
    if (.not. success) then
      message = 'outer Zhao zero-current root was not found for the observed photoelectron outflow: '// &
                trim(root_message)
      return
    end if
    call finish_observed_closure( &
      app, local_outer, outer_flux_m2_s, outer_mean_energy_ev, result, success, message, &
      outer_spectrum_m2_s, secondary_share, active_cell &
      )
    if (present(outer)) outer = local_outer
  end subroutine solve_zhao_outflow_closure

  !> 外部解を観測したH外向き流出で更新する。前回の根を初期値にし、同じbranchを先に試す。
  !! 解けなければ outer は前回のまま残り、success=.false. を返す。
  subroutine refresh_zhao_outflow_closure( &
    app, outer, outer_flux_m2_s, outer_mean_energy_ev, outer_spectrum_m2_s, secondary_share, cell, result, &
    success, message &
    )
    type(app_config), intent(in) :: app
    type(sheath_matching_iteration), intent(inout) :: outer
    real(dp), intent(in) :: outer_flux_m2_s, outer_mean_energy_ev
    real(dp), intent(in) :: outer_spectrum_m2_s(outflow_spectrum_bins)
    real(dp), intent(in) :: secondary_share
    type(irradiation_cell_contract_type), intent(in) :: cell
    type(surface_current_model_result_type), intent(out) :: result
    logical, intent(out) :: success
    character(len=*), intent(out) :: message
    type(emission_source) :: source
    integer(i32) :: status
    character(len=256) :: root_message

    call outer_emission(app, outer_flux_m2_s, outer_mean_energy_ev, source, success, message, outer_spectrum_m2_s)
    if (.not. success) return
    call outer%update(source, status, root_message)
    success = status == SHEATH_OK
    if (.not. success) then
      message = 'outer Zhao zero-current root was not found for the observed photoelectron outflow: '// &
                trim(root_message)
      return
    end if
    call finish_observed_closure( &
      app, outer, outer_flux_m2_s, outer_mean_energy_ev, result, success, message, &
      outer_spectrum_m2_s, secondary_share, cell &
      )
  end subroutine refresh_zhao_outflow_closure

  !> 観測した束・平均法線エネルギー・スペクトルから外部の放出源を作る。
  subroutine outer_emission(app, outer_flux_m2_s, outer_mean_energy_ev, source, success, message, outer_spectrum_m2_s)
    type(app_config), intent(in) :: app
    real(dp), intent(in) :: outer_flux_m2_s, outer_mean_energy_ev
    type(emission_source), intent(out) :: source
    logical, intent(out) :: success
    character(len=*), intent(out) :: message
    real(dp), intent(in), optional :: outer_spectrum_m2_s(outflow_spectrum_bins)
    integer(i32) :: status

    success = .false.
    message = ''
    if (.not. all(ieee_is_finite([outer_flux_m2_s, outer_mean_energy_ev])) .or. &
        outer_flux_m2_s <= 0.0_dp .or. outer_mean_energy_ev <= 0.0_dp) then
      message = 'outer photoelectron source must have positive finite flux and mean normal energy.'
      return
    end if
    if (present(outer_spectrum_m2_s)) then
      call observed_outflow_source(outflow_spectrum_edges_ev(), outer_spectrum_m2_s, source, status, message)
      if (status /= SHEATH_OK) then
        message = 'outer photoelectron spectrum must be finite, nonnegative, and nonzero.'
        return
      end if
    else
      source = maxwellian_emission(2.0_dp*sqrt(pi)*outer_flux_m2_s/thermal_speed(app, outer_mean_energy_ev), &
                                   outer_mean_energy_ev)
    end if
    success = .true.
  end subroutine outer_emission

  !> 観測した外部源で解いた外部解から closure を作り、外部源の受領値を記録する。
  subroutine finish_observed_closure( &
    app, outer, outer_flux_m2_s, outer_mean_energy_ev, result, success, message, outer_spectrum_m2_s, &
    secondary_share, cell &
    )
    type(app_config), intent(in) :: app
    type(sheath_matching_iteration), intent(in) :: outer
    real(dp), intent(in) :: outer_flux_m2_s, outer_mean_energy_ev
    type(surface_current_model_result_type), intent(out) :: result
    logical, intent(out) :: success
    character(len=*), intent(out) :: message
    real(dp), intent(in), optional :: outer_spectrum_m2_s(outflow_spectrum_bins)
    real(dp), intent(in), optional :: secondary_share
    type(irradiation_cell_contract_type), intent(in) :: cell

    call allocate_closure_channels(app, result)
    result%model = 'zhao_stationary'
    call finish_zhao_closure(app, outer%conditions(), result, success, message, secondary_share, cell=cell)
    result%outer_photoelectron_flux_m2_s = outer_flux_m2_s
    result%outer_photoelectron_mean_energy_ev = outer_mean_energy_ev
    if (present(outer_spectrum_m2_s)) result%outer_photoelectron_spectrum_m2_s = outer_spectrum_m2_s
  end subroutine finish_observed_closure

  !> 設定の太陽高度・基準密度・倍率から、表面放出するhalf-MaxwellianのMaxwell規格化密度と温度を返す。
  subroutine surface_photoelectron_source(app, density_m3, temperature_ev)
    type(app_config), intent(in) :: app
    real(dp), intent(out) :: density_m3, temperature_ev
    integer(i32) :: electron_idx, photo_idx

    electron_idx = species_index(app, app%surface_current%electron_species)
    density_m3 = 0.0_dp
    ! PEなしでは源密度0のMaxwell源を渡す。温度は解に寄与しない。
    temperature_ev = species_temperature_k(app%particle_species(electron_idx))*k_boltzmann/qe
    if (app%surface_current%photoelectron_source_scale <= 0.0_dp) return
    photo_idx = species_index(app, app%surface_current%photoelectron_species)
    temperature_ev = species_temperature_k(app%particle_species(photo_idx))*k_boltzmann/qe
    density_m3 = app%surface_current%photoelectron_source_scale*app%surface_current%photoelectron_ref_density_m3* &
                 sin(app%surface_current%solar_elevation_deg*pi/180.0_dp)
  end subroutine surface_photoelectron_source

  !> 設定のbranchで外部根を解く。auto では order の順に試し、明示したbranchではそれだけを解く。
  subroutine start_outer(app, cell, emission, order, outer, status, message)
    type(app_config), intent(in) :: app
    type(irradiation_cell_contract_type), intent(in) :: cell
    type(emission_source), intent(in) :: emission
    character(len=3), intent(in) :: order
    type(sheath_matching_iteration), intent(out) :: outer
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=:), allocatable :: branch

    branch = trim(lower_ascii(app%surface_current%zhao_branch))
    if (branch == 'auto') then
      call outer%initialize(zhao_outer_plasma(app, cell, emission), zhao_solver(app), status, message, &
                            initial_order=order)
    else
      call outer%initialize(zhao_outer_plasma(app, cell, emission), zhao_solver(app), status, message, branch=branch)
    end if
  end subroutine start_outer

  !> 太陽高度20度未満ではC、A、B、それ以外はA、B、Cの順に試す。
  pure character(len=3) function elevation_order(solar_elevation_deg) result(order)
    real(dp), intent(in) :: solar_elevation_deg

    order = 'ABC'
    if (solar_elevation_deg < 20.0_dp) order = 'CAB'
  end function elevation_order

  !> preferred_branch を先頭にし、残りをA、B、Cの順に並べる。不明なbranchではA、B、C。
  pure character(len=3) function preferred_order(preferred_branch) result(order)
    character(len=1), intent(in) :: preferred_branch
    character(len=3), parameter :: branches = 'ABC'
    character(len=1) :: preferred
    integer :: i

    preferred = preferred_branch
    if (preferred >= 'a' .and. preferred <= 'z') preferred = achar(iachar(preferred) - 32)
    order = branches
    if (preferred == ' ' .or. index(branches, preferred) == 0) return
    order = preferred
    do i = 1, 3
      if (branches(i:i) /= preferred) order = trim(order)//branches(i:i)
    end do
  end function preferred_order

  !> 外部シースの上流状態。ionは ion species の温度と設定の圧力係数を持つ上流接続の流体ビーム
  !! (Sana & Mishra 2026)、電子は設定の内向きdriftを持つ上流Maxwell分布とする。太陽風では T_i が
  !! 入口エネルギーよりはるかに小さく、根への影響は 0.1 mV 以下である。
  function zhao_outer_plasma(app, cell, emission) result(plasma)
    type(app_config), intent(in) :: app
    type(irradiation_cell_contract_type), intent(in) :: cell
    type(emission_source), intent(in) :: emission
    type(plasma_input) :: plasma
    integer(i32) :: electron_idx, ion_idx

    electron_idx = species_index(app, app%surface_current%electron_species)
    ion_idx = species_index(app, app%surface_current%ion_species)
    plasma%ion_density_m3 = species_number_density_m3(app%particle_species(ion_idx))
    plasma%electron_temperature_ev = species_temperature_k(app%particle_species(electron_idx))*k_boltzmann/qe
    plasma%ion_temperature_ev = species_temperature_k(app%particle_species(ion_idx))*k_boltzmann/qe
    plasma%ion_pressure_factor = app%surface_current%zhao_ion_pressure_factor
    plasma%electron_drift_mps = -cell%outward_normal_velocity(app%particle_species(electron_idx)%drift_velocity)
    plasma%ion_entry_speed_mps = -cell%outward_normal_velocity(app%particle_species(ion_idx)%drift_velocity)
    plasma%ion_mass_kg = app%particle_species(ion_idx)%m_particle
    plasma%electron_mass_kg = app%particle_species(electron_idx)%m_particle
    plasma%emission = emission
  end function zhao_outer_plasma

  !> 内向きdriftのある電子で反射集団を持つA/Cは無限遠へ接続できず、
  !! upstream_band_tolerance と upstream_field_floor_v_m が共に 0 なら sheath-model が棄却する。
  function zhao_solver(app) result(solver)
    type(app_config), intent(in) :: app
    type(sheath_solver) :: solver

    solver%search%upstream_band_tolerance = app%surface_current%zhao_upstream_band_tolerance
    solver%search%upstream_field_floor_v_m = app%surface_current%zhao_upstream_field_floor_v_m
  end function zhao_solver

  !> Maxwell規格化密度と束の換算に使う、温度 [eV] の電子の熱速度 sqrt(2T/m) [m/s]。
  real(dp) function thermal_speed(app, temperature_ev) result(speed)
    type(app_config), intent(in) :: app
    real(dp), intent(in) :: temperature_ev

    speed = sqrt(2.0_dp*qe*temperature_ev/ &
                 app%particle_species(species_index(app, app%surface_current%electron_species))%m_particle)
  end function thermal_speed

  !> 静的closureの外部源（設定の表面放出）を受領値として記録する。PE・SEEのない closure では 0。
  subroutine record_outer_source(app, result, source_density_m3, source_temperature_ev)
    type(app_config), intent(in) :: app
    type(surface_current_model_result_type), intent(inout) :: result
    real(dp), intent(in) :: source_density_m3, source_temperature_ev

    result%outer_photoelectron_flux_m2_s = 0.0_dp
    result%outer_photoelectron_mean_energy_ev = 0.0_dp
    if (.not. (result%photoelectron_active .or. result%secondary_electron_species_idx > 0_i32)) return
    result%outer_photoelectron_flux_m2_s = source_density_m3*thermal_speed(app, source_temperature_ev)/(2.0_dp*sqrt(pi))
    result%outer_photoelectron_mean_energy_ev = source_temperature_ev
  end subroutine record_outer_source

  !> sheath-modelの面上条件から、表面放出を別に与えて species 別の固定電流targetと境界写像を作る。
  !!
  !! 外部シースが見る放出源と表面放出の差は、周期セル内で再吸収されたPEであり、return target に含まれる。
  subroutine finish_zhao_closure(app, plane, result, success, message, secondary_share, cell)
    type(app_config), intent(in) :: app
    type(sheath_matching_conditions), intent(in) :: plane
    type(surface_current_model_result_type), intent(inout) :: result
    logical, intent(out) :: success
    character(len=*), intent(out) :: message
    real(dp), intent(in), optional :: secondary_share
    type(irradiation_cell_contract_type), intent(in) :: cell
    integer(i32) :: electron_idx, ion_idx, photo_idx, secondary_idx
    real(dp) :: area, budget_scale, budget_tolerance
    real(dp) :: surface_density_m3, surface_temperature_ev, emission_current_density
    logical :: photoelectron_active, secondary_active

    success = .false.
    message = ''
    electron_idx = species_index(app, app%surface_current%electron_species)
    ion_idx = species_index(app, app%surface_current%ion_species)
    photoelectron_active = app%surface_current%photoelectron_source_scale > 0.0_dp
    secondary_active = len_trim(app%surface_current%secondary_electron_species) > 0
    secondary_idx = 0_i32
    if (secondary_active) secondary_idx = species_index(app, app%surface_current%secondary_electron_species)
    result%secondary_share = 0.0_dp
    if (secondary_active .and. .not. photoelectron_active) result%secondary_share = 1.0_dp
    if (present(secondary_share)) result%secondary_share = secondary_share
    photo_idx = 0_i32
    if (photoelectron_active) photo_idx = species_index(app, app%surface_current%photoelectron_species)
    area = cell%wall_area(app%sim%box_min, app%sim%box_max)
    if (app%surface_current%has_reference_area_m2) area = app%surface_current%reference_area_m2
    if (.not. ieee_is_finite(area) .or. area <= 0.0_dp) then
      message = 'surface-current reference area must be finite and positive.'
      return
    end if

    result%zhao_branch = plane%branch
    result%phi0_v = plane%plane_potential_v
    result%phi_m_v = plane%plane_potential_v
    if (plane%branch == 'A') result%phi_m_v = plane%electron_path_minimum_v
    result%upstream_negative_band_v = plane%upstream_negative_band_v
    result%upstream_field_deficit_v_m = plane%upstream_field_deficit_v_m
    result%ambient_electron_density_m3 = plane%electron_reservoir_density_m3
    call surface_photoelectron_source(app, surface_density_m3, surface_temperature_ev)
    emission_current_density = qe*surface_density_m3*thermal_speed(app, surface_temperature_ev)/(2.0_dp*sqrt(pi))
    result%electron_current_density_a_m2 = -qe*plane%electron_inward_flux_m2_s
    result%ion_current_density_a_m2 = qe*plane%ion_inward_flux_m2_s
    result%photoelectron_emission_current_density_a_m2 = emission_current_density
    result%photoelectron_escape_current_density_a_m2 = 0.0_dp
    if (photoelectron_active .or. secondary_active) then
      result%photoelectron_escape_current_density_a_m2 = qe*plane%emitted_escape_flux_m2_s
    end if
    result%photoelectron_return_current_density_a_m2 = &
      result%photoelectron_escape_current_density_a_m2 - result%photoelectron_emission_current_density_a_m2
    result%net_current_density_a_m2 = result%electron_current_density_a_m2 + &
                                      result%ion_current_density_a_m2 + &
                                      result%photoelectron_escape_current_density_a_m2
    result%photoelectron_budget_residual_current_density_a_m2 = &
      result%photoelectron_emission_current_density_a_m2 + &
      result%photoelectron_return_current_density_a_m2 - &
      result%photoelectron_escape_current_density_a_m2
    result%surface_budget_residual_current_density_a_m2 = &
      result%electron_current_density_a_m2 + result%ion_current_density_a_m2 + &
      result%photoelectron_emission_current_density_a_m2 + &
      result%photoelectron_return_current_density_a_m2
    if (.not. all(ieee_is_finite([ &
                                 result%electron_current_density_a_m2, result%ion_current_density_a_m2, &
                                 result%photoelectron_emission_current_density_a_m2, &
                                 result%photoelectron_escape_current_density_a_m2, &
                                 result%photoelectron_return_current_density_a_m2, result%net_current_density_a_m2, &
                                 result%photoelectron_budget_residual_current_density_a_m2, &
                                 result%surface_budget_residual_current_density_a_m2 &
                                 ]))) then
      message = 'surface-current evaluation produced non-finite currents.'
      return
    end if
    if (.not. ieee_is_finite(result%ambient_electron_density_m3) .or. &
        result%ambient_electron_density_m3 <= 0.0_dp) then
      message = 'surface-current root produced an invalid upstream electron density.'
      return
    end if
    if (result%electron_current_density_a_m2 >= 0.0_dp .or. &
        result%ion_current_density_a_m2 <= 0.0_dp .or. &
        (.not. secondary_active .and. result%photoelectron_return_current_density_a_m2 > 0.0_dp) .or. &
        result%photoelectron_escape_current_density_a_m2 < 0.0_dp) then
      message = 'surface-current evaluation produced invalid channel signs.'
      return
    end if
    if (photoelectron_active) then
      if (result%photoelectron_emission_current_density_a_m2 <= 0.0_dp) then
        message = 'photoelectron closure requires a positive emission current.'
        return
      end if
    else if (.not. secondary_active .and. any([ &
                                              result%photoelectron_emission_current_density_a_m2, &
                                              result%photoelectron_escape_current_density_a_m2, &
                                              result%photoelectron_return_current_density_a_m2 &
                                              ] /= 0.0_dp)) then
      message = 'no-photoelectron closure produced a nonzero photoelectron current.'
      return
    end if
    budget_scale = max( &
                   abs(result%electron_current_density_a_m2), abs(result%ion_current_density_a_m2), &
                   abs(result%photoelectron_emission_current_density_a_m2), &
                   abs(result%photoelectron_escape_current_density_a_m2), &
                   abs(result%photoelectron_return_current_density_a_m2), tiny(1.0_dp) &
                   )
    budget_tolerance = sqrt(epsilon(1.0_dp))*budget_scale
    if (abs(result%photoelectron_budget_residual_current_density_a_m2) > budget_tolerance) then
      message = 'PE current budget does not close.'
      return
    end if
    if (abs(result%surface_budget_residual_current_density_a_m2) > budget_tolerance) then
      message = 'surface current budget does not close.'
      return
    end if

    result%active = .true.
    result%reference_area_m2 = area
    result%electron_species_idx = electron_idx
    result%ion_species_idx = ion_idx
    result%photoelectron_species_idx = photo_idx
    result%photoelectron_active = photoelectron_active
    result%secondary_electron_species_idx = secondary_idx
    result%has_absorbed_target([electron_idx, ion_idx]) = .true.
    result%absorbed_current_a(electron_idx) = &
      checked_area_current(area, result%electron_current_density_a_m2)
    result%absorbed_current_a(ion_idx) = checked_area_current(area, result%ion_current_density_a_m2)
    if (photoelectron_active) then
      result%has_absorbed_target(photo_idx) = .true.
      result%has_emission_target(photo_idx) = .true.
      result%has_escape_target(photo_idx) = .true.
      result%emission_current_a(photo_idx) = &
        checked_area_current(area, result%photoelectron_emission_current_density_a_m2)
    end if

    call split_emission_escape_target(result, result%secondary_share)

    ! Zhao の1-D外部シースを、z-high interfaceに対するkinetic boundary mapへ縮約する。
    ! 経路の電位は面の電位とは別に越える電位で、Type Aの電子はphi_m、それ以外は上流の0 V。
    result%kinetic_contract = 'zhao_barrier_v1'
    result%has_inflow_kinetic_map([electron_idx, ion_idx]) = .true.
    result%inflow_reservoir_potential_v([electron_idx, ion_idx]) = 0.0_dp
    ! 上流電子のMaxwellian密度は、壁で吸われて戻らない高速電子を見込んで無限遠の準中性から決まる。
    ! 設定の太陽風密度をそのまま電子源に使うと流入fluxが根の電子電流からずれる。
    result%inflow_reservoir_density_m3(electron_idx) = plane%electron_reservoir_density_m3
    result%inflow_reservoir_density_m3(ion_idx) = plane%ion_reservoir_density_m3
    result%inflow_access_potential_v(electron_idx) = plane%electron_path_minimum_v
    result%inflow_access_potential_v(ion_idx) = plane%ion_path_maximum_v
    result%inflow_kinetic_face([electron_idx, ion_idx]) = cell%wall_face
    result%has_outflow_kinetic_barrier([electron_idx, ion_idx]) = .true.
    result%outflow_barrier_potential_v(electron_idx) = plane%electron_path_minimum_v
    result%outflow_barrier_potential_v(ion_idx) = plane%ion_path_maximum_v
    result%outflow_barrier_face([electron_idx, ion_idx]) = cell%wall_face
    if (photoelectron_active) then
      result%has_outflow_kinetic_barrier(photo_idx) = .true.
      result%outflow_barrier_potential_v(photo_idx) = plane%electron_path_minimum_v
      result%outflow_barrier_face(photo_idx) = cell%wall_face
    end if
    if (secondary_active) then
      result%has_outflow_kinetic_barrier(secondary_idx) = .true.
      result%outflow_barrier_potential_v(secondary_idx) = plane%electron_path_minimum_v
      result%outflow_barrier_face(secondary_idx) = cell%wall_face
    end if
    ! 外部障壁と流入写像は上流0 Vからの電位差なので、z-high面の平均電位を壁電位phi0へ固定する。
    ! 外部での横移動はセル幅よりはるかに大きいので、x/y周期セルでは戻り位置を面内一様にする。
    result%has_plane_gauge = .true.
    result%plane_gauge_potential_v = plane%plane_potential_v
    ! e_top_from_sheath は外部根の面直上の法線電場を零モードへ渡す。
    result%top_field_v_m = plane%plane_field_v_m
    result%outer_return_cell_uniform = cell%wall_tangents_periodic()
    success = .true.
  end subroutine finish_zhao_closure

  !> Phase 1 uses the latest observed window even when the outer root is retained.
  subroutine split_emission_escape_target(result, secondary_share)
    type(surface_current_model_result_type), intent(inout) :: result
    real(dp), value :: secondary_share
    integer(i32) :: photo_idx
    real(dp) :: photo_escape_current

    result%secondary_share = secondary_share
    if (.not. result%photoelectron_active) return
    photo_idx = result%photoelectron_species_idx
    photo_escape_current = (1.0_dp - secondary_share)*result%photoelectron_escape_current_density_a_m2
    result%absorbed_current_a(photo_idx) = checked_area_current( &
                                           result%reference_area_m2, &
                                           photo_escape_current - result%photoelectron_emission_current_density_a_m2 &
                                           )
    ! Escaping particle current has the opposite sign to surface charging current.
    result%escaped_particle_current_a(photo_idx) = -checked_area_current(result%reference_area_m2, photo_escape_current)
  end subroutine split_emission_escape_target

  integer(i32) function species_index(app, species_key) result(index_value)
    type(app_config), intent(in) :: app
    character(len=*), intent(in) :: species_key
    integer(i32) :: idx

    index_value = 0_i32
    do idx = 1_i32, app%n_particle_species
      if (.not. app%particle_species(idx)%enabled) cycle
      if (trim(app%particle_species(idx)%species_key) /= trim(species_key)) cycle
      index_value = idx
      return
    end do
    error stop 'Surface current model species resolution failed: '//trim(species_key)
  end function species_index

  real(dp) function checked_area_current(area_m2, current_density_a_m2) result(current_a)
    real(dp), intent(in) :: area_m2, current_density_a_m2

    if (.not. all(ieee_is_finite([area_m2, current_density_a_m2])) .or. area_m2 <= 0.0_dp) then
      error stop 'Zhao stationary surface-current target conversion received invalid input.'
    end if
    if (area_m2 > 1.0_dp .and. abs(current_density_a_m2) > huge(current_a)/area_m2) then
      error stop 'Zhao stationary surface-current target conversion overflowed.'
    end if
    current_a = area_m2*current_density_a_m2
    if (.not. ieee_is_finite(current_a)) then
      error stop 'Zhao stationary surface-current target conversion produced a non-finite current.'
    end if
    if (current_density_a_m2 /= 0.0_dp .and. current_a == 0.0_dp) then
      error stop 'Zhao stationary surface-current target conversion underflowed.'
    end if
  end function checked_area_current

end module bem_surface_current_model
