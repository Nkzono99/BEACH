!> 外部モデルから species 別の固定表面電流を解決する。
module bem_surface_current_model
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann, pi, qe
  use bem_app_config_types, only: app_config
  use bem_types, only: bc_periodic
  use bem_surface_closure_contract, only: surface_closure_contract_type
  use bem_config_helpers, only: species_number_density_m3, species_temperature_k
  use sheath_model, only: fixed_entry_equilibrium_input, sheath_equilibrium_result, sheath_solver, &
                          maxwellian_photoelectrons, SHEATH_OK
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
    real(dp) :: reference_area_m2 = 0.0_dp
    real(dp) :: phi0_v = 0.0_dp
    real(dp) :: phi_m_v = 0.0_dp
    real(dp) :: ambient_electron_density_m3 = 0.0_dp
    real(dp) :: electron_current_density_a_m2 = 0.0_dp
    real(dp) :: ion_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_emission_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_escape_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_return_current_density_a_m2 = 0.0_dp
    real(dp) :: net_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_budget_residual_current_density_a_m2 = 0.0_dp
    real(dp) :: surface_budget_residual_current_density_a_m2 = 0.0_dp
    !> 外部Zhaoシースが壁放出として見るPE束 [1/(m2 s)] とhalf-Maxwellian温度 [eV]。
    real(dp) :: outer_photoelectron_flux_m2_s = 0.0_dp
    real(dp) :: outer_photoelectron_mean_energy_ev = 0.0_dp
    character(len=32) :: kinetic_contract = 'none'
  end type surface_current_model_result_type

  public :: evaluate_surface_current_model
  public :: evaluate_surface_closure
  public :: solve_zhao_outflow_closure

contains

  !> モデル固有の診断値をシミュレータへ漏らさず、境界契約だけを返す。
  subroutine evaluate_surface_closure(app, contract)
    type(app_config), intent(in) :: app
    type(surface_closure_contract_type), intent(out) :: contract
    type(surface_current_model_result_type) :: detailed_result

    call evaluate_surface_current_model(app, detailed_result)
    contract = detailed_result%surface_closure_contract_type
  end subroutine evaluate_surface_closure

  !> 設定されたmodelをdispatchし、固定電流closure用のtarget配列を返す。
  subroutine evaluate_surface_current_model(app, result)
    type(app_config), intent(in) :: app
    type(surface_current_model_result_type), intent(out) :: result

    call allocate_closure_channels(app, result)
    result%model = trim(lower_ascii(app%surface_current%model))

    select case (trim(result%model))
    case ('none')
      return
    case ('zhao_stationary')
      call evaluate_zhao_stationary_current(app, result)
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
  subroutine evaluate_zhao_stationary_current(app, result)
    type(app_config), intent(in) :: app
    type(surface_current_model_result_type), intent(inout) :: result
    type(sheath_equilibrium_result) :: root
    character(len=1), allocatable :: order(:)
    character(len=512) :: message
    real(dp) :: source_density_m3, source_temperature_ev
    logical :: success

    call surface_photoelectron_source(app, source_density_m3, source_temperature_ev)
    if (app%surface_current%photoelectron_source_scale > 0.0_dp) then
      call outer_branch_order(app, ' ', app%surface_current%solar_elevation_deg, order)
    else
      select case (trim(lower_ascii(app%surface_current%zhao_branch)))
      case ('auto', 'c')
        order = ['C']
      case default
        error stop 'photoelectron_source_scale=0 requires surface_current_model.zhao_branch="auto" or "c".'
      end select
    end if
    call solve_outer_root(app, source_density_m3, source_temperature_ev, order, root, success, message)
    if (.not. success) error stop 'Zhao stationary surface-current root solve failed: '//trim(message)
    call finish_zhao_closure(app, root, source_density_m3, source_temperature_ev, result, success, message)
    if (.not. success) error stop 'Zhao stationary '//trim(message)
  end subroutine evaluate_zhao_stationary_current

  !> 表面放出はそのままに、matching planeで観測したPE流出を外部Zhaoシースの放出源として零電流根を解き直す。
  !!
  !! 周期セルは外部シースから見て一様な壁なので、外部の1-D解は壁から出るPEとしてH通過流出だけを見る。
  !! 流出は同じ束と平均法線エネルギーを持つhalf-Maxwellianへ縮約する。previousは同じbranchの初期値に使う。
  subroutine solve_zhao_outflow_closure( &
    app, outer_flux_m2_s, outer_mean_energy_ev, preferred_branch, result, success, message, previous &
    )
    type(app_config), intent(in) :: app
    real(dp), intent(in) :: outer_flux_m2_s, outer_mean_energy_ev
    character(len=1), intent(in) :: preferred_branch
    type(surface_current_model_result_type), intent(out) :: result
    logical, intent(out) :: success
    character(len=*), intent(out) :: message
    type(surface_current_model_result_type), intent(in), optional :: previous
    type(sheath_equilibrium_result) :: root
    character(len=1), allocatable :: order(:)
    character(len=512) :: root_message
    real(dp) :: outer_density_m3

    success = .false.
    message = ''
    if (.not. all(ieee_is_finite([outer_flux_m2_s, outer_mean_energy_ev])) .or. &
        outer_flux_m2_s <= 0.0_dp .or. outer_mean_energy_ev <= 0.0_dp) then
      message = 'outer photoelectron source must have positive finite flux and mean normal energy.'
      return
    end if
    call allocate_closure_channels(app, result)
    result%model = 'zhao_stationary'
    outer_density_m3 = 2.0_dp*sqrt(pi)*outer_flux_m2_s/ &
                       thermal_speed(app, outer_mean_energy_ev)
    ! 外部根の更新では表面の太陽高度ではなく、前回のbranchを優先して連続性を保つ。
    call outer_branch_order(app, preferred_branch, 90.0_dp, order)
    if (present(previous)) then
      call solve_outer_root( &
        app, outer_density_m3, outer_mean_energy_ev, order, root, success, root_message, previous=previous &
        )
    else
      call solve_outer_root(app, outer_density_m3, outer_mean_energy_ev, order, root, success, root_message)
    end if
    if (.not. success) then
      message = 'outer Zhao zero-current root was not found for the observed photoelectron outflow: '// &
                trim(root_message)
      return
    end if
    call finish_zhao_closure(app, root, outer_density_m3, outer_mean_energy_ev, result, success, message)
    result%outer_photoelectron_flux_m2_s = outer_flux_m2_s
    result%outer_photoelectron_mean_energy_ev = outer_mean_energy_ev
  end subroutine solve_zhao_outflow_closure

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

  !> 明示branchはそれだけ、autoは優先branch、続いて高度20度未満ならC/A/B、他はA/B/Cの順に試す。
  subroutine outer_branch_order(app, preferred_branch, solar_elevation_deg, order)
    type(app_config), intent(in) :: app
    character(len=1), intent(in) :: preferred_branch
    real(dp), intent(in) :: solar_elevation_deg
    character(len=1), allocatable, intent(out) :: order(:)
    character(len=1) :: base(3)
    character(len=1) :: preferred
    integer :: i

    if (trim(lower_ascii(app%surface_current%zhao_branch)) /= 'auto') then
      order = [upper_branch(app%surface_current%zhao_branch(1:1))]
      return
    end if
    base = ['A', 'B', 'C']
    if (solar_elevation_deg < 20.0_dp) base = ['C', 'A', 'B']
    preferred = upper_branch(preferred_branch)
    if (index('ABC', preferred) == 0 .or. preferred == ' ') then
      order = base
      return
    end if
    order = [preferred]
    do i = 1, size(base)
      if (base(i) /= preferred) order = [order, base(i)]
    end do
  end subroutine outer_branch_order

  pure character(len=1) function upper_branch(branch) result(upper)
    character(len=1), intent(in) :: branch

    upper = branch
    if (branch >= 'a' .and. branch <= 'z') upper = achar(iachar(branch) - 32)
  end function upper_branch

  !> sheath-modelで指定順にbranchを試し、プロファイル全域の成立条件を満たす最初の零電流根を返す。
  !!
  !! ionは従来のZhao closureと同じ冷たいビーム、電子は設定の内向きdriftを持つ上流Maxwell分布とする。
  !! 内向きdriftのある電子で反射集団を持つA/Cは無限遠へ接続できず、sheath-modelが棄却する。
  subroutine solve_outer_root(app, source_density_m3, source_temperature_ev, order, root, success, message, previous)
    type(app_config), intent(in) :: app
    real(dp), intent(in) :: source_density_m3, source_temperature_ev
    character(len=1), intent(in) :: order(:)
    type(sheath_equilibrium_result), intent(out) :: root
    logical, intent(out) :: success
    character(len=*), intent(out) :: message
    type(surface_current_model_result_type), intent(in), optional :: previous
    type(fixed_entry_equilibrium_input) :: input
    type(sheath_solver) :: solver
    type(sheath_equilibrium_result) :: guesses(1)
    character(len=256) :: attempt_message
    integer(i32) :: electron_idx, ion_idx, status
    integer :: i
    logical :: use_guess

    success = .false.
    message = ''
    electron_idx = species_index(app, app%surface_current%electron_species)
    ion_idx = species_index(app, app%surface_current%ion_species)
    input%plasma%ion_density_m3 = species_number_density_m3(app%particle_species(ion_idx))
    input%plasma%electron_temperature_ev = &
      species_temperature_k(app%particle_species(electron_idx))*k_boltzmann/qe
    input%plasma%ion_temperature_ev = 0.0_dp
    input%plasma%ion_pressure_factor = 1.0_dp
    input%plasma%electron_drift_mps = -app%particle_species(electron_idx)%drift_velocity(3)
    input%plasma%ion_entry_speed_mps = -app%particle_species(ion_idx)%drift_velocity(3)
    input%plasma%ion_mass_kg = app%particle_species(ion_idx)%m_particle
    input%plasma%electron_mass_kg = app%particle_species(electron_idx)%m_particle
    input%plasma%photoelectrons = maxwellian_photoelectrons(source_density_m3, source_temperature_ev)
    use_guess = .false.
    if (present(previous)) use_guess = previous%active .and. index('ABC', previous%zhao_branch) > 0
    if (use_guess) then
      guesses(1)%valid = .true.
      guesses(1)%branch = previous%zhao_branch
      guesses(1)%surface_potential_v = previous%phi0_v
      guesses(1)%minimum_potential_v = previous%phi_m_v
      guesses(1)%ambient_electron_density_m3 = previous%ambient_electron_density_m3
    end if

    do i = 1, size(order)
      input%branch = order(i)
      if (use_guess) then
        call solver%solve_equilibrium(input, root, status, attempt_message, initial_guesses=guesses)
      else
        call solver%solve_equilibrium(input, root, status, attempt_message)
      end if
      if (status == SHEATH_OK .and. root%valid) then
        success = .true.
        message = ''
        return
      end if
      if (len_trim(message) > 0) message = trim(message)//' | '
      message = trim(message)//order(i)//': '//trim(attempt_message)
    end do
  end subroutine solve_outer_root

  !> Maxwell規格化密度と束の換算に使う、温度 [eV] の電子の熱速度 sqrt(2T/m) [m/s]。
  real(dp) function thermal_speed(app, temperature_ev) result(speed)
    type(app_config), intent(in) :: app
    real(dp), intent(in) :: temperature_ev

    speed = sqrt(2.0_dp*qe*temperature_ev/ &
                 app%particle_species(species_index(app, app%surface_current%electron_species))%m_particle)
  end function thermal_speed

  !> sheath-modelの零電流根から、表面放出を別に与えて species 別の固定電流targetと境界写像を作る。
  !!
  !! source_density_m3/source_temperature_ev は外部シースが見る放出源。表面放出との差は
  !! 周期セル内で再吸収されたPEであり、return target に含まれる。
  subroutine finish_zhao_closure(app, root, source_density_m3, source_temperature_ev, result, success, message)
    type(app_config), intent(in) :: app
    type(sheath_equilibrium_result), intent(in) :: root
    real(dp), intent(in) :: source_density_m3, source_temperature_ev
    type(surface_current_model_result_type), intent(inout) :: result
    logical, intent(out) :: success
    character(len=*), intent(out) :: message
    integer(i32) :: electron_idx, ion_idx, photo_idx
    real(dp) :: area, budget_scale, budget_tolerance, electron_bottleneck_potential_v
    real(dp) :: surface_density_m3, surface_temperature_ev, emission_current_density
    logical :: photoelectron_active

    success = .false.
    message = ''
    electron_idx = species_index(app, app%surface_current%electron_species)
    ion_idx = species_index(app, app%surface_current%ion_species)
    photoelectron_active = app%surface_current%photoelectron_source_scale > 0.0_dp
    photo_idx = 0_i32
    if (photoelectron_active) photo_idx = species_index(app, app%surface_current%photoelectron_species)
    area = (app%sim%box_max(1) - app%sim%box_min(1))*(app%sim%box_max(2) - app%sim%box_min(2))
    if (app%surface_current%has_reference_area_m2) area = app%surface_current%reference_area_m2
    if (.not. ieee_is_finite(area) .or. area <= 0.0_dp) then
      message = 'surface-current reference area must be finite and positive.'
      return
    end if
    if (index('ABC', root%branch) == 0 .or. root%branch == ' ') then
      message = 'surface-current root returned an unknown branch.'
      return
    end if

    result%zhao_branch = root%branch
    result%phi0_v = root%surface_potential_v
    result%phi_m_v = root%surface_potential_v
    if (root%branch == 'A') result%phi_m_v = root%minimum_potential_v
    result%ambient_electron_density_m3 = root%ambient_electron_density_m3
    call surface_photoelectron_source(app, surface_density_m3, surface_temperature_ev)
    emission_current_density = qe*surface_density_m3*thermal_speed(app, surface_temperature_ev)/(2.0_dp*sqrt(pi))
    result%electron_current_density_a_m2 = -qe*root%electron_inward_flux_m2_s
    result%ion_current_density_a_m2 = qe*root%ion_inward_flux_m2_s
    result%photoelectron_emission_current_density_a_m2 = emission_current_density
    result%photoelectron_escape_current_density_a_m2 = 0.0_dp
    if (photoelectron_active) result%photoelectron_escape_current_density_a_m2 = qe*root%photoelectron_escape_flux_m2_s
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
        result%photoelectron_return_current_density_a_m2 > 0.0_dp .or. &
        result%photoelectron_escape_current_density_a_m2 < 0.0_dp) then
      message = 'surface-current evaluation produced invalid channel signs.'
      return
    end if
    if (photoelectron_active) then
      if (result%photoelectron_emission_current_density_a_m2 <= 0.0_dp) then
        message = 'photoelectron closure requires a positive emission current.'
        return
      end if
    else if (any([ &
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
    result%outer_photoelectron_flux_m2_s = source_density_m3*thermal_speed(app, source_temperature_ev)/(2.0_dp*sqrt(pi))
    result%outer_photoelectron_mean_energy_ev = source_temperature_ev
    if (.not. photoelectron_active) then
      result%outer_photoelectron_flux_m2_s = 0.0_dp
      result%outer_photoelectron_mean_energy_ev = 0.0_dp
    end if
    result%has_absorbed_target([electron_idx, ion_idx]) = .true.
    result%absorbed_current_a(electron_idx) = &
      checked_area_current(area, result%electron_current_density_a_m2)
    result%absorbed_current_a(ion_idx) = checked_area_current(area, result%ion_current_density_a_m2)
    if (photoelectron_active) then
      result%has_absorbed_target(photo_idx) = .true.
      result%has_emission_target(photo_idx) = .true.
      result%has_escape_target(photo_idx) = .true.
      result%absorbed_current_a(photo_idx) = &
        checked_area_current(area, result%photoelectron_return_current_density_a_m2)
      result%emission_current_a(photo_idx) = &
        checked_area_current(area, result%photoelectron_emission_current_density_a_m2)
      ! escaped_to_infinity は粒子電荷の外向きfluxなので、正の表面帯電電流とは符号が逆。
      result%escaped_particle_current_a(photo_idx) = &
        -checked_area_current(area, result%photoelectron_escape_current_density_a_m2)
    end if

    ! Zhao の1-D外部シースを、z-high interfaceに対するkinetic boundary mapへ縮約する。
    ! Type Aの電子はphi_mがaccess bottleneckであり、Type B/Cはphi_infinity=0を使う。
    electron_bottleneck_potential_v = 0.0_dp
    if (result%zhao_branch == 'A') electron_bottleneck_potential_v = result%phi_m_v
    result%kinetic_contract = 'zhao_barrier_v1'
    result%has_inflow_kinetic_map([electron_idx, ion_idx]) = .true.
    result%inflow_reservoir_potential_v([electron_idx, ion_idx]) = 0.0_dp
    ! 上流電子のMaxwellian密度は、壁で吸われて戻らない高速電子を見込んで無限遠の準中性から決まる。
    ! 設定の太陽風密度をそのまま電子源に使うと流入fluxが根の電子電流からずれる。
    result%inflow_reservoir_density_m3(electron_idx) = result%ambient_electron_density_m3
    result%inflow_reservoir_density_m3(ion_idx) = species_number_density_m3(app%particle_species(ion_idx))
    result%inflow_access_potential_v(electron_idx) = electron_bottleneck_potential_v
    result%inflow_access_potential_v(ion_idx) = 0.0_dp
    result%inflow_kinetic_face([electron_idx, ion_idx]) = 6_i32
    result%has_outflow_kinetic_barrier([electron_idx, ion_idx]) = .true.
    result%outflow_barrier_potential_v(electron_idx) = electron_bottleneck_potential_v
    result%outflow_barrier_potential_v(ion_idx) = 0.0_dp
    result%outflow_barrier_face([electron_idx, ion_idx]) = 6_i32
    if (photoelectron_active) then
      result%has_outflow_kinetic_barrier(photo_idx) = .true.
      result%outflow_barrier_potential_v(photo_idx) = electron_bottleneck_potential_v
      result%outflow_barrier_face(photo_idx) = 6_i32
    end if
    ! 外部障壁と流入写像は上流0 Vからの電位差なので、z-high面の平均電位を壁電位phi0へ固定する。
    ! 外部での横移動はセル幅よりはるかに大きいので、x/y周期セルでは戻り位置を面内一様にする。
    result%has_plane_gauge = .true.
    result%plane_gauge_potential_v = result%phi0_v
    result%outer_return_cell_uniform = all(app%sim%bc_low(1:2) == bc_periodic) .and. &
                                       all(app%sim%bc_high(1:2) == bc_periodic)
    success = .true.
  end subroutine finish_zhao_closure

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
