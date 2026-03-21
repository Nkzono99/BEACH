!> Zhao 系シース条件から注入フラックス補正を解決するモジュール。
module bem_sheath_injection_model
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_types, only: sim_config
  use bem_app_config_types, only: app_config, particle_species_spec
  use bem_app_config_parser, only: &
    lower, resolve_inject_face, resolve_inward_normal, species_number_density_m3, species_temperature_k
  use bem_injection, only: compute_inflow_flux_from_drifting_maxwellian
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  real(dp), parameter :: pi = 3.1415926535897932384626433832795d0
  real(dp), parameter :: qe = 1.602176634d-19
  real(dp), parameter :: nonlinear_tol = 1.0d-5
  integer, parameter :: nonlinear_max_iter = 60
  integer, parameter :: nonlinear_max_backtrack = 20

  type :: sheath_injection_context
    logical :: enabled = .false.
    logical :: has_photo_species = .false.
    character(len=32) :: model = 'none'
    character(len=1) :: branch = ' '
    integer(i32) :: electron_species = 0_i32
    integer(i32) :: ion_species = 0_i32
    integer(i32) :: photo_species = 0_i32
    character(len=16) :: reference_face = ''
    integer(i32) :: reference_axis = 0_i32
    real(dp) :: reference_coordinate = 0.0d0
    real(dp) :: reference_inward_normal(3) = 0.0d0
    real(dp) :: phi0_v = 0.0d0
    real(dp) :: phi_m_v = 0.0d0
    real(dp) :: n_swe_inf_m3 = 0.0d0
    real(dp) :: electron_number_density_m3 = 0.0d0
    real(dp) :: electron_vmin_normal = 0.0d0
    real(dp) :: photo_emit_current_density_a_m2 = 0.0d0
    real(dp) :: photo_vmin_normal = 0.0d0
  end type sheath_injection_context

  public :: sheath_injection_context
  public :: resolve_sheath_injection_context

  abstract interface
    subroutine nonlinear_residual(x, f)
      import :: dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f(:)
    end subroutine nonlinear_residual
  end interface

  type :: zhao_params_type
    real(dp) :: alpha_rad = 0.0d0
    real(dp) :: n_swi_inf_m3 = 0.0d0
    real(dp) :: n_phe_ref_m3 = 0.0d0
    real(dp) :: n_phe0_m3 = 0.0d0
    real(dp) :: t_swe_ev = 0.0d0
    real(dp) :: t_phe_ev = 0.0d0
    real(dp) :: v_d_electron_mps = 0.0d0
    real(dp) :: v_d_ion_mps = 0.0d0
    real(dp) :: m_i_kg = 0.0d0
    real(dp) :: m_e_kg = 0.0d0
    real(dp) :: v_swe_th_mps = 0.0d0
    real(dp) :: v_phe_th_mps = 0.0d0
    real(dp) :: cs_mps = 0.0d0
    real(dp) :: mach = 0.0d0
    real(dp) :: u = 0.0d0
    real(dp) :: tau = 0.0d0
  end type zhao_params_type

contains

  !> 設定から Zhao/無光電子シース注入コンテキストを解決する。
  subroutine resolve_sheath_injection_context(cfg, ctx)
    type(app_config), intent(in) :: cfg
    type(sheath_injection_context), intent(out) :: ctx

    integer(i32) :: electron_idx, ion_idx, photo_idx
    type(particle_species_spec) :: spec_e, spec_i, spec_p
    real(dp) :: phi0_v, phi_m_v, n_swe_inf_m3
    real(dp) :: n_e_inf_m3, n_i_inf_m3, n_phe_ref_m3, n_phe0_m3
    real(dp) :: t_swe_ev, t_phe_ev, delta_phi_swe, delta_phi_phe
    real(dp) :: v_d_electron_mps, v_d_ion_mps
    integer :: reference_axis
    real(dp) :: reference_coordinate, reference_inward_normal(3)
    type(zhao_params_type) :: params
    character(len=32) :: model
    character(len=1) :: branch

    ctx = sheath_injection_context()
    model = trim(lower(cfg%sim%sheath_injection_model))
    ctx%model = model
    if (model == 'none') return

    call detect_sheath_species(cfg, electron_idx, ion_idx, photo_idx)
    spec_e = cfg%particle_species(electron_idx)
    spec_i = cfg%particle_species(ion_idx)

    ctx%enabled = .true.
    ctx%electron_species = electron_idx
    ctx%ion_species = ion_idx
    ctx%reference_face = trim(lower(spec_e%inject_face))
    call resolve_sheath_reference_plane( &
      cfg%sim, ctx%reference_face, reference_axis, reference_coordinate, reference_inward_normal &
      )
    ctx%reference_axis = int(reference_axis, i32)
    ctx%reference_coordinate = reference_coordinate
    ctx%reference_inward_normal = reference_inward_normal

    n_e_inf_m3 = species_number_density_m3(spec_e)
    n_i_inf_m3 = species_number_density_m3(spec_i)
    t_swe_ev = temperature_ev_from_spec(spec_e)
    v_d_electron_mps = resolve_species_drift_speed( &
                       spec_e, cfg%sim%sheath_electron_drift_mode, reference_inward_normal &
                       )
    v_d_ion_mps = resolve_species_drift_speed( &
                  spec_i, cfg%sim%sheath_ion_drift_mode, reference_inward_normal &
                  )

    if (v_d_electron_mps <= 0.0d0) error stop 'sheath electron drift must point inward.'
    if (v_d_ion_mps <= 0.0d0) error stop 'sheath ion drift must point inward.'

    select case (model)
    case ('floating_no_photo')
      call solve_no_photo_floating_potential(spec_e, spec_i, reference_inward_normal, phi0_v)
      ctx%branch = 'N'
      ctx%phi0_v = phi0_v
      ctx%phi_m_v = phi0_v
      ctx%n_swe_inf_m3 = n_e_inf_m3
      ctx%electron_number_density_m3 = n_e_inf_m3
      ctx%electron_vmin_normal = sqrt(max(0.0d0, 2.0d0*qe*max(0.0d0, -phi0_v))/abs(spec_e%m_particle))
      return
    case ('zhao_auto', 'zhao_a', 'zhao_b', 'zhao_c')
      if (photo_idx <= 0_i32) then
        error stop 'Zhao sheath injection requires one enabled negative-q photo_raycast species.'
      end if
      spec_p = cfg%particle_species(photo_idx)
      if (abs(spec_p%q_particle) <= 0.0d0) error stop 'Zhao photoelectron species must have non-zero q_particle.'
      ctx%photo_species = photo_idx
      ctx%has_photo_species = .true.
      t_phe_ev = temperature_ev_from_spec(spec_p)
      n_phe_ref_m3 = cfg%sim%sheath_photoelectron_ref_density_cm3*1.0d6
      call build_zhao_params( &
        cfg%sim%sheath_alpha_deg, n_i_inf_m3, n_phe_ref_m3, t_swe_ev, t_phe_ev, v_d_electron_mps, v_d_ion_mps, &
        spec_i%m_particle, abs(spec_e%m_particle), params &
        )
      call solve_zhao_unknowns(model, params, phi0_v, phi_m_v, n_swe_inf_m3, branch)

      ctx%branch = branch
      ctx%phi0_v = phi0_v
      ctx%phi_m_v = phi_m_v
      ctx%n_swe_inf_m3 = n_swe_inf_m3
      ctx%electron_number_density_m3 = n_swe_inf_m3

      select case (branch)
      case ('A')
        delta_phi_swe = -phi_m_v
        delta_phi_phe = phi0_v - phi_m_v
      case ('B')
        delta_phi_swe = 0.0d0
        delta_phi_phe = phi0_v
      case ('C')
        delta_phi_swe = -phi0_v
        delta_phi_phe = 0.0d0
      case default
        error stop 'Unexpected Zhao sheath branch.'
      end select

      ctx%electron_vmin_normal = sqrt(max(0.0d0, 2.0d0*qe*max(0.0d0, delta_phi_swe))/abs(spec_e%m_particle))
      ctx%photo_vmin_normal = sqrt(max(0.0d0, 2.0d0*qe*max(0.0d0, delta_phi_phe))/abs(spec_p%m_particle))

      n_phe0_m3 = params%n_phe0_m3
      select case (branch)
      case ('A')
        ctx%photo_emit_current_density_a_m2 = abs(spec_p%q_particle)* &
                                              n_phe0_m3*exp((phi_m_v - phi0_v)/t_phe_ev)* &
                                              params%v_phe_th_mps/(2.0d0*sqrt(pi))
      case ('B')
        ctx%photo_emit_current_density_a_m2 = abs(spec_p%q_particle)* &
                                              n_phe0_m3*exp(-phi0_v/t_phe_ev)* &
                                              params%v_phe_th_mps/(2.0d0*sqrt(pi))
      case ('C')
        ctx%photo_emit_current_density_a_m2 = abs(spec_p%q_particle)* &
                                              n_phe0_m3*params%v_phe_th_mps/(2.0d0*sqrt(pi))
      end select
    case default
      error stop 'Unknown sim.sheath_injection_model in runtime.'
    end select
  end subroutine resolve_sheath_injection_context

  subroutine detect_sheath_species(cfg, electron_idx, ion_idx, photo_idx)
    type(app_config), intent(in) :: cfg
    integer(i32), intent(out) :: electron_idx, ion_idx, photo_idx

    integer(i32) :: s
    character(len=16) :: mode

    electron_idx = 0_i32
    ion_idx = 0_i32
    photo_idx = 0_i32
    do s = 1_i32, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      mode = trim(lower(cfg%particle_species(s)%source_mode))
      select case (mode)
      case ('reservoir_face')
        if (cfg%particle_species(s)%q_particle < 0.0d0) then
          if (electron_idx == 0_i32) electron_idx = s
        else if (cfg%particle_species(s)%q_particle > 0.0d0) then
          if (ion_idx == 0_i32) ion_idx = s
        end if
      case ('photo_raycast')
        if (cfg%particle_species(s)%q_particle < 0.0d0 .and. photo_idx == 0_i32) photo_idx = s
      end select
    end do

    if (electron_idx <= 0_i32) error stop 'sheath injection requires one enabled negative-q reservoir_face species.'
    if (ion_idx <= 0_i32) error stop 'sheath injection requires one enabled positive-q reservoir_face species.'
    if (trim(lower(cfg%particle_species(electron_idx)%inject_face)) /= trim(lower(cfg%particle_species(ion_idx)%inject_face))) then
      error stop 'sheath electron/ion reservoir species must share the same inject_face.'
    end if
  end subroutine detect_sheath_species

  subroutine solve_no_photo_floating_potential(spec_e, spec_i, inward_normal, phi0_v)
    type(particle_species_spec), intent(in) :: spec_e, spec_i
    real(dp), intent(in) :: inward_normal(3)
    real(dp), intent(out) :: phi0_v

    real(dp) :: n_e_inf_m3, n_i_inf_m3, gamma_i, f_low, f_high, f_mid, phi_low, phi_high, phi_mid
    integer :: iter

    n_e_inf_m3 = species_number_density_m3(spec_e)
    n_i_inf_m3 = species_number_density_m3(spec_i)
    gamma_i = compute_inflow_flux_from_drifting_maxwellian( &
              n_i_inf_m3, species_temperature_k(spec_i), spec_i%m_particle, spec_i%drift_velocity, inward_normal &
              )
    if (gamma_i <= 0.0d0) error stop 'floating_no_photo requires a positive ion inflow flux.'

    phi_low = -128.0d0*max(1.0d0, temperature_ev_from_spec(spec_e))
    phi_high = 0.0d0
    f_low = no_photo_current_balance(phi_low, n_e_inf_m3, spec_e, inward_normal, gamma_i)
    f_high = no_photo_current_balance(phi_high, n_e_inf_m3, spec_e, inward_normal, gamma_i)
    if (f_high < 0.0d0) then
      error stop 'floating_no_photo could not bracket a negative sheath potential.'
    end if
    do while (f_low > 0.0d0)
      phi_low = 2.0d0*phi_low
      f_low = no_photo_current_balance(phi_low, n_e_inf_m3, spec_e, inward_normal, gamma_i)
      if (phi_low < -1.0d6) error stop 'floating_no_photo bracket search failed.'
    end do

    do iter = 1, 80
      phi_mid = 0.5d0*(phi_low + phi_high)
      f_mid = no_photo_current_balance(phi_mid, n_e_inf_m3, spec_e, inward_normal, gamma_i)
      if (abs(f_mid) <= 1.0d-12*max(1.0d0, gamma_i)) exit
      if (f_mid > 0.0d0) then
        phi_high = phi_mid
      else
        phi_low = phi_mid
      end if
    end do
    phi0_v = 0.5d0*(phi_low + phi_high)
  end subroutine solve_no_photo_floating_potential

  real(dp) function no_photo_current_balance(phi0_v, n_e_inf_m3, spec_e, inward_normal, gamma_i) result(balance)
    real(dp), intent(in) :: phi0_v
    real(dp), intent(in) :: n_e_inf_m3
    type(particle_species_spec), intent(in) :: spec_e
    real(dp), intent(in) :: inward_normal(3)
    real(dp), intent(in) :: gamma_i
    real(dp) :: vmin_normal

    vmin_normal = sqrt(max(0.0d0, 2.0d0*qe*max(0.0d0, -phi0_v))/abs(spec_e%m_particle))
    balance = compute_inflow_flux_from_drifting_maxwellian( &
              n_e_inf_m3, species_temperature_k(spec_e), spec_e%m_particle, spec_e%drift_velocity, inward_normal, &
              vmin_normal=vmin_normal &
              ) - gamma_i
  end function no_photo_current_balance

  subroutine build_zhao_params( &
    alpha_deg, n_swi_inf_m3, n_phe_ref_m3, t_swe_ev, t_phe_ev, v_d_electron_mps, v_d_ion_mps, m_i_kg, m_e_kg, p &
    )
    real(dp), intent(in) :: alpha_deg, n_swi_inf_m3, n_phe_ref_m3, t_swe_ev, t_phe_ev
    real(dp), intent(in) :: v_d_electron_mps, v_d_ion_mps, m_i_kg, m_e_kg
    type(zhao_params_type), intent(out) :: p

    if (t_swe_ev <= 0.0d0) error stop 'Zhao sheath requires electron temperature > 0.'
    if (t_phe_ev <= 0.0d0) error stop 'Zhao sheath requires photoelectron temperature > 0.'
    if (n_swi_inf_m3 <= 0.0d0) error stop 'Zhao sheath requires ion density > 0.'
    if (n_phe_ref_m3 <= 0.0d0) error stop 'Zhao sheath requires sheath_photoelectron_ref_density_cm3 > 0.'
    if (v_d_ion_mps <= 0.0d0) error stop 'Zhao sheath requires positive ion drift.'
    if (m_i_kg <= 0.0d0 .or. m_e_kg <= 0.0d0) error stop 'Zhao sheath requires positive particle masses.'

    p%alpha_rad = alpha_deg*pi/180.0d0
    p%n_swi_inf_m3 = n_swi_inf_m3
    p%n_phe_ref_m3 = n_phe_ref_m3
    p%n_phe0_m3 = n_phe_ref_m3*sin(p%alpha_rad)
    p%t_swe_ev = t_swe_ev
    p%t_phe_ev = t_phe_ev
    p%v_d_electron_mps = v_d_electron_mps
    p%v_d_ion_mps = v_d_ion_mps
    p%m_i_kg = m_i_kg
    p%m_e_kg = m_e_kg
    p%v_swe_th_mps = sqrt(2.0d0*qe*p%t_swe_ev/p%m_e_kg)
    p%v_phe_th_mps = sqrt(2.0d0*qe*p%t_phe_ev/p%m_e_kg)
    p%cs_mps = sqrt(qe*p%t_swe_ev/p%m_i_kg)
    p%mach = p%v_d_ion_mps/p%cs_mps
    p%u = p%v_d_electron_mps/p%v_swe_th_mps
    p%tau = p%t_swe_ev/p%t_phe_ev

    if (.not. ieee_is_finite(p%mach) .or. p%mach <= 0.0d0) then
      error stop 'Zhao sheath produced an invalid Mach number.'
    end if
  end subroutine build_zhao_params

  subroutine solve_zhao_unknowns(model, p, phi0_v, phi_m_v, n_swe_inf_m3, branch)
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(out) :: phi0_v, phi_m_v, n_swe_inf_m3
    character(len=1), intent(out) :: branch

    real(dp) :: x3(3), x2(2)
    logical :: success
    character(len=1), dimension(3) :: order
    integer :: i

    select case (trim(model))
    case ('zhao_a')
      call solve_zhao_branch_a(p, phi0_v, phi_m_v, n_swe_inf_m3)
      branch = 'A'
      return
    case ('zhao_b')
      call solve_zhao_branch_b(p, phi0_v, n_swe_inf_m3)
      phi_m_v = phi0_v
      branch = 'B'
      return
    case ('zhao_c')
      call solve_zhao_branch_c(p, phi0_v, n_swe_inf_m3)
      phi_m_v = phi0_v
      branch = 'C'
      return
    case ('zhao_auto')
      if (p%alpha_rad*180.0d0/pi < 20.0d0) then
        order = ['C', 'A', 'B']
      else
        order = ['A', 'B', 'C']
      end if
    case default
      error stop 'Unknown Zhao sheath model.'
    end select

    do i = 1, size(order)
      select case (order(i))
      case ('A')
        call try_solve_zhao_branch_a(p, x3, success)
        if (success) then
          phi0_v = x3(1)
          phi_m_v = x3(2)
          n_swe_inf_m3 = x3(3)
          branch = 'A'
          return
        end if
      case ('B')
        call try_solve_zhao_branch_b(p, x2, success)
        if (success) then
          phi0_v = x2(1)
          phi_m_v = x2(1)
          n_swe_inf_m3 = x2(2)
          branch = 'B'
          return
        end if
      case ('C')
        call try_solve_zhao_branch_c(p, x2, success)
        if (success) then
          phi0_v = x2(1)
          phi_m_v = x2(1)
          n_swe_inf_m3 = x2(2)
          branch = 'C'
          return
        end if
      end select
    end do

    error stop 'Zhao sheath auto branch selection failed.'
  end subroutine solve_zhao_unknowns

  subroutine solve_zhao_branch_a(p, phi0_v, phi_m_v, n_swe_inf_m3)
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(out) :: phi0_v, phi_m_v, n_swe_inf_m3
    real(dp) :: x(3)
    logical :: success

    call try_solve_zhao_branch_a(p, x, success)
    if (.not. success) error stop 'Zhao Type-A root solve failed.'
    phi0_v = x(1)
    phi_m_v = x(2)
    n_swe_inf_m3 = x(3)
  end subroutine solve_zhao_branch_a

  subroutine solve_zhao_branch_b(p, phi0_v, n_swe_inf_m3)
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(out) :: phi0_v, n_swe_inf_m3
    real(dp) :: x(2)
    logical :: success

    call try_solve_zhao_branch_b(p, x, success)
    if (.not. success) error stop 'Zhao Type-B root solve failed.'
    phi0_v = x(1)
    n_swe_inf_m3 = x(2)
  end subroutine solve_zhao_branch_b

  subroutine solve_zhao_branch_c(p, phi0_v, n_swe_inf_m3)
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(out) :: phi0_v, n_swe_inf_m3
    real(dp) :: x(2)
    logical :: success

    call try_solve_zhao_branch_c(p, x, success)
    if (.not. success) error stop 'Zhao Type-C root solve failed.'
    phi0_v = x(1)
    n_swe_inf_m3 = x(2)
  end subroutine solve_zhao_branch_c

  subroutine try_solve_zhao_branch_a(p, x, success)
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(out) :: x(3)
    logical, intent(out) :: success

    real(dp) :: guesses(3, 3)

    guesses(:, 1) = [3.6d0, -0.5d0, 8.2d6]
    guesses(:, 2) = [2.8d0, -0.3d0, 8.0d6]
    guesses(:, 3) = [4.5d0, -0.8d0, 8.4d6]
    call solve_nonlinear_system(3, guesses, residual_a, x, success)

  contains

    subroutine residual_a(xa, fa)
      real(dp), intent(in) :: xa(:)
      real(dp), intent(out) :: fa(:)

      call zhao_residuals_type_a(p, xa, fa)
    end subroutine residual_a

  end subroutine try_solve_zhao_branch_a

  subroutine try_solve_zhao_branch_b(p, x, success)
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(out) :: x(2)
    logical, intent(out) :: success

    real(dp) :: guesses(2, 3)

    guesses(:, 1) = [1.3d0, 7.0d6]
    guesses(:, 2) = [0.8d0, 6.5d6]
    guesses(:, 3) = [2.0d0, 7.8d6]
    call solve_nonlinear_system(2, guesses, residual_b, x, success)

  contains

    subroutine residual_b(xb, fb)
      real(dp), intent(in) :: xb(:)
      real(dp), intent(out) :: fb(:)

      call zhao_residuals_type_b(p, xb, fb)
    end subroutine residual_b

  end subroutine try_solve_zhao_branch_b

  subroutine try_solve_zhao_branch_c(p, x, success)
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(out) :: x(2)
    logical, intent(out) :: success

    real(dp) :: guesses(2, 5)

    guesses(:, 1) = [-0.5d0, 6.0d6]
    guesses(:, 2) = [-2.0d0, 7.0d6]
    guesses(:, 3) = [-5.0d0, 8.0d6]
    guesses(:, 4) = [-10.0d0, 8.2d6]
    guesses(:, 5) = [-15.0d0, 8.5d6]
    call solve_nonlinear_system(2, guesses, residual_c, x, success)

  contains

    subroutine residual_c(xc, fc)
      real(dp), intent(in) :: xc(:)
      real(dp), intent(out) :: fc(:)

      call zhao_residuals_type_c(p, xc, fc)
    end subroutine residual_c

  end subroutine try_solve_zhao_branch_c

  subroutine zhao_residuals_type_a(p, x, f)
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)

    real(dp) :: phi0_v, phi_m_v, n_swe_inf_m3, a_swe, a_phe, ion_term

    phi0_v = x(1)
    phi_m_v = x(2)
    n_swe_inf_m3 = x(3)
    if (phi0_v <= 0.0d0 .or. phi_m_v >= 0.0d0 .or. phi_m_v >= phi0_v .or. n_swe_inf_m3 <= 0.0d0) then
      f = 1.0d6
      return
    end if

    a_swe = sqrt(max(0.0d0, -phi_m_v/p%t_swe_ev)) - p%u
    a_phe = sqrt(max(0.0d0, -phi_m_v/p%t_phe_ev))
    ion_term = p%n_swi_inf_m3*sqrt(2.0d0*pi*p%t_swe_ev/p%t_phe_ev*p%m_e_kg/p%m_i_kg)*p%mach

    f(1) = 0.5d0*n_swe_inf_m3*(1.0d0 + 2.0d0*erf(p%u) + erf(a_swe)) + &
           0.5d0*p%n_phe0_m3*exp(-phi0_v/p%t_phe_ev)*(1.0d0 - erf(a_phe)) - p%n_swi_inf_m3
    f(2) = p%n_phe0_m3*exp((phi_m_v - phi0_v)/p%t_phe_ev) - swe_free_current_term(p, n_swe_inf_m3, a_swe) + ion_term
    f(3) = type_a_e2_sum_at_infinity(p, phi0_v, phi_m_v, n_swe_inf_m3)
  end subroutine zhao_residuals_type_a

  subroutine zhao_residuals_type_b(p, x, f)
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)

    real(dp) :: phi0_v, n_swe_inf_m3, ion_term

    phi0_v = x(1)
    n_swe_inf_m3 = x(2)
    if (phi0_v <= 0.0d0 .or. n_swe_inf_m3 <= 0.0d0) then
      f = 1.0d6
      return
    end if

    ion_term = p%n_swi_inf_m3*sqrt(2.0d0*pi*p%t_swe_ev/p%t_phe_ev*p%m_e_kg/p%m_i_kg)*p%mach
    f(1) = 0.5d0*n_swe_inf_m3*(1.0d0 + erf(p%u)) + 0.5d0*p%n_phe0_m3*exp(-phi0_v/p%t_phe_ev) - p%n_swi_inf_m3
    f(2) = p%n_phe0_m3*exp(-phi0_v/p%t_phe_ev) - swe_free_current_term(p, n_swe_inf_m3, -p%u) + ion_term
  end subroutine zhao_residuals_type_b

  subroutine zhao_residuals_type_c(p, x, f)
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)

    real(dp) :: phi0_v, n_swe_inf_m3, a_swe, a_phe, ion_term

    phi0_v = x(1)
    n_swe_inf_m3 = x(2)
    if (phi0_v >= 0.0d0 .or. n_swe_inf_m3 <= 0.0d0) then
      f = 1.0d6
      return
    end if

    a_swe = sqrt(max(0.0d0, -phi0_v/p%t_swe_ev)) - p%u
    a_phe = sqrt(max(0.0d0, -phi0_v/p%t_phe_ev))
    ion_term = p%n_swi_inf_m3*sqrt(2.0d0*pi*p%t_swe_ev/p%t_phe_ev*p%m_e_kg/p%m_i_kg)*p%mach

    f(1) = 0.5d0*n_swe_inf_m3*(1.0d0 + 2.0d0*erf(p%u) + erf(a_swe)) + &
           0.5d0*p%n_phe0_m3*exp(-phi0_v/p%t_phe_ev)*erfc(a_phe) - p%n_swi_inf_m3
    f(2) = p%n_phe0_m3 - swe_free_current_term(p, n_swe_inf_m3, a_swe) + ion_term
  end subroutine zhao_residuals_type_c

  real(dp) function swe_free_current_term(p, n_swe_inf_m3, a_swe) result(term)
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(in) :: n_swe_inf_m3, a_swe

    term = n_swe_inf_m3*(sqrt(p%t_swe_ev/p%t_phe_ev)*exp(-(a_swe*a_swe)) + &
                         sqrt(pi)*(p%v_d_electron_mps/p%v_phe_th_mps)*erfc(a_swe))
  end function swe_free_current_term

  real(dp) function type_a_e2_sum_at_infinity(p, phi0_v, phi_m_v, n_swe_inf_m3) result(e2_sum)
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(in) :: phi0_v, phi_m_v, n_swe_inf_m3

    real(dp) :: phi, s_swe, s_phe, e2_swe_f, e2_swe_r, e2_phe_f, e2_swi, arg_phi, arg_m

    if (abs(p%u) <= 1.0d-12) then
      e2_sum = 1.0d30
      return
    end if

    phi = 0.0d0
    s_swe = sqrt(max(0.0d0, (phi - phi_m_v)/p%t_swe_ev))
    s_phe = sqrt(max(0.0d0, (phi - phi_m_v)/p%t_phe_ev))

    e2_swe_f = (p%t_swe_ev/p%t_phe_ev)*(n_swe_inf_m3/p%n_phe_ref_m3)*( &
               exp(phi/p%t_swe_ev)*(1.0d0 - erf(s_swe - p%u)) - &
               exp(phi_m_v/p%t_swe_ev)*(1.0d0 - erf(-p%u)) + &
               (1.0d0/(sqrt(pi)*p%u))*exp(phi_m_v/p%t_swe_ev - p%u*p%u)*(exp(2.0d0*p%u*s_swe) - 1.0d0) &
               )

    e2_swe_r = 2.0d0*(p%t_swe_ev/p%t_phe_ev)*(n_swe_inf_m3/p%n_phe_ref_m3)*( &
               exp(phi/p%t_swe_ev)*(erf(s_swe - p%u) + erf(p%u)) - &
               (1.0d0/(sqrt(pi)*p%u))*exp(phi_m_v/p%t_swe_ev - p%u*p%u)*(exp(2.0d0*p%u*s_swe) - 1.0d0) &
               )

    e2_phe_f = (p%n_phe0_m3/p%n_phe_ref_m3)*( &
               exp((phi - phi0_v)/p%t_phe_ev)*(1.0d0 - erf(s_phe)) - &
               exp((phi_m_v - phi0_v)/p%t_phe_ev)*(1.0d0 - 2.0d0*s_phe/sqrt(pi)) &
               )

    arg_phi = 1.0d0 - 2.0d0*phi/(p%t_swe_ev*p%mach*p%mach)
    arg_m = 1.0d0 - 2.0d0*phi_m_v/(p%t_swe_ev*p%mach*p%mach)
    if (arg_phi <= 0.0d0 .or. arg_m <= 0.0d0) then
      e2_sum = 1.0d30
      return
    end if

    e2_swi = 2.0d0*(p%t_swe_ev/p%t_phe_ev)*(p%n_swi_inf_m3/p%n_phe_ref_m3)*p%mach*p%mach*( &
             sqrt(arg_phi) - sqrt(arg_m) &
             )
    e2_sum = e2_swe_f + e2_swe_r + e2_phe_f + e2_swi
  end function type_a_e2_sum_at_infinity

  subroutine solve_nonlinear_system(n, guesses, residual_fn, x_best, success)
    integer, intent(in) :: n
    real(dp), intent(in) :: guesses(:, :)
    procedure(nonlinear_residual) :: residual_fn
    real(dp), intent(out) :: x_best(n)
    logical, intent(out) :: success

    integer :: guess_idx
    real(dp) :: x_trial(n), best_norm, trial_norm
    logical :: trial_success

    success = .false.
    if (size(guesses, 1) /= n) error stop 'solve_nonlinear_system guess dimension mismatch.'
    x_best = guesses(:, 1)
    best_norm = huge(1.0d0)
    do guess_idx = 1, size(guesses, 2)
      call try_newton_solve(n, guesses(:, guess_idx), residual_fn, x_trial, trial_norm, trial_success)
      if (trial_norm < best_norm) then
        best_norm = trial_norm
        x_best = x_trial
      end if
      if (trial_success .and. trial_norm < nonlinear_tol) then
        success = .true.
        x_best = x_trial
        return
      end if
    end do

    success = best_norm < nonlinear_tol
  end subroutine solve_nonlinear_system

  subroutine try_newton_solve(n, x0, residual_fn, x_out, final_norm, success)
    integer, intent(in) :: n
    real(dp), intent(in) :: x0(n)
    procedure(nonlinear_residual) :: residual_fn
    real(dp), intent(out) :: x_out(n)
    real(dp), intent(out) :: final_norm
    logical, intent(out) :: success

    integer :: iter, backtrack
    real(dp) :: x(n), f(n), jac(n, n), dx(n), x_trial(n), f_trial(n), step_scale, fnorm, trial_norm
    logical :: linear_ok, improved

    x = x0
    call residual_fn(x, f)
    fnorm = residual_norm(f)
    do iter = 1, nonlinear_max_iter
      if (fnorm < nonlinear_tol) exit
      call numerical_jacobian(n, x, f, residual_fn, jac)
      call solve_small_linear_system(n, jac, -f, dx, linear_ok)
      if (.not. linear_ok) exit

      step_scale = 1.0d0
      improved = .false.
      do backtrack = 1, nonlinear_max_backtrack
        x_trial = x + step_scale*dx
        call residual_fn(x_trial, f_trial)
        trial_norm = residual_norm(f_trial)
        if (trial_norm < fnorm) then
          x = x_trial
          f = f_trial
          fnorm = trial_norm
          improved = .true.
          exit
        end if
        step_scale = 0.5d0*step_scale
      end do
      if (.not. improved) exit
    end do

    x_out = x
    final_norm = fnorm
    success = fnorm < nonlinear_tol
  end subroutine try_newton_solve

  subroutine numerical_jacobian(n, x, f0, residual_fn, jac)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), f0(n)
    procedure(nonlinear_residual) :: residual_fn
    real(dp), intent(out) :: jac(n, n)

    integer :: j
    real(dp) :: h, xh(n), fh(n)

    do j = 1, n
      h = 1.0d-6*max(1.0d0, abs(x(j)))
      xh = x
      xh(j) = xh(j) + h
      call residual_fn(xh, fh)
      jac(:, j) = (fh - f0)/h
    end do
  end subroutine numerical_jacobian

  subroutine solve_small_linear_system(n, a_in, b_in, x, ok)
    integer, intent(in) :: n
    real(dp), intent(in) :: a_in(n, n), b_in(n)
    real(dp), intent(out) :: x(n)
    logical, intent(out) :: ok

    integer :: i, j, k, pivot_row
    real(dp) :: a(n, n), b(n), factor, pivot_abs, tmp_row(n), tmp_val

    a = a_in
    b = b_in
    ok = .true.

    do k = 1, n
      pivot_row = k
      pivot_abs = abs(a(k, k))
      do i = k + 1, n
        if (abs(a(i, k)) > pivot_abs) then
          pivot_abs = abs(a(i, k))
          pivot_row = i
        end if
      end do
      if (pivot_abs <= 1.0d-18) then
        ok = .false.
        x = 0.0d0
        return
      end if
      if (pivot_row /= k) then
        tmp_row = a(k, :)
        a(k, :) = a(pivot_row, :)
        a(pivot_row, :) = tmp_row
        tmp_val = b(k)
        b(k) = b(pivot_row)
        b(pivot_row) = tmp_val
      end if
      do i = k + 1, n
        factor = a(i, k)/a(k, k)
        a(i, k:n) = a(i, k:n) - factor*a(k, k:n)
        b(i) = b(i) - factor*b(k)
      end do
    end do

    x = 0.0d0
    do i = n, 1, -1
      x(i) = b(i)
      do j = i + 1, n
        x(i) = x(i) - a(i, j)*x(j)
      end do
      x(i) = x(i)/a(i, i)
    end do
  end subroutine solve_small_linear_system

  real(dp) function residual_norm(f) result(norm2)
    real(dp), intent(in) :: f(:)

    if (.not. all(ieee_is_finite(f))) then
      norm2 = huge(1.0d0)
      return
    end if
    norm2 = sqrt(sum(f*f))
  end function residual_norm

  real(dp) function temperature_ev_from_spec(spec) result(temp_ev)
    type(particle_species_spec), intent(in) :: spec

    temp_ev = species_temperature_k(spec)*k_boltzmann/qe
  end function temperature_ev_from_spec

  subroutine resolve_sheath_reference_plane(sim, inject_face, axis, reference_coordinate, inward_normal)
    type(sim_config), intent(in) :: sim
    character(len=*), intent(in) :: inject_face
    integer, intent(out) :: axis
    real(dp), intent(out) :: reference_coordinate
    real(dp), intent(out) :: inward_normal(3)
    real(dp) :: boundary_value

    call resolve_inward_normal(inject_face, inward_normal)
    call resolve_inject_face(sim%box_min, sim%box_max, inject_face, axis, boundary_value)
    if (sim%has_sheath_reference_coordinate) then
      reference_coordinate = sim%sheath_reference_coordinate
    else
      reference_coordinate = boundary_value
    end if
  end subroutine resolve_sheath_reference_plane

  real(dp) function resolve_species_drift_speed(spec, drift_mode, inward_normal) result(speed)
    type(particle_species_spec), intent(in) :: spec
    character(len=*), intent(in) :: drift_mode
    real(dp), intent(in) :: inward_normal(3)

    select case (trim(lower(drift_mode)))
    case ('full')
      speed = sqrt(sum(spec%drift_velocity*spec%drift_velocity))
    case ('normal')
      speed = dot_product(spec%drift_velocity, inward_normal)
    case default
      error stop 'Unknown sheath drift mode.'
    end select
  end function resolve_species_drift_speed

end module bem_sheath_injection_model
