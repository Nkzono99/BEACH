!> Zhao 系シース数値モデルの core 実装。
module bem_sheath_model_core
  use bem_kinds, only: dp
  use bem_constants, only: k_boltzmann
  use bem_injection, only: compute_inflow_flux_from_drifting_maxwellian
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  real(dp), parameter :: pi = 3.1415926535897932384626433832795d0
  real(dp), parameter :: eps0 = 8.8541878128d-12
  real(dp), parameter :: qe = 1.602176634d-19
  real(dp), parameter :: nonlinear_tol = 1.0d-5
  integer, parameter :: nonlinear_max_iter = 60
  integer, parameter :: nonlinear_max_backtrack = 20
  real(dp), parameter :: zhao_profile_phi_tol_hat = 1.0d-3
  real(dp), parameter :: zhao_type_a_phi_m_eps_hat = 1.0d-5
  integer, parameter :: zhao_profile_grid = 8000

  type :: sheath_model_species
    real(dp) :: q_particle = 0.0d0
    real(dp) :: m_particle = 0.0d0
    real(dp) :: number_density_m3 = 0.0d0
    real(dp) :: temperature_k = 0.0d0
    real(dp) :: drift_velocity(3) = 0.0d0
  end type sheath_model_species

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
    real(dp) :: lambda_d_phe_ref_m = 0.0d0
  end type zhao_params_type

  type :: zhao_local_state_type
    character(len=1) :: branch = ' '
    character(len=16) :: side = 'monotonic'
    real(dp) :: z_m_m = 0.0d0
    real(dp) :: phi_hat = 0.0d0
    real(dp) :: phi_v = 0.0d0
    real(dp) :: n_swi_m3 = 0.0d0
    real(dp) :: n_swe_f_m3 = 0.0d0
    real(dp) :: n_swe_r_m3 = 0.0d0
    real(dp) :: n_phe_f_m3 = 0.0d0
    real(dp) :: n_phe_c_m3 = 0.0d0
    real(dp) :: electron_source_density_m3 = 0.0d0
    real(dp) :: vcut_swe_mps = 0.0d0
    real(dp) :: vcut_phe_mps = 0.0d0
    real(dp) :: v_i_mps = 0.0d0
    logical :: swe_reflected_active = .false.
    logical :: phe_captured_active = .false.
  end type zhao_local_state_type

  public :: sheath_model_species
  public :: zhao_params_type
  public :: zhao_local_state_type
  public :: solve_no_photo_floating_potential
  public :: build_zhao_params
  public :: solve_zhao_unknowns
  public :: sample_zhao_state_at_z
  public :: resolve_species_drift_speed
  public :: zhao_electron_vmin_normal
  public :: zhao_photo_vmin_normal
  public :: zhao_photo_emit_current_density

contains

  subroutine solve_no_photo_floating_potential(spec_e, spec_i, inward_normal, phi0_v)
    type(sheath_model_species), intent(in) :: spec_e, spec_i
    real(dp), intent(in) :: inward_normal(3)
    real(dp), intent(out) :: phi0_v

    real(dp) :: n_e_inf_m3, n_i_inf_m3, gamma_i, f_low, f_high, f_mid, phi_low, phi_high, phi_mid
    integer :: iter

    n_e_inf_m3 = spec_e%number_density_m3
    n_i_inf_m3 = spec_i%number_density_m3
    gamma_i = compute_inflow_flux_from_drifting_maxwellian( &
              n_i_inf_m3, spec_i%temperature_k, spec_i%m_particle, spec_i%drift_velocity, inward_normal &
              )
    if (gamma_i <= 0.0d0) error stop 'floating_no_photo requires a positive ion inflow flux.'

    phi_low = -128.0d0*max(1.0d0, temperature_ev_from_species(spec_e))
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
    type(sheath_model_species), intent(in) :: spec_e
    real(dp), intent(in) :: inward_normal(3)
    real(dp), intent(in) :: gamma_i
    real(dp) :: vmin_normal

    vmin_normal = sqrt(max(0.0d0, 2.0d0*qe*max(0.0d0, -phi0_v))/abs(spec_e%m_particle))
    balance = compute_inflow_flux_from_drifting_maxwellian( &
              n_e_inf_m3, spec_e%temperature_k, spec_e%m_particle, spec_e%drift_velocity, inward_normal, &
              vmin_normal=vmin_normal &
              ) - gamma_i
  end function no_photo_current_balance

  subroutine sample_zhao_state_at_z(p, branch, phi0_v, phi_m_v, n_swe_inf_m3, z_m, state)
    type(zhao_params_type), intent(in) :: p
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: phi0_v, phi_m_v, n_swe_inf_m3, z_m
    type(zhao_local_state_type), intent(out) :: state

    real(dp) :: phi0_hat, phi_m_hat, phi_hat, z_m_hat
    character(len=16) :: side

    phi0_hat = phi0_v/p%t_phe_ev
    phi_m_hat = phi_m_v/p%t_phe_ev
    if (z_m <= 0.0d0) then
      phi_hat = phi0_hat
      if (branch == 'A') then
        side = 'lower'
      else
        side = 'monotonic'
      end if
      z_m_hat = 0.0d0
    else
      select case (branch)
      case ('A')
        call sample_type_a_phi_hat_at_z(p, phi0_hat, phi_m_hat, n_swe_inf_m3/p%n_phe_ref_m3, z_m, phi_hat, side, z_m_hat)
      case ('B', 'C')
        call sample_monotonic_phi_hat_at_z( &
          p, branch, phi0_hat, phi_m_hat, n_swe_inf_m3/p%n_phe_ref_m3, z_m, phi_hat, side, z_m_hat &
          )
      case default
        error stop 'Unknown Zhao branch in local state reconstruction.'
      end select
    end if

    call evaluate_zhao_state_from_phi_hat( &
      p, branch, side, phi_hat, phi0_hat, phi_m_hat, n_swe_inf_m3, z_m_hat*p%lambda_d_phe_ref_m, state &
      )
  end subroutine sample_zhao_state_at_z

  subroutine sample_type_a_phi_hat_at_z(p, phi0_hat, phi_m_hat, n_swe_inf_hat, z_m, phi_hat, side, z_m_hat)
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(in) :: phi0_hat, phi_m_hat, n_swe_inf_hat, z_m
    real(dp), intent(out) :: phi_hat, z_m_hat
    character(len=16), intent(out) :: side

    integer :: i, ngrid_lower, ngrid_upper
    real(dp) :: z_hat_target, phi_m_eps_hat, phi_end_hat, x
    real(dp), allocatable :: phi_lower_asc(:), s_lower_asc(:), phi_upper_asc(:), s_upper_asc(:)
    real(dp), allocatable :: z_lower_desc(:), phi_lower_desc(:), z_upper_asc(:)

    z_hat_target = z_m/p%lambda_d_phe_ref_m
    phi_m_eps_hat = min(zhao_type_a_phi_m_eps_hat, 5.0d-2*max(1.0d-8, abs(phi_m_hat)))
    phi_m_eps_hat = max(phi_m_eps_hat, 1.0d-8)
    phi_end_hat = -abs(zhao_profile_phi_tol_hat)
    if (phi_end_hat <= phi_m_hat) phi_end_hat = 0.5d0*phi_m_hat

    ngrid_upper = max(2000, zhao_profile_grid)
    ngrid_lower = max(ngrid_upper/2, 1500)
    allocate (phi_lower_asc(ngrid_lower), s_lower_asc(ngrid_lower))
    allocate (phi_upper_asc(ngrid_upper), s_upper_asc(ngrid_upper), z_upper_asc(ngrid_upper))
    allocate (z_lower_desc(ngrid_lower), phi_lower_desc(ngrid_lower))

    do i = 1, ngrid_lower
      x = real(i - 1, dp)/real(ngrid_lower - 1, dp)
      phi_lower_asc(i) = (phi_m_hat + phi_m_eps_hat) + (phi0_hat - (phi_m_hat + phi_m_eps_hat))*x*x
    end do
    call build_type_a_branch_from_minimum(p, phi_lower_asc, phi0_hat, n_swe_inf_hat, phi_m_hat, 'lower', s_lower_asc)
    z_m_hat = s_lower_asc(ngrid_lower)
    do i = 1, ngrid_lower
      phi_lower_desc(i) = phi_lower_asc(ngrid_lower - i + 1)
      z_lower_desc(i) = z_m_hat - s_lower_asc(ngrid_lower - i + 1)
    end do

    do i = 1, ngrid_upper
      x = real(i - 1, dp)/real(ngrid_upper - 1, dp)
      phi_upper_asc(i) = (phi_m_hat + phi_m_eps_hat) + (phi_end_hat - (phi_m_hat + phi_m_eps_hat))*(2.0d0*x - x*x)
    end do
    call build_type_a_branch_from_minimum(p, phi_upper_asc, phi0_hat, n_swe_inf_hat, phi_m_hat, 'upper', s_upper_asc)
    z_upper_asc = z_m_hat + s_upper_asc

    if (z_hat_target <= z_m_hat) then
      side = 'lower'
      phi_hat = interpolate_profile_value(z_lower_desc, phi_lower_desc, z_hat_target)
    else if (z_hat_target <= z_upper_asc(ngrid_upper)) then
      side = 'upper'
      phi_hat = interpolate_profile_value(z_upper_asc, phi_upper_asc, z_hat_target)
    else
      side = 'upper'
      phi_hat = 0.0d0
    end if
  end subroutine sample_type_a_phi_hat_at_z

  subroutine sample_monotonic_phi_hat_at_z(p, branch, phi0_hat, phi_m_hat, n_swe_inf_hat, z_m, phi_hat, side, z_m_hat)
    type(zhao_params_type), intent(in) :: p
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: phi0_hat, phi_m_hat, n_swe_inf_hat, z_m
    real(dp), intent(out) :: phi_hat, z_m_hat
    character(len=16), intent(out) :: side

    integer :: i, ngrid
    real(dp) :: z_hat_target, x, phi_end_hat, dphi, integral
    real(dp), allocatable :: phi_nodes(:), rho_nodes(:), e2_nodes(:), z_nodes(:), cumulative(:)

    z_hat_target = z_m/p%lambda_d_phe_ref_m
    z_m_hat = 0.0d0
    side = 'monotonic'
    ngrid = max(2000, zhao_profile_grid)
    allocate (phi_nodes(ngrid), rho_nodes(ngrid), e2_nodes(ngrid), z_nodes(ngrid))

    select case (branch)
    case ('B')
      phi_end_hat = abs(zhao_profile_phi_tol_hat)
      do i = 1, ngrid
        x = real(i - 1, dp)/real(ngrid - 1, dp)
        phi_nodes(i) = phi0_hat + (phi_end_hat - phi0_hat)*(2.0d0*x - x*x)
      end do
      do i = 1, ngrid
        call evaluate_zhao_rho_hat(p, branch, side, phi_nodes(i), phi0_hat, phi_m_hat, n_swe_inf_hat, rho_nodes(i))
      end do
      e2_nodes(ngrid) = 0.0d0
      integral = 0.0d0
      do i = ngrid - 1, 1, -1
        dphi = phi_nodes(i + 1) - phi_nodes(i)
        integral = integral + 0.5d0*(rho_nodes(i) + rho_nodes(i + 1))*dphi
        e2_nodes(i) = max(2.0d0*integral, 0.0d0)
      end do
    case ('C')
      phi_end_hat = -abs(zhao_profile_phi_tol_hat)
      allocate (cumulative(ngrid))
      do i = 1, ngrid
        x = real(i - 1, dp)/real(ngrid - 1, dp)
        phi_nodes(i) = phi0_hat + (phi_end_hat - phi0_hat)*(2.0d0*x - x*x)
      end do
      do i = 1, ngrid
        call evaluate_zhao_rho_hat(p, branch, side, phi_nodes(i), phi0_hat, phi_m_hat, n_swe_inf_hat, rho_nodes(i))
      end do
      cumulative(1) = 0.0d0
      do i = 2, ngrid
        dphi = phi_nodes(i) - phi_nodes(i - 1)
        cumulative(i) = cumulative(i - 1) + 0.5d0*(rho_nodes(i - 1) + rho_nodes(i))*dphi
      end do
      do i = 1, ngrid
        e2_nodes(i) = max(2.0d0*(cumulative(ngrid) - cumulative(i)), 0.0d0)
      end do
      deallocate (cumulative)
    case default
      error stop 'Unknown monotonic Zhao branch.'
    end select

    z_nodes(1) = 0.0d0
    do i = 2, ngrid
      dphi = abs(phi_nodes(i) - phi_nodes(i - 1))
      z_nodes(i) = z_nodes(i - 1) + dphi/sqrt(max(0.5d0*(e2_nodes(i - 1) + e2_nodes(i)), 1.0d-14))
    end do

    if (z_hat_target <= z_nodes(ngrid)) then
      phi_hat = interpolate_profile_value(z_nodes, phi_nodes, z_hat_target)
    else
      phi_hat = 0.0d0
    end if
  end subroutine sample_monotonic_phi_hat_at_z

  subroutine build_type_a_branch_from_minimum(p, phi_nodes_asc, phi0_hat, n_swe_inf_hat, phi_m_hat, side, s_nodes)
    type(zhao_params_type), intent(in) :: p
    real(dp), intent(in) :: phi_nodes_asc(:), phi0_hat, n_swe_inf_hat, phi_m_hat
    character(len=*), intent(in) :: side
    real(dp), intent(out) :: s_nodes(:)

    integer :: i
    real(dp) :: rho_m_hat, dphi0, integral_node, e2_mid, rho_m_neg, dphi
    real(dp), allocatable :: rho_nodes(:)

    if (size(phi_nodes_asc) /= size(s_nodes)) error stop 'Type-A profile work array size mismatch.'
    allocate (rho_nodes(size(phi_nodes_asc)))
    do i = 1, size(phi_nodes_asc)
      call evaluate_zhao_rho_hat(p, 'A', side, phi_nodes_asc(i), phi0_hat, phi_m_hat, n_swe_inf_hat, rho_nodes(i))
    end do
    call evaluate_zhao_rho_hat(p, 'A', side, phi_m_hat, phi0_hat, phi_m_hat, n_swe_inf_hat, rho_m_hat)

    rho_m_neg = max(-rho_m_hat, 1.0d-14)
    dphi0 = phi_nodes_asc(1) - phi_m_hat
    integral_node = 0.5d0*(rho_m_hat + rho_nodes(1))*dphi0
    s_nodes(1) = sqrt(max(0.0d0, 2.0d0*dphi0/rho_m_neg))

    do i = 2, size(phi_nodes_asc)
      dphi = phi_nodes_asc(i) - phi_nodes_asc(i - 1)
      integral_node = integral_node + 0.5d0*(rho_nodes(i - 1) + rho_nodes(i))*dphi
      e2_mid = max(-2.0d0*(integral_node - 0.25d0*(rho_nodes(i - 1) + rho_nodes(i))*dphi), 1.0d-14)
      s_nodes(i) = s_nodes(i - 1) + dphi/sqrt(e2_mid)
    end do
  end subroutine build_type_a_branch_from_minimum

  subroutine evaluate_zhao_state_from_phi_hat(p, branch, side, phi_hat, phi0_hat, phi_m_hat, n_swe_inf_m3, z_m_m, state)
    type(zhao_params_type), intent(in) :: p
    character(len=1), intent(in) :: branch
    character(len=*), intent(in) :: side
    real(dp), intent(in) :: phi_hat, phi0_hat, phi_m_hat, n_swe_inf_m3, z_m_m
    type(zhao_local_state_type), intent(out) :: state

    real(dp) :: n_swi_hat, n_swe_f_hat, n_swe_r_hat, n_phe_f_hat, n_phe_c_hat
    real(dp) :: arg_ion, a_swe, a_phe

    call evaluate_zhao_density_hat( &
      p, branch, side, phi_hat, phi0_hat, phi_m_hat, n_swe_inf_m3/p%n_phe_ref_m3, &
      n_swi_hat, n_swe_f_hat, n_swe_r_hat, n_phe_f_hat, n_phe_c_hat &
      )

    arg_ion = 1.0d0 - 2.0d0*phi_hat/(p%tau*p%mach*p%mach)
    if (arg_ion <= 0.0d0) error stop 'Zhao local ion energy argument became non-positive.'

    select case (branch)
    case ('A')
      a_swe = sqrt(max(0.0d0, (phi_hat - phi_m_hat)/p%tau))
      a_phe = sqrt(max(0.0d0, phi_hat - phi_m_hat))
      state%swe_reflected_active = trim(side) == 'upper'
      state%phe_captured_active = trim(side) == 'lower'
    case ('B')
      a_swe = 0.0d0
      a_phe = sqrt(max(0.0d0, phi_hat))
      state%swe_reflected_active = .false.
      state%phe_captured_active = .true.
    case ('C')
      a_swe = sqrt(max(0.0d0, (phi_hat - phi0_hat)/p%tau))
      a_phe = sqrt(max(0.0d0, phi_hat - phi0_hat))
      state%swe_reflected_active = .true.
      state%phe_captured_active = .false.
    case default
      error stop 'Unknown Zhao branch in local state evaluation.'
    end select

    state%branch = branch
    state%side = side
    state%z_m_m = z_m_m
    state%phi_hat = phi_hat
    state%phi_v = phi_hat*p%t_phe_ev
    state%n_swi_m3 = n_swi_hat*p%n_phe_ref_m3
    state%n_swe_f_m3 = n_swe_f_hat*p%n_phe_ref_m3
    state%n_swe_r_m3 = n_swe_r_hat*p%n_phe_ref_m3
    state%n_phe_f_m3 = n_phe_f_hat*p%n_phe_ref_m3
    state%n_phe_c_m3 = n_phe_c_hat*p%n_phe_ref_m3
    state%electron_source_density_m3 = n_swe_inf_m3*exp(phi_hat/p%tau)
    state%vcut_swe_mps = a_swe*p%v_swe_th_mps
    state%vcut_phe_mps = a_phe*p%v_phe_th_mps
    state%v_i_mps = p%v_d_ion_mps*sqrt(arg_ion)
  end subroutine evaluate_zhao_state_from_phi_hat

  subroutine evaluate_zhao_rho_hat(p, branch, side, phi_hat, phi0_hat, phi_m_hat, n_swe_inf_hat, rho_hat)
    type(zhao_params_type), intent(in) :: p
    character(len=1), intent(in) :: branch
    character(len=*), intent(in) :: side
    real(dp), intent(in) :: phi_hat, phi0_hat, phi_m_hat, n_swe_inf_hat
    real(dp), intent(out) :: rho_hat

    real(dp) :: n_swi_hat, n_swe_f_hat, n_swe_r_hat, n_phe_f_hat, n_phe_c_hat

    call evaluate_zhao_density_hat( &
      p, branch, side, phi_hat, phi0_hat, phi_m_hat, n_swe_inf_hat, &
      n_swi_hat, n_swe_f_hat, n_swe_r_hat, n_phe_f_hat, n_phe_c_hat &
      )
    rho_hat = n_swi_hat - n_swe_f_hat - n_swe_r_hat - n_phe_f_hat - n_phe_c_hat
  end subroutine evaluate_zhao_rho_hat

  subroutine evaluate_zhao_density_hat( &
    p, branch, side, phi_hat, phi0_hat, phi_m_hat, n_swe_inf_hat, &
    n_swi_hat, n_swe_f_hat, n_swe_r_hat, n_phe_f_hat, n_phe_c_hat &
    )
    type(zhao_params_type), intent(in) :: p
    character(len=1), intent(in) :: branch
    character(len=*), intent(in) :: side
    real(dp), intent(in) :: phi_hat, phi0_hat, phi_m_hat, n_swe_inf_hat
    real(dp), intent(out) :: n_swi_hat, n_swe_f_hat, n_swe_r_hat, n_phe_f_hat, n_phe_c_hat

    real(dp) :: arg_ion, s_swe, s_phe, sin_alpha

    sin_alpha = p%n_phe0_m3/p%n_phe_ref_m3
    arg_ion = 1.0d0 - 2.0d0*phi_hat/(p%tau*p%mach*p%mach)
    if (arg_ion <= 0.0d0) error stop 'Zhao ion density argument became non-positive.'
    n_swi_hat = (p%n_swi_inf_m3/p%n_phe_ref_m3)*arg_ion**(-0.5d0)

    select case (branch)
    case ('A')
      s_swe = sqrt(max(0.0d0, (phi_hat - phi_m_hat)/p%tau))
      s_phe = sqrt(max(0.0d0, phi_hat - phi_m_hat))
      n_swe_f_hat = 0.5d0*n_swe_inf_hat*exp(phi_hat/p%tau)*(1.0d0 - erf(s_swe - p%u))
      n_phe_f_hat = 0.5d0*sin_alpha*exp(phi_hat - phi0_hat)*(1.0d0 - erf(s_phe))
      if (trim(side) == 'lower') then
        n_swe_r_hat = 0.0d0
        n_phe_c_hat = sin_alpha*exp(phi_hat - phi0_hat)*erf(s_phe)
      else if (trim(side) == 'upper') then
        n_swe_r_hat = n_swe_inf_hat*exp(phi_hat/p%tau)*(erf(s_swe - p%u) + erf(p%u))
        n_phe_c_hat = 0.0d0
      else
        error stop 'Unknown Type-A Zhao side.'
      end if
    case ('B')
      s_phe = sqrt(max(0.0d0, phi_hat))
      n_swe_f_hat = 0.5d0*n_swe_inf_hat*exp(phi_hat/p%tau)*(1.0d0 + erf(p%u))
      n_swe_r_hat = 0.0d0
      n_phe_f_hat = 0.5d0*sin_alpha*exp(phi_hat - phi0_hat)*(1.0d0 - erf(s_phe))
      n_phe_c_hat = sin_alpha*exp(phi_hat - phi0_hat)*erf(s_phe)
    case ('C')
      s_swe = sqrt(max(0.0d0, (phi_hat - phi0_hat)/p%tau))
      s_phe = sqrt(max(0.0d0, phi_hat - phi0_hat))
      n_swe_f_hat = 0.5d0*n_swe_inf_hat*exp(phi_hat/p%tau)*(1.0d0 - erf(s_swe - p%u))
      n_swe_r_hat = n_swe_inf_hat*exp(phi_hat/p%tau)*(erf(s_swe - p%u) + erf(p%u))
      n_phe_f_hat = 0.5d0*sin_alpha*exp(phi_hat - phi0_hat)*erfc(s_phe)
      n_phe_c_hat = 0.0d0
    case default
      error stop 'Unknown Zhao branch in density evaluation.'
    end select
  end subroutine evaluate_zhao_density_hat

  real(dp) function interpolate_profile_value(x_nodes, y_nodes, x_query) result(y_query)
    real(dp), intent(in) :: x_nodes(:), y_nodes(:), x_query

    integer :: i
    real(dp) :: t

    if (size(x_nodes) /= size(y_nodes)) error stop 'Interpolation node size mismatch.'
    if (x_query <= x_nodes(1)) then
      y_query = y_nodes(1)
      return
    end if
    do i = 2, size(x_nodes)
      if (x_query <= x_nodes(i)) then
        t = (x_query - x_nodes(i - 1))/max(x_nodes(i) - x_nodes(i - 1), tiny(1.0d0))
        y_query = (1.0d0 - t)*y_nodes(i - 1) + t*y_nodes(i)
        return
      end if
    end do
    y_query = y_nodes(size(y_nodes))
  end function interpolate_profile_value

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
    p%lambda_d_phe_ref_m = sqrt(eps0*qe*p%t_phe_ev/(p%n_phe_ref_m3*qe*qe))

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

  real(dp) function temperature_ev_from_species(spec) result(temp_ev)
    type(sheath_model_species), intent(in) :: spec

    temp_ev = spec%temperature_k*k_boltzmann/qe
  end function temperature_ev_from_species

  real(dp) function resolve_species_drift_speed(spec, drift_mode, inward_normal) result(speed)
    type(sheath_model_species), intent(in) :: spec
    character(len=*), intent(in) :: drift_mode
    real(dp), intent(in) :: inward_normal(3)

    select case (trim(lower_ascii(drift_mode)))
    case ('full')
      speed = sqrt(sum(spec%drift_velocity*spec%drift_velocity))
    case ('normal')
      speed = dot_product(spec%drift_velocity, inward_normal)
    case default
      error stop 'Unknown sheath drift mode.'
    end select
  end function resolve_species_drift_speed

  real(dp) function zhao_electron_vmin_normal(branch, phi0_v, phi_m_v, m_e_kg) result(vmin_normal)
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: phi0_v, phi_m_v, m_e_kg
    real(dp) :: delta_phi_swe

    select case (branch)
    case ('A')
      delta_phi_swe = -phi_m_v
    case ('B')
      delta_phi_swe = 0.0d0
    case ('C')
      delta_phi_swe = -phi0_v
    case default
      error stop 'Unexpected Zhao sheath branch.'
    end select
    vmin_normal = sqrt(max(0.0d0, 2.0d0*qe*max(0.0d0, delta_phi_swe))/abs(m_e_kg))
  end function zhao_electron_vmin_normal

  real(dp) function zhao_photo_vmin_normal(branch, phi0_v, phi_m_v, m_phe_kg) result(vmin_normal)
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: phi0_v, phi_m_v, m_phe_kg
    real(dp) :: delta_phi_phe

    select case (branch)
    case ('A')
      delta_phi_phe = phi0_v - phi_m_v
    case ('B')
      delta_phi_phe = phi0_v
    case ('C')
      delta_phi_phe = 0.0d0
    case default
      error stop 'Unexpected Zhao sheath branch.'
    end select
    vmin_normal = sqrt(max(0.0d0, 2.0d0*qe*max(0.0d0, delta_phi_phe))/abs(m_phe_kg))
  end function zhao_photo_vmin_normal

  real(dp) function zhao_photo_emit_current_density(branch, phi0_v, phi_m_v, photo_charge_c, p) result(current_a_m2)
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: phi0_v, phi_m_v, photo_charge_c
    type(zhao_params_type), intent(in) :: p

    select case (branch)
    case ('A')
      current_a_m2 = abs(photo_charge_c)*p%n_phe0_m3*exp((phi_m_v - phi0_v)/p%t_phe_ev)*p%v_phe_th_mps/(2.0d0*sqrt(pi))
    case ('B')
      current_a_m2 = abs(photo_charge_c)*p%n_phe0_m3*exp(-phi0_v/p%t_phe_ev)*p%v_phe_th_mps/(2.0d0*sqrt(pi))
    case ('C')
      current_a_m2 = abs(photo_charge_c)*p%n_phe0_m3*p%v_phe_th_mps/(2.0d0*sqrt(pi))
    case default
      error stop 'Unexpected Zhao sheath branch.'
    end select
  end function zhao_photo_emit_current_density

  pure character(len=len(text)) function lower_ascii(text) result(lowered)
    character(len=*), intent(in) :: text

    integer :: i, code

    lowered = text
    do i = 1, len(text)
      code = iachar(lowered(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) lowered(i:i) = achar(code + 32)
    end do
  end function lower_ascii

end module bem_sheath_model_core
