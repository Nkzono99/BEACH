!> zhao_stationary の外部シース源を、matching plane で観測した PE 流出から定期的に解き直す弱連成。
!!
!! 周期セルは Debye 長よりはるかに小さいため、外部の 1-D シースから見ると一様な壁であり、
!! その壁が外部へ出す PE は H を外向きに通過した流出だけである。表面放出は設定値のまま保ち、
!! `outflow_refresh_batches` 個の accepted batch で平均した透過率と平均法線エネルギーを
!! 外部 Zhao の放出源へ与えて零電流根を解き直す。総電荷は毎 batch の固定電流 target が浮遊条件へ拘束する。
module bem_zhao_outflow_refresh
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe
  use bem_types, only: mesh_type, sim_stats
  use bem_app_config, only: app_config
  use bem_string_utils, only: lower_ascii
  use bem_charge_ledger, only: finite_charge_sum
  use bem_surface_closure_contract, only: surface_closure_contract_type
  use bem_surface_current_model, only: surface_current_model_result_type, evaluate_surface_current_model, &
                                       solve_zhao_outflow_closure
  use bem_mpi, only: mpi_context, mpi_is_root, mpi_allreduce_sum_real_dp_array
  implicit none
  private

  type, public :: zhao_outflow_refresh_type
    private
    logical :: active = .false.
    integer(i32) :: refresh_batches = 0_i32
    integer(i32) :: window_batches = 0_i32
    integer(i32) :: photo_idx = 0_i32
    integer(i32) :: accepted_solves = 0_i32
    real(dp) :: area = 0.0_dp
    real(dp) :: photo_abs_charge = 0.0_dp
    real(dp) :: last_change = 0.0_dp
    !> 窓内の H 外向き通過数、その法線運動エネルギー和 [J]、表面放出数。raw weight 単位。
    real(dp) :: window_moments(3) = 0.0_dp
    type(surface_current_model_result_type) :: state
  contains
    procedure :: initialize
    procedure :: is_active
    procedure :: commit_batch
  end type zhao_outflow_refresh_type

contains

  !> 新規 run では設定の表面放出を外部源とする静的根から始め、再開では保存した外部源で根を再構成する。
  subroutine initialize(self, app, mesh, stats, surface_closure)
    class(zhao_outflow_refresh_type), intent(out) :: self
    type(app_config), intent(in) :: app
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(inout) :: stats
    type(surface_closure_contract_type), intent(inout) :: surface_closure
    character(len=256) :: message
    logical :: success

    self%active = trim(lower_ascii(app%surface_current%model)) == 'zhao_stationary' .and. &
                  app%surface_current%outflow_refresh_batches > 0_i32
    if (.not. self%active) return
    self%refresh_batches = app%surface_current%outflow_refresh_batches
    self%area = product(app%sim%box_max(1:2) - app%sim%box_min(1:2))
    if (app%surface_current%has_reference_area_m2) self%area = app%surface_current%reference_area_m2

    if (stats%matching_plane_state_valid) then
      call solve_zhao_outflow_closure( &
        app, stats%matching_plane_feedback(1), stats%matching_plane_feedback(2), &
        branch_from_response(stats%matching_plane_response), self%state, success, message &
        )
      if (.not. success) then
        error stop 'zhao_stationary outflow refresh could not reconstruct the checkpoint outer state: '// &
          trim(message)
      end if
      self%accepted_solves = stats%matching_plane_iterations
      self%last_change = stats%matching_plane_residual
    else
      call evaluate_surface_current_model(app, self%state)
      self%accepted_solves = 1_i32
      self%last_change = 0.0_dp
    end if
    self%photo_idx = self%state%photoelectron_species_idx
    self%photo_abs_charge = abs(app%particle_species(self%photo_idx)%q_particle)
    surface_closure = self%state%surface_closure_contract_type
    call stage_state( &
      self, stats, finite_charge_sum(mesh%q_elem, 'zhao outflow refresh initial charge')/self%area &
      )
  end subroutine initialize

  logical function is_active(self) result(active)
    class(zhao_outflow_refresh_type), intent(in) :: self
    active = self%active
  end function is_active

  !> accepted batch の H 通過モーメントを窓へ加え、窓が満ちたら外部根を解き直して次 batch の closure を返す。
  subroutine commit_batch( &
    self, app, mesh, mpi_ctx, moments_thread, fixed_current_charge_values, batch_idx, stats_candidate, &
    surface_closure, refreshed &
    )
    class(zhao_outflow_refresh_type), intent(inout) :: self
    type(app_config), intent(in) :: app
    type(mesh_type), intent(in) :: mesh
    type(mpi_context), intent(in) :: mpi_ctx
    real(dp), intent(in) :: moments_thread(:, :, :)
    !> 固定電流closureがrank合計したraw電荷。n+species番目が表面放出。
    real(dp), intent(in) :: fixed_current_charge_values(:)
    integer(i32), intent(in) :: batch_idx
    type(sim_stats), intent(inout) :: stats_candidate
    type(surface_closure_contract_type), intent(inout) :: surface_closure
    logical, intent(out) :: refreshed
    type(surface_current_model_result_type) :: trial
    real(dp) :: batch_moments(3), transmission, outer_flux, outer_energy_ev, displacement
    character(len=256) :: message
    logical :: success

    refreshed = .false.
    if (.not. self%active) return
    batch_moments(1:2) = sum(moments_thread(1:2, self%photo_idx, :), dim=2)
    batch_moments(3) = 0.0_dp
    call mpi_allreduce_sum_real_dp_array(mpi_ctx, batch_moments)
    batch_moments(3) = fixed_current_charge_values(app%n_particle_species + self%photo_idx)/self%photo_abs_charge
    if (.not. all(ieee_is_finite(batch_moments)) .or. any(batch_moments < 0.0_dp)) then
      error stop 'zhao_stationary outflow refresh received invalid matching-plane moments.'
    end if
    self%window_moments = self%window_moments + batch_moments
    self%window_batches = self%window_batches + 1_i32

    if (self%window_batches >= self%refresh_batches) then
      if (self%window_moments(1) > 0.0_dp .and. self%window_moments(3) > 0.0_dp) then
        transmission = self%window_moments(1)/self%window_moments(3)
        outer_flux = transmission*self%state%photoelectron_emission_current_density_a_m2/qe
        outer_energy_ev = self%window_moments(2)/(self%window_moments(1)*qe)
        call solve_zhao_outflow_closure( &
          app, outer_flux, outer_energy_ev, self%state%zhao_branch, trial, success, message &
          )
        if (success) then
          self%last_change = max( &
                             abs(outer_flux/self%state%outer_photoelectron_flux_m2_s - 1.0_dp), &
                             abs(outer_energy_ev/self%state%outer_photoelectron_mean_energy_ev - 1.0_dp) &
                             )
          self%state = trial
          self%accepted_solves = self%accepted_solves + 1_i32
          surface_closure = self%state%surface_closure_contract_type
          refreshed = .true.
          if (mpi_is_root(mpi_ctx)) then
            write (output_unit, '(a,i0,a,a1,4(a,es13.5))') &
              'zhao outflow refresh: batch=', batch_idx, ' branch=', self%state%zhao_branch, &
              ' transmission=', transmission, ' mean_normal_energy_eV=', outer_energy_ev, &
              ' phi0_V=', self%state%phi0_v, ' relative_change=', self%last_change
            flush (output_unit)
          end if
        else if (mpi_is_root(mpi_ctx)) then
          write (error_unit, '(a,i0,a,es24.16,a,es24.16,a)') &
            'WARNING: zhao outflow refresh kept the previous outer state: batch=', batch_idx, &
            ', outer_flux_m2_s=', outer_flux, ', mean_normal_energy_eV=', outer_energy_ev, &
            ': '//trim(message)
          flush (error_unit)
        end if
      else if (mpi_is_root(mpi_ctx)) then
        write (error_unit, '(a,i0)') &
          'WARNING: zhao outflow refresh observed no photoelectron outflow; kept the previous outer state: batch=', &
          batch_idx
        flush (error_unit)
      end if
      self%window_moments = 0.0_dp
      self%window_batches = 0_i32
    end if

    displacement = finite_charge_sum(mesh%q_elem, 'zhao outflow refresh committed charge')/self%area
    call stage_state(self, stats_candidate, displacement)
  end subroutine commit_batch

  !> 現在の外部状態を matching-plane state として保存し、履歴・summary・checkpoint へ渡す。
  subroutine stage_state(self, stats, displacement)
    type(zhao_outflow_refresh_type), intent(in) :: self
    type(sim_stats), intent(inout) :: stats
    real(dp), intent(in) :: displacement
    real(dp) :: escape_flux, return_flux, access_potential

    access_potential = 0.0_dp
    if (self%state%zhao_branch == 'A') access_potential = self%state%phi_m_v
    escape_flux = self%state%photoelectron_escape_current_density_a_m2/qe
    return_flux = max(0.0_dp, self%state%outer_photoelectron_flux_m2_s - escape_flux)
    escape_flux = self%state%outer_photoelectron_flux_m2_s - return_flux
    stats%matching_plane_state_valid = .true.
    stats%matching_plane_displacement_c_m2 = displacement
    stats%matching_plane_phi_v = self%state%phi0_v
    stats%matching_plane_response = [ &
                                    self%state%phi0_v, -self%state%electron_current_density_a_m2/qe, &
                                    self%state%ion_current_density_a_m2/qe, access_potential, 0.0_dp, &
                                    access_potential &
                                    ]
    stats%matching_plane_feedback = [ &
                                    self%state%outer_photoelectron_flux_m2_s, &
                                    self%state%outer_photoelectron_mean_energy_ev, 0.0_dp, 0.0_dp &
                                    ]
    stats%matching_plane_photoelectron_return_flux_m2_s = return_flux
    stats%matching_plane_photoelectron_escape_flux_m2_s = escape_flux
    stats%matching_plane_iterations = self%accepted_solves
    stats%matching_plane_residual = self%last_change
  end subroutine stage_state

  !> 保存した外部障壁と壁電位の符号から Zhao branch を復元する。
  pure character(len=1) function branch_from_response(response) result(branch)
    real(dp), intent(in) :: response(6)

    if (response(6) < 0.0_dp) then
      branch = 'A'
    else if (response(1) >= 0.0_dp) then
      branch = 'B'
    else
      branch = 'C'
    end if
  end function branch_from_response

end module bem_zhao_outflow_refresh
