!> zhao_stationary の外部シース源を、matching plane で観測した PE 流出から定期的に解き直す弱連成。
!!
!! 周期セルは Debye 長よりはるかに小さいため、外部の 1-D シースから見ると一様な壁であり、
!! その壁が外部へ出す PE は H を外向きに通過した流出だけである。表面放出は設定値のまま保ち、
!! `outflow_refresh_batches` 個の accepted batch で数えた H 外向きの法線運動エネルギーのスペクトルを
!! 外部 Zhao の bin 放出源へ与えて零電流根を解き直す。セル内の加速や電位極小で H のスペクトルは
!! 表面放出の半 Maxwell 分布から変形するため、平均だけの半 Maxwell 分布へは縮約しない。
!! 総電荷は毎 batch の固定電流 target が浮遊条件へ拘束する。
module bem_zhao_outflow_refresh
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use bem_irradiation_cell_contract, only: &
    irradiation_cell_contract_type, resolve_irradiation_cell_contract, irradiation_cell_ok
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe
  use bem_types, only: mesh_type, sim_stats
  use bem_app_config, only: app_config
  use bem_string_utils, only: lower_ascii
  use bem_charge_ledger, only: finite_charge_sum
  use bem_surface_closure_contract, only: surface_closure_contract_type
  use bem_surface_current_model, only: surface_current_model_result_type, evaluate_surface_current_model, &
                                       solve_zhao_outflow_closure, split_emission_escape_target
  use bem_mpi, only: mpi_context, mpi_is_root, mpi_allreduce_sum_real_dp_array
  use bem_outflow_spectrum, only: outflow_spectrum_bins
  implicit none
  private

  type, public :: zhao_outflow_refresh_type
    private
    logical :: active = .false.
    integer(i32) :: refresh_batches = 0_i32
    integer(i32) :: window_batches = 0_i32
    integer(i32) :: photo_idx = 0_i32
    integer(i32) :: secondary_idx = 0_i32
    real(dp) :: window_duration = 0.0_dp
    integer(i32) :: accepted_solves = 0_i32
    real(dp) :: area = 0.0_dp
    real(dp) :: photo_abs_charge = 0.0_dp
    real(dp) :: last_change = 0.0_dp
    !> 窓内の H 外向き通過数、その法線運動エネルギー和 [J]、表面放出数。raw weight 単位。
    real(dp) :: window_moments(3, 2) = 0.0_dp
    !> 窓内の H 外向き通過数の法線運動エネルギー bin 別内訳。raw weight 単位。
    real(dp) :: window_spectrum(outflow_spectrum_bins, 2) = 0.0_dp
    type(irradiation_cell_contract_type) :: cell
    type(surface_current_model_result_type) :: state
  contains
    procedure :: initialize
    procedure :: is_active
    procedure :: commit_batch
  end type zhao_outflow_refresh_type

contains

  !> 新規 run では設定の表面放出を外部源とする静的根から始め、再開では保存した外部源で根を再構成する。
  subroutine initialize(self, app, mesh, stats, surface_closure, cell)
    class(zhao_outflow_refresh_type), intent(out) :: self
    type(app_config), intent(in) :: app
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(inout) :: stats
    type(surface_closure_contract_type), intent(inout) :: surface_closure
    type(irradiation_cell_contract_type), intent(in), optional :: cell
    integer(i32) :: status
    character(len=256) :: message
    logical :: success

    self%active = trim(lower_ascii(app%surface_current%model)) == 'zhao_stationary' .and. &
                  app%surface_current%outflow_refresh_batches > 0_i32
    if (.not. self%active) then
      ! 外部根を更新しない run は、checkpoint に残る外部状態を引き継がない。
      stats%matching_plane_state_valid = .false.
      return
    end if
    if (present(cell)) then
      self%cell = cell
    else
      call resolve_irradiation_cell_contract(app, self%cell, status, message)
      if (status /= irradiation_cell_ok) error stop trim(message)
    end if
    self%refresh_batches = app%surface_current%outflow_refresh_batches
    self%area = self%cell%wall_area(app%sim%box_min, app%sim%box_max)
    if (app%surface_current%has_reference_area_m2) self%area = app%surface_current%reference_area_m2

    if (stats%matching_plane_state_valid) then
      if (stats%matching_plane_feedback(1) == 0.0_dp) then
        ! Before the first SEE-only refresh, the saved source is zero and the root is Type C.
        call evaluate_surface_current_model(app, self%state, self%cell)
        success = .true.
      else if (any(stats%matching_plane_photoelectron_spectrum_m2_s > 0.0_dp)) then
        call solve_zhao_outflow_closure( &
          app, stats%matching_plane_feedback(1), stats%matching_plane_feedback(2), &
          branch_from_response(stats%matching_plane_response), self%state, success, message, &
          outer_spectrum_m2_s=stats%matching_plane_photoelectron_spectrum_m2_s, &
          secondary_share=stats%matching_plane_secondary_share, cell=self%cell &
          )
      else
        ! 外部根をまだ更新していない checkpoint の外部源は、表面放出の半 Maxwell 分布である。
        call solve_zhao_outflow_closure( &
          app, stats%matching_plane_feedback(1), stats%matching_plane_feedback(2), &
          branch_from_response(stats%matching_plane_response), self%state, success, message, &
          secondary_share=stats%matching_plane_secondary_share, cell=self%cell &
          )
      end if
      if (.not. success) then
        error stop 'zhao_stationary outflow refresh could not reconstruct the checkpoint outer state: '// &
          trim(message)
      end if
      self%accepted_solves = stats%matching_plane_iterations
      self%last_change = stats%matching_plane_residual
    else
      call evaluate_surface_current_model(app, self%state, self%cell)
      self%accepted_solves = 1_i32
      self%last_change = 0.0_dp
    end if
    self%photo_idx = self%state%photoelectron_species_idx
    self%secondary_idx = self%state%secondary_electron_species_idx
    if (self%photo_idx > 0_i32) self%photo_abs_charge = abs(app%particle_species(self%photo_idx)%q_particle)
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
    real(dp) :: batch_moments(3, 2), transmission, outer_flux, outer_energy_ev, displacement
    real(dp) :: photo_flux, secondary_flux, secondary_share, tracked_secondary_current, expected_secondary_current
    real(dp) :: packed_values(2*(2 + outflow_spectrum_bins)), denominator
    integer :: group, idx
    real(dp) :: batch_values(2 + outflow_spectrum_bins, 2), outer_spectrum(outflow_spectrum_bins)
    character(len=256) :: message
    logical :: success

    refreshed = .false.
    if (.not. self%active) return
    batch_values = 0.0_dp
    batch_moments = 0.0_dp
    do group = 1, 2
      idx = self%photo_idx
      if (group == 2) idx = self%secondary_idx
      if (idx == 0) cycle
      batch_values(1:2, group) = sum(moments_thread(1:2, idx, :), dim=2)
      batch_values(3:, group) = sum(moments_thread(5:, idx, :), dim=2)
    end do
    packed_values = reshape(batch_values, [size(packed_values)])
    call mpi_allreduce_sum_real_dp_array(mpi_ctx, packed_values)
    batch_values = reshape(packed_values, shape(batch_values))
    batch_moments(1:2, :) = batch_values(1:2, :)
    if (self%photo_idx > 0) then
      batch_moments(3, 1) = fixed_current_charge_values(app%n_particle_species + self%photo_idx)/self%photo_abs_charge
    end if
    self%window_moments = self%window_moments + batch_moments
    self%window_spectrum = self%window_spectrum + batch_values(3:, :)
    self%window_batches = self%window_batches + 1_i32
    self%window_duration = self%window_duration + app%sim%batch_duration
    if (self%window_batches >= self%refresh_batches) then
      photo_flux = 0.0_dp
      secondary_flux = self%window_moments(1, 2)/(self%area*self%window_duration)
      transmission = 0.0_dp
      if (self%window_moments(3, 1) > 0.0_dp) then
        transmission = self%window_moments(1, 1)/self%window_moments(3, 1)
        photo_flux = transmission*self%state%photoelectron_emission_current_density_a_m2/qe
      end if
      outer_flux = photo_flux + secondary_flux
      secondary_share = 0.0_dp
      if (outer_flux > 0.0_dp) then
        outer_spectrum = 0.0_dp
        outer_energy_ev = 0.0_dp
        if (photo_flux > 0.0_dp) then
          outer_spectrum = photo_flux*self%window_spectrum(:, 1)/sum(self%window_spectrum(:, 1))
          outer_energy_ev = photo_flux*self%window_moments(2, 1)/(self%window_moments(1, 1)*qe)
        end if
        if (secondary_flux > 0.0_dp) then
          outer_spectrum = outer_spectrum + secondary_flux*self%window_spectrum(:, 2)/sum(self%window_spectrum(:, 2))
          outer_energy_ev = outer_energy_ev + secondary_flux*self%window_moments(2, 2)/(self%window_moments(1, 2)*qe)
        end if
        outer_energy_ev = outer_energy_ev/outer_flux
        secondary_share = secondary_flux/outer_flux
        call solve_zhao_outflow_closure( &
          app, outer_flux, outer_energy_ev, self%state%zhao_branch, trial, success, message, previous=self%state, &
          outer_spectrum_m2_s=outer_spectrum, secondary_share=secondary_share, cell=self%cell &
          )
        if (success) then
          denominator = self%state%outer_photoelectron_flux_m2_s
          if (denominator <= 0.0_dp) denominator = outer_flux
          self%last_change = abs(outer_flux - self%state%outer_photoelectron_flux_m2_s)/denominator
          denominator = self%state%outer_photoelectron_mean_energy_ev
          if (denominator <= 0.0_dp) denominator = outer_energy_ev
          self%last_change = max(self%last_change, &
                                 abs(outer_energy_ev - self%state%outer_photoelectron_mean_energy_ev)/denominator)
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
          'WARNING: zhao outflow refresh observed no surface-electron outflow; kept the previous outer state: batch=', &
          batch_idx
        flush (error_unit)
      end if
      if (.not. refreshed) then
        call split_emission_escape_target(self%state, secondary_share)
        surface_closure = self%state%surface_closure_contract_type
      end if
      self%window_moments = 0.0_dp
      self%window_spectrum = 0.0_dp
      self%window_batches = 0_i32
      self%window_duration = 0.0_dp
    end if

    ! Compare this batch's tracked escape with the root and split published in the same receipt.
    tracked_secondary_current = 0.0_dp
    expected_secondary_current = 0.0_dp
    if (self%secondary_idx > 0) then
      tracked_secondary_current = -fixed_current_charge_values(2*app%n_particle_species + self%secondary_idx)/ &
                                  (self%area*app%sim%batch_duration)
      expected_secondary_current = self%state%photoelectron_escape_current_density_a_m2*self%state%secondary_share
    end if
    stats_candidate%secondary_budget_residual_a_m2 = tracked_secondary_current - expected_secondary_current
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
    stats%matching_plane_photoelectron_spectrum_m2_s = self%state%outer_photoelectron_spectrum_m2_s
    stats%matching_plane_photoelectron_return_flux_m2_s = return_flux
    stats%matching_plane_photoelectron_escape_flux_m2_s = escape_flux
    stats%matching_plane_iterations = self%accepted_solves
    stats%matching_plane_residual = self%last_change
    stats%matching_plane_secondary_share = self%state%secondary_share
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
