!> Matching-plane coupling state and batch lifecycle; particle replay is owned by the simulator.
module bem_matching_plane_coupling
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe
  use bem_types, only: mesh_type, sim_stats
  use bem_app_config, only: app_config
  use bem_string_utils, only: lower_ascii
  use bem_charge_ledger, only: finite_charge_sum
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_surface_closure_contract, only: surface_closure_contract_type
  use bem_matching_plane_response_provider, only: matching_plane_response_provider_type, matching_plane_provider_ok
  use bem_matching_plane_zhao, only: matching_plane_zhao_root_seed_type
  use bem_matching_plane_implicit, only: solve_matching_implicit_zero_mode
  use bem_mpi, only: mpi_context, mpi_is_root, mpi_allreduce_sum_real_dp_array
  implicit none
  private

  type, public :: matching_plane_coupling_type
    private
    logical :: active = .false.
    logical :: implicit_zero_mode = .false.
    logical :: continuation_active = .false.
    logical :: photoelectron_active = .false.
    logical :: displacement_bounded = .false.
    integer(i32) :: electron_idx = 0_i32
    integer(i32) :: ion_idx = 0_i32
    integer(i32) :: photoelectron_idx = 0_i32
    integer(i32) :: search_direction = 0_i32
    integer(i32) :: iteration = 0_i32
    real(dp) :: plane_z
    real(dp) :: area
    real(dp) :: photoelectron_charge
    real(dp) :: photoelectron_emission_current_density
    real(dp) :: displacement_min
    real(dp) :: displacement_max
    real(dp) :: displacement_scale
    real(dp) :: feedback_reference(4)
    real(dp) :: displacement
    real(dp) :: displacement_before
    real(dp) :: guess(4)
    real(dp) :: observed(4)
    real(dp) :: response(6)
    real(dp) :: residual
    real(dp) :: return_flux
    real(dp) :: escape_flux
    type(matching_plane_response_provider_type) :: provider
    type(matching_plane_zhao_root_seed_type) :: root_committed
    type(matching_plane_zhao_root_seed_type) :: root_trial
    type(matching_plane_zhao_root_seed_type) :: root_candidate
    real(dp), allocatable :: moments(:, :)
    real(dp), allocatable :: reduce(:)
  contains
    procedure :: initialize
    procedure :: is_active
    procedure :: restore_gauge
    procedure :: begin_trial
    procedure :: prepare_iteration
    procedure :: finish_iteration
    procedure :: stage_stats
    procedure :: commit
  end type matching_plane_coupling_type

contains

  subroutine initialize(self, app, mesh, stats, mpi_ctx)
    class(matching_plane_coupling_type), intent(out) :: self
    type(app_config), intent(in) :: app
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(inout) :: stats
    type(mpi_context), intent(in) :: mpi_ctx
    integer(i32) :: matching_response_status
    character(len=512) :: matching_response_message
    logical :: implicit_zero_mode_supported
    real(dp) :: matching_response_input(5)

    self%active = trim(lower_ascii(app%surface_current%model)) == 'matching_plane_quasistatic'
    self%implicit_zero_mode = self%active .and. app%surface_current%implicit_zero_mode
    self%continuation_active = self%active .and. &
                               trim(lower_ascii(app%surface_current%zhao_root_selection)) == 'continuation'
    implicit_zero_mode_supported = .false.
    self%displacement_min = 0.0_dp
    self%displacement_max = 0.0_dp
    self%displacement_scale = 0.0_dp
    self%feedback_reference = 0.0_dp
    call self%provider%initialize( &
      app, mpi_ctx, matching_response_status, matching_response_message &
      )
    if (matching_response_status /= matching_plane_provider_ok) then
      error stop 'matching-plane response preflight failed: '//trim(matching_response_message)
    end if
    self%photoelectron_charge = 0.0_dp
    self%photoelectron_emission_current_density = 0.0_dp
    if (self%active) then
      call self%provider%get_matching_plane_z( &
        self%plane_z, matching_response_status, matching_response_message &
        )
      if (matching_response_status /= matching_plane_provider_ok) then
        error stop 'matching-plane response metadata failed: '//trim(matching_response_message)
      end if
      self%electron_idx = matching_species_index(app, app%surface_current%electron_species)
      self%ion_idx = matching_species_index(app, app%surface_current%ion_species)
      self%photoelectron_active = len_trim(app%surface_current%photoelectron_species) > 0
      if (self%photoelectron_active) then
        self%photoelectron_idx = matching_species_index(app, app%surface_current%photoelectron_species)
        self%photoelectron_charge = app%particle_species(self%photoelectron_idx)%q_particle
        self%photoelectron_emission_current_density = &
          app%particle_species(self%photoelectron_idx)%emit_current_density_a_m2
      end if
      self%area = product(app%sim%box_max(1:2) - app%sim%box_min(1:2))
      if (.not. ieee_is_finite(self%area) .or. self%area <= 0.0_dp) then
        error stop 'matching-plane area must be finite and positive.'
      end if
      if (self%implicit_zero_mode) then
        call self%provider%get_implicit_zero_mode_contract( &
          implicit_zero_mode_supported, self%displacement_bounded, self%displacement_min, &
          self%displacement_max, self%displacement_scale, self%feedback_reference &
          )
        if (.not. implicit_zero_mode_supported) then
          error stop 'implicit matching-plane zero mode requires a compatible table or Zhao online response.'
        end if
        if (.not. all(ieee_is_finite([ &
                                     app%sim%batch_duration, self%displacement_scale, self%feedback_reference, &
                                     app%particle_species(self%electron_idx)%q_particle, &
                                     app%particle_species(self%ion_idx)%q_particle, self%photoelectron_charge &
                                     ])) .or. app%sim%batch_duration <= 0.0_dp .or. &
            self%displacement_scale <= 0.0_dp) then
          error stop 'implicit matching-plane zero-mode preflight inputs are invalid.'
        end if
        if (self%displacement_bounded) then
          if (.not. all(ieee_is_finite([self%displacement_min, self%displacement_max])) .or. &
              self%displacement_min >= self%displacement_max) then
            error stop 'implicit matching-plane bounded displacement domain is invalid.'
          end if
        else
          select case (trim(lower_ascii(app%surface_current%zhao_branch)))
          case ('a', 'b')
            self%search_direction = 1_i32
          case ('c')
            self%search_direction = -1_i32
          case default
            self%search_direction = 0_i32
          end select
        end if
        if (self%photoelectron_active) then
          if (self%feedback_reference(1) < 0.0_dp .or. self%feedback_reference(2) <= 0.0_dp) then
            error stop 'implicit matching-plane PE mode requires nonnegative PE flux and positive PE energy.'
          end if
        else
          if (any(self%feedback_reference(1:2) /= 0.0_dp)) then
            error stop 'implicit matching-plane no-PE mode requires zero singleton PE flux and energy.'
          end if
        end if
      end if
      if (mesh%nelem > 0_i32) then
        if (app%sim%box_max(3) <= max( &
            maxval(mesh%v0(3, :)), maxval(mesh%v1(3, :)), maxval(mesh%v2(3, :)) &
            )) then
          error stop 'matching-plane gauge must lie strictly above every mesh vertex.'
        end if
      end if
      allocate (self%moments(4_i32, app%n_particle_species))
      allocate (self%reduce(4_i32*app%n_particle_species))
      if (stats%batches > 0_i32 .and. .not. stats%matching_plane_state_valid) then
        error stop 'matching-plane resume requires a schema-v9 committed coupling state.'
      end if
      if (self%continuation_active .and. stats%matching_plane_state_valid .and. mpi_is_root(mpi_ctx)) then
        matching_response_input = [ &
                                  stats%matching_plane_displacement_c_m2, stats%matching_plane_feedback &
                                  ]
        call self%provider%reconstruct_continuation_seed_local( &
          matching_response_input, stats%matching_plane_response, self%root_committed, &
          matching_response_status, matching_response_message &
          )
        if (matching_response_status /= matching_plane_provider_ok) then
          self%root_committed = matching_plane_zhao_root_seed_type()
          write (error_unit, '(a)') &
            'WARNING: matching-plane continuation resume seed was not reconstructed; '// &
            'the next endpoint will use a full multistart bootstrap: '//trim(matching_response_message)
          flush (error_unit)
        end if
      end if
    else
      stats%matching_plane_state_valid = .false.
    end if

  end subroutine initialize

  logical function is_active(self) result(active)
    class(matching_plane_coupling_type), intent(in) :: self
    active = self%active
  end function is_active

  subroutine restore_gauge(self, mesh, stats, snapshot)
    class(matching_plane_coupling_type), intent(in) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    if (self%active .and. stats%matching_plane_state_valid) then
      call snapshot%set_matching_plane_gauge(mesh, self%plane_z, stats%matching_plane_phi_v)
    end if
  end subroutine restore_gauge

  subroutine begin_trial(self, mesh, snapshot, stats)
    class(matching_plane_coupling_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(electrostatic_snapshot_type), intent(in) :: snapshot
    type(sim_stats), intent(in) :: stats

    if (self%continuation_active) self%root_trial = self%root_committed
    self%iteration = 0_i32
    self%residual = 0.0_dp
    self%return_flux = 0.0_dp
    self%escape_flux = 0.0_dp
    if (self%active) then
      if (self%implicit_zero_mode) then
        self%displacement_before = finite_charge_sum( &
                                   mesh%q_elem, 'implicit matching-plane charge before trial' &
                                   )/self%area
        if (.not. ieee_is_finite(self%displacement_before)) then
          error stop 'implicit matching-plane charge density before trial overflowed.'
        end if
      else
        self%displacement_before = snapshot%get_matching_plane_displacement()
      end if
      self%displacement = self%displacement_before
      if (self%implicit_zero_mode) then
        if (.not. self%displacement_bounded .and. stats%matching_plane_state_valid) then
          self%guess = stats%matching_plane_feedback
          self%guess(3:4) = 0.0_dp
          if (self%photoelectron_active .and. self%guess(2) <= 0.0_dp) then
            self%guess(2) = self%feedback_reference(2)
          end if
        else
          self%guess = self%feedback_reference
        end if
      else if (stats%matching_plane_state_valid) then
        self%guess = stats%matching_plane_feedback
      else
        self%guess = 0.0_dp
      end if
    end if
  end subroutine begin_trial

  subroutine prepare_iteration(self, app, mesh, snapshot, surface_closure, mpi_ctx, trial_batch_duration)
    class(matching_plane_coupling_type), intent(inout) :: self
    type(app_config), intent(in) :: app
    type(mesh_type), intent(in) :: mesh
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(surface_closure_contract_type), intent(inout) :: surface_closure
    type(mpi_context), intent(in) :: mpi_ctx
    real(dp), intent(in) :: trial_batch_duration
    real(dp) :: matching_displacement_seed, matching_response_input(5)
    integer(i32) :: matching_response_status
    character(len=512) :: matching_response_message

    if (.not. self%active) return
    self%iteration = self%iteration + 1_i32
    if (self%implicit_zero_mode) then
      matching_displacement_seed = self%displacement
      call solve_matching_implicit_zero_mode( &
        self%provider, mpi_ctx, self%displacement_before, matching_displacement_seed, &
        trial_batch_duration, &
        self%displacement_bounded, self%displacement_min, self%displacement_max, &
        self%displacement_scale, self%search_direction, self%guess, &
        app%particle_species(self%electron_idx)%q_particle, &
        app%particle_species(self%ion_idx)%q_particle, &
        self%photoelectron_active, self%photoelectron_charge, &
        self%root_trial, self%root_candidate, self%displacement, self%response &
        )
      if (self%continuation_active .and. self%root_committed%valid) then
        self%root_trial = self%root_candidate
      end if
    else
      matching_response_input = [self%displacement, self%guess]
      call self%provider%evaluate( &
        matching_response_input, mpi_ctx, self%response, &
        matching_response_status, matching_response_message &
        )
      if (matching_response_status /= matching_plane_provider_ok) then
        error stop 'matching-plane response evaluation failed: '//trim(matching_response_message)
      end if
    end if
    call configure_matching_surface_closure( &
      surface_closure, self%electron_idx, self%ion_idx, self%photoelectron_idx, &
      self%photoelectron_active, self%response, self%implicit_zero_mode, self%area, &
      self%guess, &
      app%particle_species(self%electron_idx)%q_particle, &
      app%particle_species(self%ion_idx)%q_particle, &
      self%photoelectron_charge, self%photoelectron_emission_current_density &
      )
    call snapshot%set_matching_plane_gauge(mesh, self%plane_z, self%response(1))
  end subroutine prepare_iteration

  logical function finish_iteration(self, app, mpi_ctx, moments_thread, trial_batch_duration, batch_idx) result(done)
    class(matching_plane_coupling_type), intent(inout) :: self
    type(app_config), intent(in) :: app
    type(mpi_context), intent(in) :: mpi_ctx
    real(dp), intent(in) :: moments_thread(:, :, :), trial_batch_duration
    integer(i32), intent(in) :: batch_idx
    logical :: matching_converged
    integer(i32) :: matching_response_status, matching_axis
    real(dp) :: matching_feedback_scales(4), matching_absolute_defects(4), matching_component_residuals(4)
    character(len=512) :: matching_response_message

    done = .true.
    if (.not. self%active) return
    self%moments = sum(moments_thread, dim=3)
    self%reduce = reshape(self%moments, [size(self%reduce)])
    call mpi_allreduce_sum_real_dp_array(mpi_ctx, self%reduce)
    self%moments = reshape(self%reduce, shape(self%moments))
    call resolve_matching_observed_feedback( &
      self%moments, self%electron_idx, self%ion_idx, self%photoelectron_idx, &
      self%photoelectron_active, self%area, trial_batch_duration, self%observed, &
      self%return_flux, self%escape_flux &
      )
    if (self%photoelectron_active .and. self%observed(1) == 0.0_dp) then
      ! Mean energy is undefined for an empty PE sample.  Preserve the
      ! current canonical energy instead of introducing a spurious zero.
      self%observed(2) = self%guess(2)
    end if
    call self%provider%validate_feedback( &
      self%observed, matching_response_status, matching_response_message &
      )
    if (matching_response_status /= matching_plane_provider_ok) then
      error stop 'matching-plane feedback validation failed: '//trim(matching_response_message)
    end if
    self%residual = self%provider%feedback_residual( &
                    self%guess, self%observed, app%surface_current%coupling_rtol, &
                    app%surface_current%coupling_atol &
                    )
    matching_converged = self%provider%feedback_converged( &
                         self%guess, self%observed, app%surface_current%coupling_rtol, &
                         app%surface_current%coupling_atol &
                         )
    if (matching_converged) return
    if (self%iteration >= app%surface_current%coupling_max_iterations) then
      if (mpi_is_root(mpi_ctx)) then
        call self%provider%get_feedback_scales(matching_feedback_scales)
        matching_absolute_defects = abs(self%observed - self%guess)
        matching_component_residuals = 0.0_dp
        do matching_axis = 1_i32, 4_i32
          if (matching_feedback_scales(matching_axis) <= 0.0_dp) cycle
          if (app%surface_current%coupling_atol(matching_axis) > &
              app%surface_current%coupling_rtol*matching_feedback_scales(matching_axis)) then
            matching_component_residuals(matching_axis) = app%surface_current%coupling_rtol* &
                                                          (matching_absolute_defects(matching_axis)/ &
                                                           app%surface_current%coupling_atol(matching_axis))
          else
            matching_component_residuals(matching_axis) = matching_absolute_defects(matching_axis)/ &
                                                          matching_feedback_scales(matching_axis)
          end if
        end do
        write (error_unit, '(a,i0,a,i0,a,es24.16,a,es24.16)') &
          'WARNING: accepting matching-plane nonconvergence: batch=', batch_idx, ', iterations=', self%iteration, &
          ', residual=', self%residual, ', rtol=', app%surface_current%coupling_rtol
        write (error_unit, '(a,4(1x,es24.16))') 'matching-plane guess=', self%guess
        write (error_unit, '(a,4(1x,es24.16))') 'matching-plane observed=', self%observed
        write (error_unit, '(a,4(1x,es24.16))') 'matching-plane feedback scales=', matching_feedback_scales
        write (error_unit, '(a,4(1x,es24.16))') &
          'matching-plane coupling atols=', app%surface_current%coupling_atol
        write (error_unit, '(a,4(1x,es24.16))') &
          'matching-plane absolute defects=', matching_absolute_defects
        write (error_unit, '(a,4(1x,es24.16))') &
          'matching-plane effective component residuals=', matching_component_residuals
        flush (error_unit)
      end if
      ! A finite replay is still a valid batch sample.  Preserve its residual
      ! as a convergence receipt and use the observed feedback to seed the
      ! next batch instead of discarding all committed progress.
      return
    end if
    self%guess = self%guess + app%surface_current%coupling_relaxation* &
                 (self%observed - self%guess)
    if (self%implicit_zero_mode .and. .not. self%displacement_bounded) self%guess(3:4) = 0.0_dp
    done = .false.
  end function finish_iteration

  subroutine stage_stats(self, stats_candidate)
    class(matching_plane_coupling_type), intent(in) :: self
    type(sim_stats), intent(inout) :: stats_candidate

    if (self%active) then
      stats_candidate%matching_plane_state_valid = .true.
      stats_candidate%matching_plane_displacement_c_m2 = self%displacement
      stats_candidate%matching_plane_phi_v = self%response(1)
      stats_candidate%matching_plane_response = self%response
      stats_candidate%matching_plane_feedback = self%observed
      stats_candidate%matching_plane_photoelectron_return_flux_m2_s = self%return_flux
      stats_candidate%matching_plane_photoelectron_escape_flux_m2_s = self%escape_flux
      stats_candidate%matching_plane_iterations = self%iteration
      stats_candidate%matching_plane_residual = self%residual
    end if
  end subroutine stage_stats

  subroutine commit(self, mesh, stats_candidate)
    class(matching_plane_coupling_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(inout) :: stats_candidate
    real(dp) :: matching_committed_displacement

    if (self%implicit_zero_mode) then
      matching_committed_displacement = finite_charge_sum( &
                                        mesh%q_elem, 'implicit matching-plane committed charge' &
                                        )/self%area
      if (.not. ieee_is_finite(matching_committed_displacement)) then
        error stop 'implicit matching-plane committed charge density overflowed.'
      end if
      stats_candidate%matching_plane_displacement_c_m2 = matching_committed_displacement
    end if
    if (self%continuation_active) self%root_committed = self%root_candidate
  end subroutine commit

  integer(i32) function matching_species_index(app, species_key) result(index)
    type(app_config), intent(in) :: app
    character(len=*), intent(in) :: species_key
    integer(i32) :: species_idx

    index = 0_i32
    do species_idx = 1_i32, app%n_particle_species
      if (trim(app%particle_species(species_idx)%species_key) /= trim(species_key)) cycle
      index = species_idx
      return
    end do
    error stop 'matching-plane species role was not found in particle species.'
  end function matching_species_index

  subroutine configure_matching_surface_closure( &
    contract, electron_idx, ion_idx, photoelectron_idx, photoelectron_active, response, implicit_zero_mode, area_m2, &
    feedback_reference, electron_charge, ion_charge, photoelectron_charge, &
    photoelectron_emission_current_density &
    )
    type(surface_closure_contract_type), intent(inout) :: contract
    integer(i32), intent(in) :: electron_idx, ion_idx, photoelectron_idx
    real(dp), intent(in) :: response(6)
    logical, intent(in) :: implicit_zero_mode, photoelectron_active
    real(dp), intent(in) :: area_m2, feedback_reference(4)
    real(dp), intent(in) :: electron_charge, ion_charge, photoelectron_charge
    real(dp), intent(in) :: photoelectron_emission_current_density

    real(dp) :: barrier_energy_ev, emission_flux, escape_flux, return_flux

    contract%active = .true.
    contract%has_absorbed_target = .false.
    contract%has_emission_target = .false.
    contract%has_escape_target = .false.
    contract%has_inflow_kinetic_map = .false.
    contract%has_outflow_kinetic_barrier = .false.
    contract%has_inflow_number_flux = .false.
    contract%absorbed_current_a = 0.0_dp
    contract%emission_current_a = 0.0_dp
    contract%escaped_particle_current_a = 0.0_dp
    contract%inflow_reservoir_potential_v = 0.0_dp
    contract%inflow_access_potential_v = 0.0_dp
    contract%inflow_kinetic_face = 0_i32
    contract%outflow_barrier_potential_v = 0.0_dp
    contract%outflow_barrier_face = 0_i32
    contract%inflow_number_flux_m2_s = 0.0_dp

    contract%has_inflow_number_flux([electron_idx, ion_idx]) = .true.
    contract%inflow_number_flux_m2_s(electron_idx) = response(2)
    contract%inflow_number_flux_m2_s(ion_idx) = response(3)
    contract%has_inflow_kinetic_map([electron_idx, ion_idx]) = .true.
    contract%inflow_access_potential_v(electron_idx) = response(4)
    contract%inflow_access_potential_v(ion_idx) = response(5)
    contract%inflow_kinetic_face([electron_idx, ion_idx]) = 6_i32
    ! The response inward fluxes already include all outer-sheath return of
    ! ambient species. Reflecting ambient outflow locally would count it twice.
    if (photoelectron_active) then
      contract%has_outflow_kinetic_barrier(photoelectron_idx) = .true.
      contract%outflow_barrier_potential_v(photoelectron_idx) = response(6)
      contract%outflow_barrier_face(photoelectron_idx) = 6_i32
    end if
    if (implicit_zero_mode) then
      contract%has_absorbed_target([electron_idx, ion_idx]) = .true.
      contract%absorbed_current_a(electron_idx) = electron_charge*response(2)*area_m2
      contract%absorbed_current_a(ion_idx) = ion_charge*response(3)*area_m2
      if (photoelectron_active) then
        barrier_energy_ev = response(1) - response(6)
        escape_flux = 0.0_dp
        if (feedback_reference(1) > 0.0_dp) then
          if (.not. ieee_is_finite(feedback_reference(2)) .or. feedback_reference(2) <= 0.0_dp) then
            error stop 'positive implicit matching-plane PE flux requires positive mean energy.'
          end if
          escape_flux = feedback_reference(1)*exp(-barrier_energy_ev/feedback_reference(2))
        end if
        emission_flux = photoelectron_emission_current_density/(-photoelectron_charge)
        return_flux = emission_flux - escape_flux
        if (.not. all(ieee_is_finite([barrier_energy_ev, emission_flux, escape_flux, return_flux])) .or. &
            barrier_energy_ev < 0.0_dp .or. escape_flux < 0.0_dp .or. return_flux < 0.0_dp) then
          error stop 'implicit matching-plane current targets are invalid.'
        end if
        contract%has_absorbed_target(photoelectron_idx) = .true.
        contract%absorbed_current_a(photoelectron_idx) = photoelectron_charge*return_flux*area_m2
        contract%has_emission_target(photoelectron_idx) = .true.
        contract%emission_current_a(photoelectron_idx) = photoelectron_emission_current_density*area_m2
        contract%has_escape_target(photoelectron_idx) = .true.
        contract%escaped_particle_current_a(photoelectron_idx) = photoelectron_charge*escape_flux*area_m2
      end if
    end if
  end subroutine configure_matching_surface_closure

  subroutine resolve_matching_observed_feedback( &
    moments, electron_idx, ion_idx, photoelectron_idx, photoelectron_active, area_m2, duration_s, &
    feedback, return_flux, escape_flux &
    )
    real(dp), intent(in) :: moments(:, :)
    integer(i32), intent(in) :: electron_idx, ion_idx, photoelectron_idx
    logical, intent(in) :: photoelectron_active
    real(dp), intent(in) :: area_m2, duration_s
    real(dp), intent(out) :: feedback(4), return_flux, escape_flux
    real(dp) :: normalization, photoelectron_outward_number

    if (any(.not. ieee_is_finite(moments)) .or. any(moments < 0.0_dp)) then
      error stop 'matching-plane particle moments are invalid.'
    end if
    normalization = area_m2*duration_s
    if (.not. ieee_is_finite(normalization) .or. normalization <= 0.0_dp) then
      error stop 'matching-plane flux normalization must be finite and positive.'
    end if
    feedback(1:2) = 0.0_dp
    return_flux = 0.0_dp
    escape_flux = 0.0_dp
    if (photoelectron_active) then
      photoelectron_outward_number = moments(1, photoelectron_idx)
      feedback(1) = photoelectron_outward_number/normalization
      if (photoelectron_outward_number > 0.0_dp) then
        feedback(2) = moments(2, photoelectron_idx)/(photoelectron_outward_number*qe)
      end if
      return_flux = moments(3, photoelectron_idx)/normalization
      escape_flux = moments(4, photoelectron_idx)/normalization
    end if
    feedback(3) = moments(1, electron_idx)/normalization
    feedback(4) = moments(1, ion_idx)/normalization
    if (.not. all(ieee_is_finite([feedback, return_flux, escape_flux]))) then
      error stop 'matching-plane observed feedback is not finite.'
    end if
  end subroutine resolve_matching_observed_feedback

end module bem_matching_plane_coupling
