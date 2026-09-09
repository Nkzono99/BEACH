!> 応答 provider の設定解決と MPI rank 間の合意・配信を担う。
submodule(bem_matching_plane_response_provider) bem_matching_plane_response_provider_mpi
  use bem_constants, only: eps0, k_boltzmann, qe
  use bem_config_helpers, only: species_number_density_m3, species_temperature_k
  use bem_matching_plane_response, only: preflight_matching_plane_response_mpi
  use bem_mpi, only: mpi_is_root, mpi_allreduce_min_i32_scalar, &
                     mpi_allreduce_max_i32_scalar, mpi_allreduce_min_real_dp_array, &
                     mpi_allreduce_max_real_dp_array, mpi_bcast_i32_array, &
                     mpi_bcast_real_dp_array
  use bem_string_utils, only: lower_ascii
  implicit none
contains

  module procedure initialize_matching_plane_response_provider

  integer(i32) :: active_code, active_min, active_max, backend_min, backend_max
  integer(i32) :: height_valid, species_valid
  integer(i32) :: electron_idx, ion_idx, photoelectron_idx
  integer(i32) :: zhao_status, status_packet(1)
  logical :: photoelectron_active
  real(dp) :: electron_temperature_ev, photoelectron_temperature_ev
  real(dp) :: feedback_scales(4)
  character(len=512) :: backend_message
  character(len=16) :: backend_name
  real(dp) :: height_min(1), height_max(1), height_tolerance

  self%active = .false.
  self%backend = provider_backend_none
  self%matching_plane_z_m = 0.0_dp
  self%feedback_min = 0.0_dp
  self%feedback_max = 0.0_dp
  self%feedback_scale = 0.0_dp
  self%feedback_bounded = .false.
  self%implicit_displacement_min = 0.0_dp
  self%implicit_displacement_max = 0.0_dp
  self%implicit_displacement_scale = 0.0_dp
  self%implicit_feedback_reference = 0.0_dp
  self%implicit_displacement_bounded = .false.
  self%implicit_zero_mode_supported = .false.
  self%content_fingerprint = ''
  call accept_provider(status, message)
  self%active = trim(lower_ascii(cfg%surface_current%model)) == 'matching_plane_quasistatic'
  active_code = merge(1_i32, 0_i32, self%active)
  active_min = active_code
  active_max = active_code
  call mpi_allreduce_min_i32_scalar(mpi, active_min)
  call mpi_allreduce_max_i32_scalar(mpi, active_max)
  if (active_min /= active_max) then
    call reject_provider( &
      matching_plane_provider_invalid_argument, &
      'matching-plane activation differs across MPI ranks.', status, message &
      )
    return
  end if
  if (.not. self%active) return

  backend_name = trim(lower_ascii(cfg%surface_current%response_backend))
  select case (trim(backend_name))
  case ('table')
    self%backend = provider_backend_table
  case ('zhao_online')
    self%backend = provider_backend_zhao_online
  case default
    self%backend = provider_backend_none
  end select
  backend_min = self%backend
  backend_max = self%backend
  call mpi_allreduce_min_i32_scalar(mpi, backend_min)
  call mpi_allreduce_max_i32_scalar(mpi, backend_max)
  if (backend_min /= backend_max .or. backend_min == provider_backend_none) then
    call reject_provider( &
      matching_plane_provider_invalid_argument, &
      'matching-plane response backend is invalid or differs across MPI ranks.', status, message &
      )
    return
  end if

  self%matching_plane_z_m = cfg%sim%box_max(3)
  height_valid = merge(1_i32, 0_i32, ieee_is_finite(self%matching_plane_z_m))
  call mpi_allreduce_min_i32_scalar(mpi, height_valid)
  if (height_valid == 0_i32) then
    call reject_provider( &
      matching_plane_provider_invalid_argument, &
      'matching-plane height must be finite on every MPI rank.', status, message &
      )
    return
  end if
  height_min = [self%matching_plane_z_m]
  height_max = height_min
  call mpi_allreduce_min_real_dp_array(mpi, height_min)
  call mpi_allreduce_max_real_dp_array(mpi, height_max)
  height_tolerance = 128.0_dp*epsilon(1.0_dp)*max( &
                     1.0_dp, abs(height_min(1)), abs(height_max(1)) &
                     )
  if (height_max(1) - height_min(1) > height_tolerance) then
    call reject_provider( &
      matching_plane_provider_invalid_argument, &
      'matching-plane height differs across MPI ranks.', status, message &
      )
    return
  end if

  select case (self%backend)
  case (provider_backend_table)
    call preflight_matching_plane_response_mpi( &
      .true., trim(cfg%surface_current%response_table_path), mpi, &
      self%content_fingerprint, zhao_status, backend_message, table=self%table &
      )
    if (zhao_status /= matching_plane_response_ok) then
      call reject_provider( &
        matching_plane_provider_load_failure, trim(backend_message), status, message &
        )
      return
    end if
    call initialize_table_feedback_contract(self, status, message)
    if (status /= matching_plane_provider_ok) return
    call self%table%get_matching_plane_z(self%matching_plane_z_m, zhao_status, backend_message)
    if (zhao_status /= matching_plane_response_ok) then
      call reject_provider( &
        matching_plane_provider_load_failure, trim(backend_message), status, message &
        )
      return
    end if
    if (abs(self%matching_plane_z_m - cfg%sim%box_max(3)) > &
        128.0_dp*epsilon(1.0_dp)*max( &
        1.0_dp, abs(self%matching_plane_z_m), abs(cfg%sim%box_max(3)) &
        )) then
      call reject_provider( &
        matching_plane_provider_invalid_argument, &
        'matching-plane response height must equal domain.box_max[3].', status, message &
        )
      return
    end if

  case (provider_backend_zhao_online)
    electron_idx = provider_species_index(cfg, cfg%surface_current%electron_species)
    ion_idx = provider_species_index(cfg, cfg%surface_current%ion_species)
    photoelectron_active = len_trim(cfg%surface_current%photoelectron_species) > 0
    photoelectron_idx = 0_i32
    if (photoelectron_active) then
      photoelectron_idx = provider_species_index(cfg, cfg%surface_current%photoelectron_species)
    end if
    species_valid = merge( &
                    1_i32, 0_i32, min(electron_idx, ion_idx) > 0_i32 .and. &
                    (.not. photoelectron_active .or. photoelectron_idx > 0_i32) &
                    )
    call mpi_allreduce_min_i32_scalar(mpi, species_valid)
    if (species_valid == 0_i32) then
      call reject_provider( &
        matching_plane_provider_invalid_argument, &
        'matching-plane Zhao species role was not found on every MPI rank.', status, message &
        )
      return
    end if

    electron_temperature_ev = &
      species_temperature_k(cfg%particle_species(electron_idx))*k_boltzmann/qe
    photoelectron_temperature_ev = electron_temperature_ev
    if (photoelectron_active) then
      photoelectron_temperature_ev = &
        species_temperature_k(cfg%particle_species(photoelectron_idx))*k_boltzmann/qe
    end if
    zhao_status = matching_plane_zhao_ok
    backend_message = ''
    if (mpi_is_root(mpi)) then
      call self%zhao%initialize( &
        trim(cfg%surface_current%zhao_branch), &
        trim(cfg%surface_current%zhao_root_selection), &
        species_number_density_m3(cfg%particle_species(ion_idx)), &
        electron_temperature_ev, &
        -cfg%particle_species(electron_idx)%drift_velocity(3), &
        -cfg%particle_species(ion_idx)%drift_velocity(3), &
        cfg%particle_species(ion_idx)%m_particle, &
        cfg%particle_species(electron_idx)%m_particle, &
        photoelectron_temperature_ev, zhao_status, backend_message &
        )
    end if
    status_packet = [zhao_status]
    call mpi_bcast_i32_array(mpi, status_packet, 0_i32)
    zhao_status = status_packet(1)
    if (zhao_status /= matching_plane_zhao_ok) then
      if (.not. mpi_is_root(mpi)) backend_message = 'matching-plane Zhao initialization failed on MPI root.'
      call map_zhao_failure(zhao_status, backend_message, status, message)
      return
    end if

    feedback_scales = 0.0_dp
    if (mpi_is_root(mpi)) then
      call self%zhao%get_feedback_scales(feedback_scales, zhao_status, backend_message)
    end if
    status_packet = [zhao_status]
    call mpi_bcast_i32_array(mpi, status_packet, 0_i32)
    zhao_status = status_packet(1)
    if (zhao_status /= matching_plane_zhao_ok) then
      if (.not. mpi_is_root(mpi)) backend_message = 'matching-plane Zhao feedback contract failed on MPI root.'
      call map_zhao_failure(zhao_status, backend_message, status, message)
      return
    end if
    call mpi_bcast_real_dp_array(mpi, feedback_scales, 0_i32)
    if (any(.not. ieee_is_finite(feedback_scales)) .or. any(feedback_scales < 0.0_dp) .or. &
        feedback_scales(1) <= 0.0_dp .or. feedback_scales(2) <= 0.0_dp .or. &
        any(feedback_scales(3:4) /= 0.0_dp)) then
      call reject_provider( &
        matching_plane_provider_numerical_failure, &
        'matching-plane Zhao returned invalid feedback scales.', status, message &
        )
      return
    end if
    self%feedback_scale = feedback_scales
    self%feedback_bounded = .false.
    self%implicit_displacement_scale = sqrt( &
                                       eps0*species_number_density_m3(cfg%particle_species(ion_idx))* &
                                       qe*electron_temperature_ev &
                                       )
    if (.not. ieee_is_finite(self%implicit_displacement_scale) .or. &
        self%implicit_displacement_scale <= 0.0_dp) then
      call reject_provider( &
        matching_plane_provider_numerical_failure, &
        'matching-plane Zhao returned an invalid displacement scale.', status, message &
        )
      return
    end if
    self%implicit_feedback_reference = 0.0_dp
    if (photoelectron_active) then
      self%implicit_feedback_reference(1) = max( &
                                            0.0_dp, cfg%particle_species(photoelectron_idx)%emit_current_density_a_m2/ &
                                            abs(cfg%particle_species(photoelectron_idx)%q_particle) &
                                            )
      self%implicit_feedback_reference(2) = photoelectron_temperature_ev
    end if
    self%implicit_displacement_bounded = .false.
    self%implicit_zero_mode_supported = .true.
    self%content_fingerprint = 'zhao-online-v1'
  end select
  if (any(cfg%surface_current%coupling_atol > 0.0_dp .and. self%feedback_scale <= 0.0_dp)) then
    call reject_provider( &
      matching_plane_provider_invalid_argument, &
      'matching-plane coupling_atol must be zero on inactive feedback axes.', status, message &
      )
    return
  end if
  end procedure initialize_matching_plane_response_provider

  module procedure evaluate_matching_plane_response_provider

  real(dp) :: query_min(matching_plane_response_input_count)
  real(dp) :: query_max(matching_plane_response_input_count)
  real(dp) :: tolerance
  integer :: component
  integer(i32) :: local_status, query_valid, status_packet(1)
  character(len=512) :: backend_message

  output = 0.0_dp
  call accept_provider(status, message)
  query_valid = merge(1_i32, 0_i32, &
                      self%active .and. self%backend /= provider_backend_none .and. &
                      all(ieee_is_finite(input)) .and. all(input(2:5) >= 0.0_dp))
  call mpi_allreduce_min_i32_scalar(mpi, query_valid)
  if (query_valid == 0_i32) then
    call reject_provider( &
      matching_plane_provider_invalid_argument, &
      'matching-plane response provider or query is invalid on at least one MPI rank.', &
      status, message &
      )
    return
  end if

  query_min = input
  query_max = input
  call mpi_allreduce_min_real_dp_array(mpi, query_min)
  call mpi_allreduce_max_real_dp_array(mpi, query_max)
  do component = 1, matching_plane_response_input_count
    tolerance = 128.0_dp*epsilon(1.0_dp)*max( &
                1.0_dp, abs(query_min(component)), abs(query_max(component)) &
                )
    if (query_max(component) - query_min(component) > tolerance) then
      call reject_provider( &
        matching_plane_provider_invalid_argument, &
        'matching-plane response query differs across MPI ranks.', status, message &
        )
      return
    end if
  end do

  local_status = matching_plane_provider_ok
  backend_message = ''
  if (mpi_is_root(mpi)) then
    call self%evaluate_local(input, output, local_status, backend_message)
  end if
  status_packet = [local_status]
  call mpi_bcast_i32_array(mpi, status_packet, 0_i32)
  local_status = status_packet(1)
  if (local_status /= matching_plane_provider_ok) then
    if (.not. mpi_is_root(mpi)) backend_message = 'matching-plane response evaluation failed on MPI root.'
    call reject_provider(local_status, trim(backend_message), status, message)
    return
  end if
  call mpi_bcast_real_dp_array(mpi, output, 0_i32)
  end procedure evaluate_matching_plane_response_provider

  subroutine initialize_table_feedback_contract(self, status, message)
    class(matching_plane_response_provider_type), intent(inout) :: self
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    integer(i32), allocatable :: axis_sizes(:)
    real(dp), allocatable :: axis_values(:)
    integer(i32) :: axis, first, last

    call accept_provider(status, message)
    call self%table%get_fingerprint_data( &
      axis_sizes, axis_values, matching_plane_z_m=self%matching_plane_z_m &
      )
    if (size(axis_sizes) /= matching_plane_response_input_count .or. any(axis_sizes <= 0_i32) .or. &
        sum(axis_sizes) /= size(axis_values)) then
      call reject_provider( &
        matching_plane_provider_load_failure, &
        'matching-plane response returned invalid axis metadata.', status, message &
        )
      return
    end if

    self%implicit_displacement_min = minval(axis_values(:axis_sizes(1)))
    self%implicit_displacement_max = maxval(axis_values(:axis_sizes(1)))
    self%implicit_displacement_scale = max( &
                                       self%implicit_displacement_max - self%implicit_displacement_min, &
                                       abs(self%implicit_displacement_min), &
                                       abs(self%implicit_displacement_max) &
                                       )
    self%implicit_displacement_bounded = .true.

    first = 1_i32 + axis_sizes(1)
    do axis = 1_i32, 4_i32
      last = first + axis_sizes(axis + 1_i32) - 1_i32
      self%feedback_min(axis) = minval(axis_values(first:last))
      self%feedback_max(axis) = maxval(axis_values(first:last))
      self%implicit_feedback_reference(axis) = &
        0.5_dp*(self%feedback_min(axis) + self%feedback_max(axis))
      if (axis_sizes(axis + 1_i32) > 1_i32) then
        self%feedback_scale(axis) = self%feedback_max(axis) - self%feedback_min(axis)
        if (.not. ieee_is_finite(self%feedback_scale(axis)) .or. self%feedback_scale(axis) <= 0.0_dp) then
          call reject_provider( &
            matching_plane_provider_load_failure, &
            'matching-plane active feedback axis must have positive finite span.', status, message &
            )
          return
        end if
        self%feedback_bounded(axis) = .true.
      end if
      first = last + 1_i32
    end do
    self%implicit_zero_mode_supported = axis_sizes(1) > 1_i32 .and. &
                                        all(axis_sizes(2:5) == 1_i32) .and. &
                                        ( &
                                        (self%implicit_feedback_reference(1) > 0.0_dp .and. &
                                         self%implicit_feedback_reference(2) > 0.0_dp) .or. &
                                        all(self%implicit_feedback_reference(1:2) == 0.0_dp) &
                                        ) .and. &
                                        all(self%implicit_feedback_reference(3:4) == 0.0_dp)
  end subroutine initialize_table_feedback_contract

  integer(i32) function provider_species_index(cfg, species_key) result(index)
    type(app_config), intent(in) :: cfg
    character(len=*), intent(in) :: species_key
    integer(i32) :: species

    index = 0_i32
    do species = 1_i32, cfg%n_particle_species
      if (trim(cfg%particle_species(species)%species_key) /= trim(species_key)) cycle
      index = species
      return
    end do
  end function provider_species_index

  module procedure map_zhao_failure

  select case (zhao_status)
  case (matching_plane_zhao_invalid_argument)
    call reject_provider(matching_plane_provider_invalid_argument, zhao_message, status, message)
  case (matching_plane_zhao_no_physical_solution)
    call reject_provider(matching_plane_provider_no_physical_solution, zhao_message, status, message)
  case (matching_plane_zhao_ambiguous_solution)
    call reject_provider(matching_plane_provider_ambiguous_solution, zhao_message, status, message)
  case default
    call reject_provider(matching_plane_provider_numerical_failure, zhao_message, status, message)
  end select
  end procedure map_zhao_failure

  module procedure accept_provider

  status = matching_plane_provider_ok
  message = ''
  end procedure accept_provider

  module procedure reject_provider

  status = code
  message = ''
  message = trim(text)
  end procedure reject_provider

end submodule bem_matching_plane_response_provider_mpi
