!> Built-in Zhao evaluatorから既存形式のresponse CSVを生成できることを検証する。
program test_matching_plane_response_generator
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe
  use bem_app_config, only: app_config, default_app_config, species_from_defaults
  use bem_matching_plane_response, only: matching_plane_response_table_type, &
                                         matching_plane_response_query_csv_header, &
                                         matching_plane_response_ok, &
                                         get_matching_plane_response_snapshot, &
                                         reset_matching_plane_response_snapshot_cache
  use bem_matching_plane_response_generator, only: generate_matching_plane_zhao_response_table, &
                                                   matching_plane_generator_ok, &
                                                   matching_plane_generator_invalid_grid, &
                                                   matching_plane_generator_evaluation_failure
  use bem_matching_plane_response_provider, only: matching_plane_response_provider_type, &
                                                  matching_plane_provider_ok
  use bem_mpi, only: mpi_context
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, &
                          assert_equal_i32, assert_close_dp, delete_file_if_exists
  implicit none

  character(len=*), parameter :: query_path = 'test_matching_plane_zhao_queries.csv'
  character(len=*), parameter :: duplicate_query_path = 'test_matching_plane_zhao_queries_duplicate.csv'
  character(len=*), parameter :: control_query_path = 'test_matching_plane_zhao_queries_control.csv'
  character(len=*), parameter :: missing_zero_query_path = 'test_matching_plane_zhao_queries_missing_zero.csv'
  character(len=*), parameter :: evaluation_failure_query_path = 'test_matching_plane_zhao_queries_failure.csv'
  character(len=*), parameter :: output_path = 'test_matching_plane_zhao_generated.csv'
  type(app_config) :: cfg, failure_cfg
  type(matching_plane_response_table_type) :: table
  type(matching_plane_response_provider_type) :: provider
  type(mpi_context) :: serial_mpi
  real(dp) :: query(5), response(6), online_response(6), matching_plane_z_m
  real(dp) :: implicit_displacement_min, implicit_displacement_max, implicit_displacement_scale
  real(dp) :: implicit_feedback_reference(4)
  integer(i32) :: status
  integer :: unit_id, ios, query_index
  character(len=512) :: message
  character(len=64) :: preserved_line
  logical :: implicit_supported, implicit_displacement_bounded

  call cleanup_files()
  call configure_online_fixture(cfg)
  call test_init(4)

  call test_begin('generated_curve_roundtrips_through_table_backend')
  call write_nonzero_photoelectron_curve(query_path)
  call generate_matching_plane_zhao_response_table(cfg, query_path, output_path, status, message)
  call assert_equal_i32(status, matching_plane_generator_ok, 'Zhao response generation failed: '//trim(message))
  call reset_matching_plane_response_snapshot_cache()
  call get_matching_plane_response_snapshot(output_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_ok, 'generated response table did not load: '//trim(message))
  call table%get_matching_plane_z(matching_plane_z_m, status, message)
  call assert_equal_i32(status, matching_plane_response_ok, 'generated matching-plane metadata failed')
  call assert_close_dp(matching_plane_z_m, 1.25_dp, 0.0_dp, 'generated matching-plane height mismatch')
  serial_mpi = mpi_context()
  call provider%initialize(cfg, serial_mpi, status, message)
  call assert_equal_i32(status, matching_plane_provider_ok, 'online provider initialization failed: '//trim(message))
  call provider%get_implicit_zero_mode_contract( &
    implicit_supported, implicit_displacement_bounded, implicit_displacement_min, implicit_displacement_max, &
    implicit_displacement_scale, implicit_feedback_reference &
    )
  call assert_true(implicit_supported, 'online provider did not enable implicit zero mode')
  call assert_true(.not. implicit_displacement_bounded, 'online provider exposed a fixed displacement domain')
  call assert_true(implicit_displacement_scale > 0.0_dp, 'online provider displacement scale is not positive')
  call assert_true(implicit_feedback_reference(1) == 0.0_dp, 'online provider PE bootstrap flux mismatch')
  call assert_close_dp( &
    implicit_feedback_reference(2), 3.0_dp, 1.0e-9_dp, &
    'online provider PE bootstrap energy mismatch' &
    )
  do query_index = 1, 2
    query = [0.0_dp, real(query_index - 1, dp)*1.0e10_dp, 3.0_dp, 0.0_dp, 0.0_dp]
    call table%evaluate(query, response, status, message)
    call assert_equal_i32(status, matching_plane_response_ok, 'generated table node did not evaluate')
    call provider%evaluate(query, serial_mpi, online_response, status, message)
    call assert_equal_i32(status, matching_plane_provider_ok, 'online query did not evaluate')
    call assert_true( &
      all(abs(response - online_response) <= 1.0e-12_dp*max(1.0_dp, abs(online_response))), &
      'generated table node changed the online response' &
      )
  end do
  call test_end()

  call test_begin('evaluation_failure_preserves_existing_output_atomically')
  call delete_file_if_exists(output_path)
  open (newunit=unit_id, file=output_path, status='replace', action='write')
  write (unit_id, '(a)') 'preserved-response-sentinel'
  close (unit_id)
  call write_evaluation_failure_query(evaluation_failure_query_path)
  failure_cfg = cfg
  failure_cfg%surface_current%zhao_branch = 'b'
  call generate_matching_plane_zhao_response_table( &
    failure_cfg, evaluation_failure_query_path, output_path, status, message &
    )
  call assert_equal_i32( &
    status, matching_plane_generator_evaluation_failure, &
    'explicitly incompatible Zhao query did not fail table generation' &
    )
  preserved_line = ''
  open (newunit=unit_id, file=output_path, status='old', action='read', iostat=ios)
  call assert_true(ios == 0, 'pre-existing response output disappeared after evaluation failure')
  if (ios == 0) then
    read (unit_id, '(a)', iostat=ios) preserved_line
    close (unit_id)
  end if
  call assert_true( &
    ios == 0 .and. trim(preserved_line) == 'preserved-response-sentinel', &
    'evaluation failure replaced or truncated the pre-existing response output' &
    )
  call test_end()

  call test_begin('query_csv_header_and_numeric_grammar_are_strict')
  call write_invalid_header_query(control_query_path)
  call assert_invalid_query_does_not_publish( &
    control_query_path, 'non-exact query CSV header was accepted' &
    )
  call write_control_syntax_query(control_query_path)
  call assert_invalid_query_does_not_publish( &
    control_query_path, 'Fortran list-directed control syntax was accepted in a query row' &
    )
  call test_end()

  call test_begin('query_grid_requires_unique_cartesian_product_and_zero_flux')
  call write_duplicate_query_grid(duplicate_query_path)
  call assert_invalid_query_does_not_publish( &
    duplicate_query_path, 'duplicate coordinate and missing Cartesian point were accepted' &
    )
  call write_incomplete_query_grid(duplicate_query_path)
  call assert_invalid_query_does_not_publish( &
    duplicate_query_path, 'incomplete Cartesian query grid was accepted' &
    )
  call write_missing_zero_query(missing_zero_query_path)
  call assert_invalid_query_does_not_publish( &
    missing_zero_query_path, 'photoelectron-flux axis without a zero node was accepted' &
    )
  call test_end()

  call cleanup_files()
  call reset_matching_plane_response_snapshot_cache()
  call test_summary()

contains

  subroutine configure_online_fixture(fixture_cfg)
    type(app_config), intent(out) :: fixture_cfg

    call default_app_config(fixture_cfg)
    fixture_cfg%sim%box_max(3) = 1.25_dp
    fixture_cfg%surface_current%model = 'matching_plane_quasistatic'
    fixture_cfg%surface_current%response_backend = 'zhao_online'
    fixture_cfg%surface_current%zhao_branch = 'auto'
    fixture_cfg%surface_current%electron_species = 'electron'
    fixture_cfg%surface_current%ion_species = 'ion'
    fixture_cfg%surface_current%photoelectron_species = 'photoelectron'
    fixture_cfg%n_particle_species = 3_i32

    fixture_cfg%particle_species(1) = species_from_defaults()
    fixture_cfg%particle_species(1)%species_key = 'electron'
    fixture_cfg%particle_species(1)%m_particle = 9.1093837015e-31_dp
    fixture_cfg%particle_species(1)%drift_velocity(3) = -4.0529988897111727e5_dp
    fixture_cfg%particle_species(1)%temperature_ev = 12.0_dp
    fixture_cfg%particle_species(1)%has_temperature_ev = .true.

    fixture_cfg%particle_species(2) = species_from_defaults()
    fixture_cfg%particle_species(2)%species_key = 'ion'
    fixture_cfg%particle_species(2)%q_particle = qe
    fixture_cfg%particle_species(2)%m_particle = 1.67262192369e-27_dp
    fixture_cfg%particle_species(2)%number_density_m3 = 8.7e6_dp
    fixture_cfg%particle_species(2)%has_number_density_m3 = .true.
    fixture_cfg%particle_species(2)%drift_velocity(3) = -4.0529988897111727e5_dp
    fixture_cfg%particle_species(2)%temperature_ev = 0.1_dp
    fixture_cfg%particle_species(2)%has_temperature_ev = .true.

    fixture_cfg%particle_species(3) = species_from_defaults()
    fixture_cfg%particle_species(3)%species_key = 'photoelectron'
    fixture_cfg%particle_species(3)%m_particle = fixture_cfg%particle_species(1)%m_particle
    fixture_cfg%particle_species(3)%temperature_ev = 3.0_dp
    fixture_cfg%particle_species(3)%has_temperature_ev = .true.
  end subroutine configure_online_fixture

  subroutine write_duplicate_query_grid(path)
    character(len=*), intent(in) :: path
    integer :: unit_id

    open (newunit=unit_id, file=path, status='replace', action='write')
    write (unit_id, '(a)') matching_plane_response_query_csv_header
    write (unit_id, '(a)') '0,0,3.0,0,0'
    write (unit_id, '(a)') '1.0e-12,0,3.0,0,0'
    write (unit_id, '(a)') '0,1.0e10,3.0,0,0'
    write (unit_id, '(a)') '0,1.0e10,3.0,0,0'
    close (unit_id)
  end subroutine write_duplicate_query_grid

  subroutine write_incomplete_query_grid(path)
    character(len=*), intent(in) :: path
    integer :: unit_id

    open (newunit=unit_id, file=path, status='replace', action='write')
    write (unit_id, '(a)') matching_plane_response_query_csv_header
    write (unit_id, '(a)') '0,0,3.0,0,0'
    write (unit_id, '(a)') '1.0e-12,0,3.0,0,0'
    write (unit_id, '(a)') '0,1.0e10,3.0,0,0'
    close (unit_id)
  end subroutine write_incomplete_query_grid

  subroutine write_nonzero_photoelectron_curve(path)
    character(len=*), intent(in) :: path
    integer :: unit_id

    open (newunit=unit_id, file=path, status='replace', action='write')
    write (unit_id, '(a)') matching_plane_response_query_csv_header
    write (unit_id, '(a)') '0,1.0e10,3.0,0,0'
    write (unit_id, '(a)') '0,0,3.0,0,0'
    close (unit_id)
  end subroutine write_nonzero_photoelectron_curve

  subroutine write_invalid_header_query(path)
    character(len=*), intent(in) :: path
    integer :: unit_id

    open (newunit=unit_id, file=path, status='replace', action='write')
    write (unit_id, '(a)') 'invalid_'//matching_plane_response_query_csv_header
    write (unit_id, '(a)') '0,0,3.0,0,0'
    close (unit_id)
  end subroutine write_invalid_header_query

  subroutine write_control_syntax_query(path)
    character(len=*), intent(in) :: path
    integer :: unit_id

    open (newunit=unit_id, file=path, status='replace', action='write')
    write (unit_id, '(a)') matching_plane_response_query_csv_header
    write (unit_id, '(a)') '0,0,/,0,0'
    close (unit_id)
  end subroutine write_control_syntax_query

  subroutine write_missing_zero_query(path)
    character(len=*), intent(in) :: path
    integer :: unit_id

    open (newunit=unit_id, file=path, status='replace', action='write')
    write (unit_id, '(a)') matching_plane_response_query_csv_header
    write (unit_id, '(a)') '0,1.0e10,3.0,0,0'
    close (unit_id)
  end subroutine write_missing_zero_query

  subroutine write_evaluation_failure_query(path)
    character(len=*), intent(in) :: path
    integer :: unit_id

    open (newunit=unit_id, file=path, status='replace', action='write')
    write (unit_id, '(a)') matching_plane_response_query_csv_header
    write (unit_id, '(a)') '-1.0e-12,0,0,0,0'
    close (unit_id)
  end subroutine write_evaluation_failure_query

  subroutine cleanup_files()
    call delete_file_if_exists(query_path)
    call delete_file_if_exists(duplicate_query_path)
    call delete_file_if_exists(control_query_path)
    call delete_file_if_exists(missing_zero_query_path)
    call delete_file_if_exists(evaluation_failure_query_path)
    call delete_file_if_exists(output_path)
    call delete_file_if_exists(output_path//'.beach-zhao-response.tmp')
  end subroutine cleanup_files

  subroutine assert_invalid_query_does_not_publish(path, failure_message)
    character(len=*), intent(in) :: path, failure_message

    logical :: output_exists

    call delete_file_if_exists(output_path)
    call generate_matching_plane_zhao_response_table(cfg, path, output_path, status, message)
    call assert_equal_i32(status, matching_plane_generator_invalid_grid, failure_message)
    inquire (file=output_path, exist=output_exists)
    call assert_true(.not. output_exists, 'invalid query grid published a response output')
  end subroutine assert_invalid_query_does_not_publish

end program test_matching_plane_response_generator
