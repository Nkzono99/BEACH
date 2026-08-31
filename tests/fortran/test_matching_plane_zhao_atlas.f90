!> Zhao solvability atlasがbranch別の失敗理由を失わず出力することを検証する。
program test_matching_plane_zhao_atlas
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, qe
  use bem_app_config, only: app_config, default_app_config, species_from_defaults
  use bem_matching_plane_zhao_atlas, only: generate_matching_plane_zhao_atlas, &
                                           matching_plane_atlas_ok, &
                                           matching_plane_atlas_invalid_grid, &
                                           matching_plane_zhao_atlas_query_csv_header, &
                                           matching_plane_zhao_atlas_csv_header
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, &
                          assert_equal_i32, delete_file_if_exists
  implicit none

  character(len=*), parameter :: query_path = 'test_matching_plane_zhao_atlas_queries.csv'
  character(len=*), parameter :: output_path = 'test_matching_plane_zhao_atlas.csv'
  type(app_config) :: cfg
  integer(i32) :: status
  integer :: unit_id, ios
  character(len=2048) :: line
  character(len=512) :: message
  logical :: output_exists

  call cleanup_files()
  call configure_online_fixture(cfg)
  call test_init(2)

  call test_begin('known_points_keep_independent_branch_statuses')
  call write_known_queries()
  call generate_matching_plane_zhao_atlas(cfg, query_path, output_path, status, message)
  call assert_equal_i32(status, matching_plane_atlas_ok, 'Zhao atlas generation failed: '//trim(message))
  open (newunit=unit_id, file=output_path, status='old', action='read', iostat=ios)
  call assert_true(ios == 0, 'Zhao atlas output was not published')
  if (ios == 0) then
    call read_and_assert(matching_plane_zhao_atlas_csv_header, 'atlas header changed')
    call read_and_assert(',A,no_physical_solution,', 'zero-field Zhao-A status changed')
    call read_and_assert(',B,ok,', 'zero-field Zhao-B was not certified')
    call read_and_assert(',C,no_physical_solution,', 'zero-field Zhao-C status changed')
    call read_and_assert(',A,', 'negative-field Zhao-A row is missing')
    call read_and_assert(',B,no_physical_solution,', 'negative-field Zhao-B silently fell back')
    call read_and_assert(',C,ok,', 'negative-field Zhao-C was not certified')
    call read_and_assert(',A,ok,', 'known Zhao-A point was not certified')
    call read_and_assert(',B,', 'known Zhao-A query lost the independent Zhao-B row')
    call read_and_assert(',C,', 'known Zhao-A query lost the independent Zhao-C row')
    call read_and_assert(',A,invalid_input,', 'invalid PE moment was not retained for Zhao-A')
    call read_and_assert(',B,invalid_input,', 'invalid PE moment was not retained for Zhao-B')
    call read_and_assert(',C,invalid_input,', 'invalid PE moment was not retained for Zhao-C')
    read (unit_id, '(a)', iostat=ios) line
    call assert_true(ios < 0, 'Zhao atlas wrote an unexpected extra row')
    close (unit_id)
  end if
  call test_end()

  call test_begin('query_csv_contract_is_strict')
  open (newunit=unit_id, file=query_path, status='replace', action='write')
  write (unit_id, '(a)') 'invalid_'//matching_plane_zhao_atlas_query_csv_header
  write (unit_id, '(a)') '0,0,0'
  close (unit_id)
  call delete_file_if_exists(output_path)
  call generate_matching_plane_zhao_atlas(cfg, query_path, output_path, status, message)
  call assert_equal_i32(status, matching_plane_atlas_invalid_grid, 'invalid atlas header was accepted')
  inquire (file=output_path, exist=output_exists)
  call assert_true(.not. output_exists, 'invalid atlas query published an output')
  call test_end()

  call cleanup_files()
  call test_summary()

contains

  subroutine configure_online_fixture(fixture_cfg)
    type(app_config), intent(out) :: fixture_cfg

    call default_app_config(fixture_cfg)
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
    fixture_cfg%particle_species(3)%temperature_ev = 2.2_dp
    fixture_cfg%particle_species(3)%has_temperature_ev = .true.
  end subroutine configure_online_fixture

  subroutine write_known_queries()
    real(dp), parameter :: type_a_displacement = 1.4187346568707933e-11_dp
    real(dp), parameter :: type_a_flux = 1.3754433596232731e13_dp

    open (newunit=unit_id, file=query_path, status='replace', action='write')
    write (unit_id, '(a)') matching_plane_zhao_atlas_query_csv_header
    write (unit_id, '(a)') '0,0,0'
    write (unit_id, '(*(g0,:,","))') - 0.02_dp*eps0, 0.0_dp, 0.0_dp
    write (unit_id, '(*(g0,:,","))') type_a_displacement, type_a_flux, 2.2_dp
    write (unit_id, '(a)') '0,1.0e10,0'
    close (unit_id)
  end subroutine write_known_queries

  subroutine read_and_assert(expected, failure_message)
    character(len=*), intent(in) :: expected, failure_message

    read (unit_id, '(a)', iostat=ios) line
    call assert_true(ios == 0 .and. index(trim(line), expected) > 0, failure_message)
  end subroutine read_and_assert

  subroutine cleanup_files()
    call delete_file_if_exists(query_path)
    call delete_file_if_exists(output_path)
    call delete_file_if_exists(output_path//'.beach-zhao-atlas.tmp')
  end subroutine cleanup_files

end program test_matching_plane_zhao_atlas
