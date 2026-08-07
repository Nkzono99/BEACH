!> 現行の局所 reservoir / closed PE 設定をFortran parserで検証する。
program test_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_types, only: bc_open, bc_reflect, bc_redistributed_reflect
  use bem_app_config_types, only: particle_inflow_reservoir
  use bem_app_config, only: app_config, default_app_config, load_app_config, &
                            particles_per_batch_from_config
  use bem_config_helpers, only: resolve_particle_boundaries
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, delete_file_if_exists
  implicit none

  type(app_config) :: cfg
  integer(i32) :: effective_boundary_low(3), effective_boundary_high(3)
  integer :: i
  character(len=64) :: run_mode
  character(len=512) :: probe_config_path
  character(len=*), parameter :: zhao_magnetized_path = 'test_zhao_magnetized_tmp.toml'
  character(len=*), parameter :: zhao_generic_barrier_path = 'test_zhao_generic_barrier_tmp.toml'
  character(len=*), parameter :: config_failure_path = 'test_zhao_config_failure_tmp.log'

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == '--config-failure-probe') then
    call get_command_argument(2, probe_config_path)
    call default_app_config(cfg)
    call load_app_config(trim(probe_config_path), cfg)
    error stop 'invalid Zhao config probe unexpectedly completed'
  end if

  call test_init(8)

  call test_begin('default_config')
  call default_app_config(cfg)
  call assert_true(trim(cfg%sim%field_solver) == 'auto', 'default field solver mismatch')
  call assert_true(trim(cfg%sim%field_bc_mode) == 'free', 'default field boundary mismatch')
  call assert_true(trim(cfg%sim%reservoir_potential_model) == 'none', 'default inflow model mismatch')
  call assert_true(trim(cfg%sim%open_boundary_model) == 'escape', 'default open model mismatch')
  call assert_equal_i32(cfg%checkpoint_stride, 0_i32, 'default checkpoint stride mismatch')
  call test_end()

  call test_begin('zhao_rejects_magnetized_closure')
  call write_zhao_variant(zhao_magnetized_path, 'b0 = [0.0, 0.0, 1.0e-9]', .false.)
  call assert_config_rejected(zhao_magnetized_path, 'requires sim.b0=[0,0,0]')
  call delete_file_if_exists(zhao_magnetized_path)
  call test_end()

  call test_begin('zhao_rejects_generic_reservoir_barrier')
  call write_zhao_variant(zhao_generic_barrier_path, '', .true.)
  call assert_config_rejected(zhao_generic_barrier_path, 'cannot be combined with the generic reservoir potential model')
  call delete_file_if_exists(zhao_generic_barrier_path)
  call test_end()

  call test_begin('zhao_fixed_current_config')
  call default_app_config(cfg)
  call load_app_config('examples/periodic2_zhao_fixed_current.toml', cfg)
  call assert_true(trim(cfg%surface_current%model) == 'zhao_stationary', 'Zhao current model mismatch')
  call assert_true(trim(cfg%surface_current%zhao_branch) == 'auto', 'Zhao branch mismatch')
  call assert_true( &
    all([(trim(cfg%particle_species(i)%surface_charge_closure) == 'fixed_current', i=1, 3)]), &
    'Zhao species must use fixed_current' &
    )
  call resolve_particle_boundaries( &
    cfg%sim, cfg%particle_boundary_low, cfg%particle_boundary_high, cfg%particle_species(3), &
    effective_boundary_low, effective_boundary_high &
    )
  call assert_equal_i32(effective_boundary_high(3), bc_open, 'Zhao PE z-high boundary must be open')
  call test_end()

  call test_begin('fixed_current_config')
  call default_app_config(cfg)
  call load_app_config('tests/fortran/fixed_current.toml', cfg)
  call assert_true( &
    trim(cfg%particle_species(1)%surface_charge_closure) == 'fixed_current', &
    'fixed-current closure mismatch' &
    )
  call assert_true(cfg%particle_species(1)%has_target_absorbed_current_a, 'fixed absorbed target presence mismatch')
  call assert_close_dp( &
    cfg%particle_species(1)%target_absorbed_current_a, -2.0_dp, 1.0e-15_dp, &
    'fixed absorbed target mismatch' &
    )
  call test_end()

  call test_begin('tutorial_config')
  call default_app_config(cfg)
  call load_app_config('examples/tutorial_insulator.toml', cfg)
  call assert_true(trim(cfg%sim%field_solver) == 'fmm', 'tutorial field solver mismatch')
  call assert_true(trim(cfg%sim%field_bc_mode) == 'periodic2', 'tutorial field boundary mismatch')
  call assert_equal_i32(cfg%n_particle_species, 1_i32, 'tutorial species count mismatch')
  call assert_equal_i32(particles_per_batch_from_config(cfg), 1_i32, 'tutorial batch particle count mismatch')
  call test_end()

  call test_begin('closed_photoelectron_config')
  call default_app_config(cfg)
  call load_app_config('examples/periodic2_closed_photoelectron.toml', cfg)
  call assert_equal_i32(cfg%n_particle_species, 3_i32, 'closed case species count mismatch')
  call assert_true(trim(cfg%sim%reservoir_potential_model) == 'none', 'closed case must use source VDF inflow')
  call assert_true(trim(cfg%sim%open_boundary_model) == 'escape', 'closed case must use ordinary escape')
  call assert_true(trim(cfg%particle_species(1)%source_mode) == 'volume_seed', 'electron source mode mismatch')
  call assert_equal_i32( &
    cfg%particle_species(1)%boundary_inflow_high(3), particle_inflow_reservoir, &
    'electron z-high boundary inflow mismatch' &
    )
  call assert_true(trim(cfg%particle_species(2)%source_mode) == 'volume_seed', 'ion source mode mismatch')
  call assert_equal_i32( &
    cfg%particle_species(2)%boundary_inflow_high(3), particle_inflow_reservoir, &
    'ion z-high boundary inflow mismatch' &
    )
  call assert_true(trim(cfg%particle_species(3)%source_mode) == 'photo_raycast', 'photoelectron source mode mismatch')
  call assert_equal_i32(cfg%particle_species(3)%boundary_high(3), bc_reflect, 'photoelectron top boundary mismatch')
  call assert_true( &
    trim(cfg%particle_species(3)%surface_charge_closure) == 'neutral_return', &
    'photoelectron closure mismatch' &
    )
  call assert_close_dp(cfg%particle_species(3)%temperature_ev, 1.5_dp, 1.0e-15_dp, 'photoelectron temperature mismatch')
  call test_end()

  call test_begin('all_particle_boundary_faces')
  call default_app_config(cfg)
  call load_app_config('tests/fortran/particle_boundary_faces.toml', cfg)
  call assert_true(all(cfg%particle_boundary_low == [bc_open, bc_reflect, bc_open]), 'global low faces mismatch')
  call assert_true( &
    all(cfg%particle_boundary_high == [bc_reflect, bc_open, bc_redistributed_reflect]), &
    'global high faces mismatch' &
    )
  call assert_equal_i32(cfg%particle_species(1)%boundary_low(1), bc_reflect, 'species x-low mismatch')
  call assert_equal_i32(cfg%particle_species(1)%boundary_high(1), bc_open, 'species x-high mismatch')
  call assert_equal_i32(cfg%particle_species(1)%boundary_high(2), bc_reflect, 'species y-high mismatch')
  call assert_equal_i32(cfg%particle_species(1)%boundary_low(3), bc_reflect, 'species z-low mismatch')
  call resolve_particle_boundaries( &
    cfg%sim, cfg%particle_boundary_low, cfg%particle_boundary_high, cfg%particle_species(1), &
    effective_boundary_low, effective_boundary_high &
    )
  call assert_true(all(effective_boundary_low == bc_reflect), 'effective low faces mismatch')
  call assert_true( &
    all(effective_boundary_high == [bc_open, bc_reflect, bc_redistributed_reflect]), &
    'effective high faces mismatch' &
    )
  call assert_equal_i32(cfg%checkpoint_stride, 2_i32, 'checkpoint stride mismatch')
  call test_end()

  call test_summary()

contains

  subroutine write_zhao_variant(path, sim_line, replace_reservoir)
    character(len=*), intent(in) :: path, sim_line
    logical, intent(in) :: replace_reservoir
    character(len=1024) :: line
    integer :: source_unit, output_unit, ios
    logical :: inserted_sim, replaced_reservoir

    inserted_sim = len_trim(sim_line) == 0
    replaced_reservoir = .not. replace_reservoir
    open (newunit=source_unit, file='examples/periodic2_zhao_fixed_current.toml', &
          status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open Zhao example fixture'
    open (newunit=output_unit, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create Zhao invalid-config fixture'
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (replace_reservoir .and. trim(line) == 'inflow_model = "source_vdf"') then
        write (output_unit, '(a)') 'inflow_model = "infinity_barrier"'
        replaced_reservoir = .true.
      else
        write (output_unit, '(a)') trim(line)
      end if
      if (.not. inserted_sim .and. trim(line) == '[sim]') then
        write (output_unit, '(a)') trim(sim_line)
        inserted_sim = .true.
      end if
    end do
    close (source_unit)
    close (output_unit)
    if (.not. inserted_sim .or. .not. replaced_reservoir) then
      error stop 'failed to specialize Zhao invalid-config fixture'
    end if
  end subroutine write_zhao_variant

  subroutine assert_config_rejected(path, expected_fragment)
    character(len=*), intent(in) :: path, expected_fragment
    character(len=1024) :: executable_path, command, child_line
    integer :: child_exit_status, child_cmd_status, child_unit, child_ios
    logical :: saw_expected

    call get_command_argument(0, executable_path)
    call delete_file_if_exists(config_failure_path)
    command = '"'//trim(executable_path)//'" --config-failure-probe "'//trim(path)//'" > '// &
              config_failure_path//' 2>&1'
    call execute_command_line(trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status)
    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'Zhao config failure probe command status mismatch')
    call assert_true(child_exit_status /= 0, 'invalid Zhao config must be rejected')

    saw_expected = .false.
    open (newunit=child_unit, file=config_failure_path, status='old', action='read', iostat=child_ios)
    if (child_ios /= 0) error stop 'failed to read Zhao config failure probe output'
    do
      read (child_unit, '(A)', iostat=child_ios) child_line
      if (child_ios /= 0) exit
      saw_expected = saw_expected .or. index(child_line, trim(expected_fragment)) > 0
    end do
    close (child_unit)
    call delete_file_if_exists(config_failure_path)
    call assert_true(saw_expected, 'Zhao config failure message mismatch: '//trim(expected_fragment))
  end subroutine assert_config_rejected
end program test_app_config_parser
