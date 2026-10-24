!> multiple_box_events の soft_discard（件数・割合・電荷の閾値、再開時の累積、初期統計の検査）を検証する。
program test_soft_discard
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config
  use bem_types, only: mesh_type, sim_stats, bc_reflect
  use bem_charge_ledger, only: charge_ledger_type
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_equal_i64, assert_close_dp, delete_file_if_exists, &
                          run_self_probe, file_contains
  implicit none

  character(len=*), parameter :: soft_discard_fraction_limit_path = 'test_soft_discard_fraction_limit_tmp.log'
  character(len=*), parameter :: soft_discard_charge_limit_path = 'test_soft_discard_charge_limit_tmp.log'
  character(len=*), parameter :: soft_discard_resume_limit_path = 'test_soft_discard_resume_limit_tmp.log'
  character(len=*), parameter :: soft_discard_invalid_state_path = 'test_soft_discard_invalid_state_tmp.log'
  character(len=64) :: run_mode

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == '--soft-discard-fraction-limit-probe') then
    call run_multiple_box_event_soft_discard_limit_probe(.false.)
    error stop 'soft-discard fraction-limit probe unexpectedly completed'
  else if (trim(run_mode) == '--soft-discard-charge-limit-probe') then
    call run_multiple_box_event_soft_discard_limit_probe(.true.)
    stop 0
  else if (trim(run_mode) == '--soft-discard-resume-limit-probe') then
    call run_multiple_box_event_soft_discard_resume_limit_probe()
    error stop 'soft-discard resume-limit probe unexpectedly completed'
  else if (trim(run_mode) == '--soft-discard-invalid-state-probe') then
    call run_multiple_box_event_soft_discard_invalid_state_probe()
    error stop 'soft-discard invalid-state probe unexpectedly completed'
  end if

  call test_init(1)

  call test_begin('multiple_box_event_soft_discard')
  call test_multiple_box_event_soft_discard()
  call test_end()

  call test_summary()

contains

  subroutine test_multiple_box_event_soft_discard()
    type(mesh_type) :: soft_mesh
    type(app_config) :: soft_cfg
    type(sim_stats) :: soft_stats
    type(charge_ledger_type) :: soft_ledger

    call configure_multiple_box_event_soft_discard_fixture( &
      soft_mesh, soft_cfg, 1_i32, 100000_i32, 1.0_dp, 1.0e9_dp &
      )
    soft_cfg%sim%multiple_box_events_retry_backend = 'upper_panel_fourier'

    call run_absorption_insulator(soft_mesh, soft_cfg, soft_stats, charge_ledger=soft_ledger)
    call assert_equal_i64(soft_stats%processed_particles, 1_i64, 'soft discard processed count mismatch')
    call assert_equal_i64(soft_stats%absorbed, 0_i64, 'soft discard absorbed count mismatch')
    call assert_equal_i64(soft_stats%escaped, 1_i64, 'soft discard escaped umbrella count mismatch')
    call assert_equal_i64( &
      soft_stats%multiple_box_events_retry_attempted, 1_i64, 'upper Fourier retry attempt count mismatch' &
      )
    call assert_equal_i64( &
      soft_stats%multiple_box_events_retry_resolved, 0_i64, 'ineligible upper Fourier retry must remain unresolved' &
      )
    call assert_equal_i64( &
      soft_stats%multiple_box_events_soft_discarded, 1_i64, 'soft discard event count mismatch' &
      )
    call assert_close_dp( &
      soft_stats%multiple_box_events_soft_discarded_abs_charge, 6.0_dp, 1.0e-12_dp, &
      'soft discard absolute charge mismatch' &
      )
    call assert_equal_i64( &
      soft_ledger%discarded_unresolved_count(1), 1_i64, 'soft discard ledger count mismatch' &
      )
    call assert_close_dp( &
      soft_ledger%discarded_unresolved(1), 6.0_dp, 1.0e-12_dp, 'soft discard ledger charge mismatch' &
      )
    call assert_close_dp(soft_ledger%residual(), 0.0_dp, 1.0e-12_dp, 'soft discard ledger residual mismatch')
    call test_multiple_box_event_soft_discard_threshold_equalities()
    call test_multiple_box_event_soft_discard_limits_and_context()
  end subroutine test_multiple_box_event_soft_discard

  subroutine test_multiple_box_event_soft_discard_threshold_equalities()
    type(mesh_type) :: threshold_mesh
    type(app_config) :: threshold_cfg
    type(sim_stats) :: threshold_stats

    call configure_multiple_box_event_soft_discard_fixture( &
      threshold_mesh, threshold_cfg, 1_i32, 1_i32, 0.5_dp, 1.0e9_dp &
      )
    call run_absorption_insulator(threshold_mesh, threshold_cfg, threshold_stats)
    call assert_equal_i64( &
      threshold_stats%multiple_box_events_soft_discarded, 1_i64, &
      'soft-discard count equal to grace must be accepted' &
      )

    call configure_multiple_box_event_soft_discard_fixture( &
      threshold_mesh, threshold_cfg, 2_i32, 0_i32, 1.0_dp, 1.0e9_dp &
      )
    call run_absorption_insulator(threshold_mesh, threshold_cfg, threshold_stats)
    call assert_equal_i64( &
      threshold_stats%multiple_box_events_soft_discarded, 2_i64, &
      'soft-discard fraction equal to limit must be accepted' &
      )
  end subroutine test_multiple_box_event_soft_discard_threshold_equalities

  subroutine test_multiple_box_event_soft_discard_limits_and_context()
    call assert_soft_discard_limit_probe( &
      '--soft-discard-fraction-limit-probe', soft_discard_fraction_limit_path, .false. &
      )
    call assert_soft_discard_limit_probe( &
      '--soft-discard-charge-limit-probe', soft_discard_charge_limit_path, .true. &
      )
    call assert_soft_discard_resume_limit_probe()
    call assert_soft_discard_invalid_state_probe()
  end subroutine test_multiple_box_event_soft_discard_limits_and_context

  subroutine assert_soft_discard_limit_probe(run_mode, output_path, charge_limit_probe)
    character(len=*), intent(in) :: run_mode, output_path
    logical, intent(in) :: charge_limit_probe
    character(len=1024) :: child_line
    integer :: child_exit_status, child_unit, child_ios
    integer(i32) :: accepted_summary_count
    logical :: saw_batch, saw_limit_summary, saw_expected_count, saw_expected_grace
    logical :: saw_processed, saw_fraction, saw_fraction_limit, saw_abs_charge, saw_abs_charge_limit

    call run_self_probe(run_mode, output_path, child_exit_status, env='OMP_NUM_THREADS=2')
    if (charge_limit_probe) then
      call assert_equal_i32(int(child_exit_status, i32), 0_i32, 'soft-discard charge warning stopped the run')
    else
      call assert_true(child_exit_status /= 0, 'soft-discard fraction-limit probe should stop the run')
    end if

    saw_batch = .false.
    saw_limit_summary = .false.
    saw_expected_count = .false.
    saw_expected_grace = .false.
    saw_processed = .false.
    saw_fraction = .false.
    saw_fraction_limit = .false.
    saw_abs_charge = .false.
    saw_abs_charge_limit = .false.
    accepted_summary_count = 0_i32
    open (newunit=child_unit, file=output_path, status='old', action='read', iostat=child_ios)
    if (child_ios /= 0) error stop 'failed to read soft-discard limit probe output'
    do
      read (child_unit, '(A)', iostat=child_ios) child_line
      if (child_ios /= 0) exit
      saw_batch = saw_batch .or. index(child_line, 'batch=1') > 0
      if (charge_limit_probe) then
        saw_limit_summary = saw_limit_summary .or. &
                            index(child_line, 'cumulative absolute-charge threshold crossed') > 0
      else
        saw_limit_summary = saw_limit_summary .or. &
                            index(child_line, 'soft-discard cumulative fraction limit exceeded') > 0
      end if
      saw_processed = saw_processed .or. index(child_line, ' processed=') > 0
      saw_fraction = saw_fraction .or. index(child_line, ' fraction=') > 0
      saw_fraction_limit = saw_fraction_limit .or. index(child_line, ' particles.tracking.events.discard_fraction_limit=') > 0
      saw_abs_charge = saw_abs_charge .or. index(child_line, 'abs_charge_C=') > 0
      saw_abs_charge_limit = saw_abs_charge_limit .or. &
                             (index(child_line, 'abs_charge_limit_C=') > 0 .or. &
                              index(child_line, 'discard_charge_warning_c=') > 0)
      if (index(child_line, 'multiple_box_events soft discard accepted:') > 0) then
        accepted_summary_count = accepted_summary_count + 1_i32
      end if
      if (charge_limit_probe) then
        saw_expected_count = saw_expected_count .or. index(child_line, ' count=1 ') > 0
        saw_expected_grace = saw_expected_grace .or. index(child_line, ' particles.tracking.events.discard_count_grace=100 ') > 0
      else
        saw_expected_count = saw_expected_count .or. index(child_line, ' count=2 ') > 0
        saw_expected_grace = saw_expected_grace .or. index(child_line, ' particles.tracking.events.discard_count_grace=1 ') > 0
      end if
    end do
    close (child_unit)
    call delete_file_if_exists(output_path)

    call assert_true(saw_batch, 'soft-discard limit summary is missing its batch index')
    if (charge_limit_probe) then
      call assert_equal_i32(accepted_summary_count, 1_i32, 'charge-warning batch was not committed')
    else
      call assert_equal_i32(accepted_summary_count, 0_i32, 'over-fraction batch emitted an accepted summary')
    end if
    call assert_true(saw_limit_summary, 'soft-discard cumulative-limit summary is missing')
    if (charge_limit_probe) then
      call assert_true( &
        saw_abs_charge .and. saw_abs_charge_limit, &
        'soft-discard absolute-charge warning is incomplete' &
        )
    else
      call assert_true( &
        saw_expected_count .and. saw_expected_grace .and. saw_processed .and. saw_fraction .and. &
        saw_fraction_limit .and. saw_abs_charge .and. saw_abs_charge_limit, &
        'soft-discard cumulative-fraction summary is incomplete' &
        )
    end if
  end subroutine assert_soft_discard_limit_probe

  subroutine assert_soft_discard_resume_limit_probe()
    character(len=1024) :: child_line
    integer :: child_exit_status, child_unit, child_ios
    logical :: saw_context, saw_accepted

    call run_self_probe('--soft-discard-resume-limit-probe', soft_discard_resume_limit_path, child_exit_status, &
                        env='OMP_NUM_THREADS=2')
    call assert_true(child_exit_status /= 0, 'soft-discard resume-limit probe should terminate with nonzero status')

    saw_context = .false.
    saw_accepted = .false.
    open (newunit=child_unit, file=soft_discard_resume_limit_path, status='old', action='read', iostat=child_ios)
    if (child_ios /= 0) error stop 'failed to read soft-discard resume-limit probe output'
    do
      read (child_unit, '(A)', iostat=child_ios) child_line
      if (child_ios /= 0) exit
      saw_context = saw_context .or. &
                    (index(child_line, 'soft-discard cumulative fraction limit exceeded: batch=5') > 0 .and. &
                     index(child_line, ' count=2 ') > 0 .and. &
                     index(child_line, ' particles.tracking.events.discard_count_grace=1 ') > 0 .and. &
                     index(child_line, ' processed=100 ') > 0 .and. index(child_line, ' fraction=') > 0 .and. &
                     index(child_line, ' particles.tracking.events.discard_fraction_limit=') > 0)
      saw_accepted = saw_accepted .or. index(child_line, 'multiple_box_events soft discard accepted:') > 0
    end do
    close (child_unit)
    call delete_file_if_exists(soft_discard_resume_limit_path)
    call assert_true(saw_context, 'resume soft-discard limit summary is incomplete')
    call assert_true(.not. saw_accepted, 'resume guard must not emit an accepted summary')
  end subroutine assert_soft_discard_resume_limit_probe

  subroutine assert_soft_discard_invalid_state_probe()
    integer :: exit_status
    logical :: saw_preflight

    call run_self_probe('--soft-discard-invalid-state-probe', soft_discard_invalid_state_path, exit_status, &
                        env='OMP_NUM_THREADS=2')
    call assert_true(exit_status /= 0, 'soft-discard invalid-state probe should terminate with nonzero status')
    saw_preflight = file_contains(soft_discard_invalid_state_path, 'soft-discard initial statistics are inconsistent')
    call delete_file_if_exists(soft_discard_invalid_state_path)
    call assert_true(saw_preflight, 'soft-discard initial-state preflight diagnostic is missing')
  end subroutine assert_soft_discard_invalid_state_probe

  subroutine run_multiple_box_event_soft_discard_limit_probe(charge_limit_probe)
    logical, intent(in) :: charge_limit_probe
    type(mesh_type) :: soft_mesh
    type(app_config) :: soft_cfg
    type(sim_stats) :: soft_stats

    if (charge_limit_probe) then
      call configure_multiple_box_event_soft_discard_fixture( &
        soft_mesh, soft_cfg, 1_i32, 100_i32, 1.0_dp, 5.0_dp &
        )
    else
      call configure_multiple_box_event_soft_discard_fixture( &
        soft_mesh, soft_cfg, 2_i32, 1_i32, 0.5_dp, 1.0e9_dp &
        )
    end if
    call run_absorption_insulator(soft_mesh, soft_cfg, soft_stats)
  end subroutine run_multiple_box_event_soft_discard_limit_probe

  subroutine run_multiple_box_event_soft_discard_resume_limit_probe()
    type(mesh_type) :: soft_mesh
    type(app_config) :: soft_cfg
    type(sim_stats) :: resume_seed, soft_stats

    call configure_multiple_box_event_soft_discard_fixture( &
      soft_mesh, soft_cfg, 1_i32, 1_i32, 0.01_dp, 1.0e9_dp &
      )
    soft_cfg%sim%batch_count = 6_i32
    resume_seed = sim_stats()
    resume_seed%batches = 5_i32
    resume_seed%processed_particles = 100_i64
    resume_seed%multiple_box_events_soft_discarded = 2_i64
    call run_absorption_insulator(soft_mesh, soft_cfg, soft_stats, initial_stats=resume_seed)
  end subroutine run_multiple_box_event_soft_discard_resume_limit_probe

  subroutine run_multiple_box_event_soft_discard_invalid_state_probe()
    type(mesh_type) :: soft_mesh
    type(app_config) :: soft_cfg
    type(sim_stats) :: invalid_seed, soft_stats

    call configure_multiple_box_event_soft_discard_fixture( &
      soft_mesh, soft_cfg, 1_i32, 100_i32, 1.0_dp, 1.0e9_dp &
      )
    invalid_seed = sim_stats()
    invalid_seed%processed_particles = 1_i64
    invalid_seed%multiple_box_events_soft_discarded = 2_i64
    call run_absorption_insulator(soft_mesh, soft_cfg, soft_stats, initial_stats=invalid_seed)
  end subroutine run_multiple_box_event_soft_discard_invalid_state_probe

  subroutine configure_multiple_box_event_soft_discard_fixture( &
    soft_mesh, soft_cfg, particle_count, count_grace, fraction_limit, abs_charge_limit &
    )
    type(mesh_type), intent(out) :: soft_mesh
    type(app_config), intent(out) :: soft_cfg
    integer(i32), intent(in) :: particle_count, count_grace
    real(dp), intent(in) :: fraction_limit, abs_charge_limit
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [10.0_dp, -1.0_dp, -1.0_dp]
    tri_v1(:, 1) = [10.0_dp, 1.0_dp, -1.0_dp]
    tri_v2(:, 1) = [10.0_dp, 0.0_dp, 1.0_dp]
    call init_mesh(soft_mesh, tri_v0, tri_v1, tri_v2)
    soft_mesh%elem_vacuum_sign = 1_i32
    soft_mesh%vacuum_normals = soft_mesh%normals
    call default_app_config(soft_cfg)
    soft_cfg%sim%rng_seed = 995_i32
    soft_cfg%sim%batch_count = 1_i32
    soft_cfg%sim%dt = 1.0_dp
    soft_cfg%sim%max_step = 1_i32
    soft_cfg%sim%use_box = .true.
    soft_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    soft_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    soft_cfg%sim%bc_low(1) = bc_reflect
    soft_cfg%sim%bc_high(1) = bc_reflect
    soft_cfg%sim%multiple_box_events_policy = 'soft_discard'
    soft_cfg%sim%multiple_box_events_soft_discard_count_grace = count_grace
    soft_cfg%sim%multiple_box_events_soft_discard_fraction_limit = fraction_limit
    soft_cfg%sim%multiple_box_events_soft_discard_abs_charge_limit = abs_charge_limit
    soft_cfg%n_particle_species = 1_i32
    soft_cfg%particle_species(1) = species_from_defaults()
    soft_cfg%particle_species(1)%source_mode = 'volume_seed'
    soft_cfg%particle_species(1)%npcls_per_step = particle_count
    soft_cfg%particle_species(1)%q_particle = 2.0_dp
    soft_cfg%particle_species(1)%m_particle = 1.0_dp
    soft_cfg%particle_species(1)%w_particle = 3.0_dp
    soft_cfg%particle_species(1)%pos_low = [0.9_dp, 0.2_dp, 0.2_dp]
    soft_cfg%particle_species(1)%pos_high = soft_cfg%particle_species(1)%pos_low
    soft_cfg%particle_species(1)%drift_velocity = [9.0_dp, 0.0_dp, 0.0_dp]
    soft_cfg%particle_species(1)%temperature_k = 0.0_dp
    call seed_particles_from_config(soft_cfg)
  end subroutine configure_multiple_box_event_soft_discard_fixture
end program test_soft_discard
