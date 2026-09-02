!> 実行サマリ・最終CSV・履歴CSVの出力を担当するモジュール。
module bem_output_writer
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_stats, surface_model_insulator, surface_model_conductor, surface_model_dielectric
  use bem_app_config_types, only: app_config
  use bem_charge_ledger, only: charge_ledger_type, finite_charge_sum
  use bem_checkpoint_contract, only: begin_checkpoint_publish, checkpoint_schema_version_current
  use bem_electrostatic_snapshot, only: electrostatic_diagnostics_type
  use bem_field_solver, only: field_solver_fmm_expansion_order, resolve_field_solver_mode, &
                              resolve_field_solver_tree_params
  use bem_external_boundary_contract, only: external_boundary_contract_type, external_boundary_ok, &
                                            external_inflow_none, external_inflow_scalar_barrier, &
                                            external_open_escape, external_open_potential_barrier, &
                                            resolve_external_boundary_contract
  use bem_model_fingerprint, only: model_fingerprint, mesh_fingerprint, species_fingerprint
  use bem_matching_plane_response, only: get_matching_plane_response_content_fingerprint, matching_plane_response_ok
  use bem_physics_config_types, only: field_physics_config, panel_kernel_config, derive_field_panel_config
  use bem_surface_current_model, only: surface_current_model_result_type, evaluate_surface_current_model
  use bem_version, only: beach_build_id, beach_source_commit, beach_version, beach_version_mode
  use bem_filesystem, only: create_directories, filesystem_empty_path, filesystem_not_directory, filesystem_os_error, &
                            filesystem_success
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  public :: open_history_writer
  public :: open_potential_history_writer
  public :: open_top_reference_history_writer
  public :: open_matching_plane_history_writer
  public :: write_top_reference_history_snapshot
  public :: write_matching_plane_history_snapshot
  public :: print_run_summary
  public :: write_result_files
  public :: write_checkpoint_state_files
  public :: ensure_output_dir

contains

  !> 履歴 CSV のオープンとヘッダ初期化を行う。
  !! @param[in] app 出力設定を含むアプリ設定。
  !! @param[in] resumed 再開実行かどうか。
  !! @param[out] history_opened 履歴ファイルを開けた場合に `.true.`。
  !! @param[out] history_unit 履歴CSVの出力ユニット番号（未使用時は `-1`）。
  subroutine open_history_writer(app, resumed, history_opened, history_unit)
    type(app_config), intent(in) :: app
    logical, intent(in) :: resumed
    logical, intent(out) :: history_opened
    integer, intent(out) :: history_unit
    character(len=1024) :: history_path
    integer :: ios
    logical :: history_exists

    history_opened = .false.
    history_unit = -1
    if (.not. app%write_output) return
    if (app%history_stride <= 0) return

    call ensure_output_dir(app%output_dir)

    history_path = trim(app%output_dir)//'/charge_history.csv'
    inquire (file=trim(history_path), exist=history_exists)
    if (resumed) then
      open (newunit=history_unit, file=trim(history_path), status='unknown', position='append', action='write', iostat=ios)
    else
      open (newunit=history_unit, file=trim(history_path), status='replace', action='write', iostat=ios)
    end if
    if (ios /= 0) error stop 'Failed to open charge history file.'

    if (.not. resumed .or. .not. history_exists) then
      ! 再開時は既存ファイルがある場合だけヘッダ追記を避ける。
      write (history_unit, '(a)') 'batch,processed_particles,rel_change,elem_idx,charge_C'
    end if
    history_opened = .true.
  end subroutine open_history_writer

  !> 電位履歴 CSV のオープンとヘッダ初期化を行う。
  !! @param[in] app 出力設定を含むアプリ設定。
  !! @param[in] resumed 再開実行かどうか。
  !! @param[out] potential_history_opened ファイルを開けた場合に `.true.`。
  !! @param[out] potential_history_unit 電位履歴CSVの出力ユニット番号（未使用時は `-1`）。
  subroutine open_potential_history_writer(app, resumed, potential_history_opened, potential_history_unit)
    type(app_config), intent(in) :: app
    logical, intent(in) :: resumed
    logical, intent(out) :: potential_history_opened
    integer, intent(out) :: potential_history_unit
    character(len=1024) :: path
    integer :: ios
    logical :: file_exists

    potential_history_opened = .false.
    potential_history_unit = -1
    if (.not. app%write_output) return
    if (.not. app%write_potential_history) return
    if (app%history_stride <= 0) return

    call ensure_output_dir(app%output_dir)

    path = trim(app%output_dir)//'/potential_history.csv'
    inquire (file=trim(path), exist=file_exists)
    if (resumed) then
      open (newunit=potential_history_unit, file=trim(path), status='unknown', position='append', action='write', iostat=ios)
    else
      open (newunit=potential_history_unit, file=trim(path), status='replace', action='write', iostat=ios)
    end if
    if (ios /= 0) error stop 'Failed to open potential history file.'

    if (.not. resumed .or. .not. file_exists) then
      write (potential_history_unit, '(a)') 'batch,elem_idx,potential_V'
    end if
    potential_history_opened = .true.
  end subroutine open_potential_history_writer

  !> z-high 面平均電位基準の履歴 CSV をオープンし、必要ならヘッダを初期化する。
  !!
  !! 要素電位履歴の companion file として `output.write_potential_history`
  !! および `output.history_stride` に連動する。
  subroutine open_top_reference_history_writer(app, resumed, history_opened, history_unit)
    type(app_config), intent(in) :: app
    logical, intent(in) :: resumed
    logical, intent(out) :: history_opened
    integer, intent(out) :: history_unit
    character(len=1024) :: path
    integer :: ios
    logical :: file_exists

    history_opened = .false.
    history_unit = -1
    if (.not. app%write_output) return
    if (.not. app%write_potential_history) return
    if (app%history_stride <= 0) return
    if (.not. app%sim%use_box) return

    call ensure_output_dir(app%output_dir)

    path = trim(app%output_dir)//'/top_reference_history.csv'
    inquire (file=trim(path), exist=file_exists)
    if (resumed) then
      open (newunit=history_unit, file=trim(path), status='unknown', position='append', action='write', iostat=ios)
    else
      open (newunit=history_unit, file=trim(path), status='replace', action='write', iostat=ios)
    end if
    if (ios /= 0) error stop 'Failed to open top-reference potential history file.'

    if (.not. resumed .or. .not. file_exists) then
      write (history_unit, '(a)') &
        'batch,simulated_time_s,z_high_m,sample_n,potential_mean_V,potential_std_V,'// &
        'potential_min_V,potential_max_V'
    end if
    history_opened = .true.
  end subroutine open_top_reference_history_writer

  !> z-high 面平均電位基準を履歴 CSV へ1行書き出す。
  subroutine write_top_reference_history_snapshot( &
    unit_id, batch_idx, simulated_time_s, z_high_m, sample_n, &
    potential_mean_v, potential_std_v, potential_min_v, potential_max_v &
    )
    integer, intent(in) :: unit_id
    integer(i32), intent(in) :: batch_idx, sample_n
    real(dp), intent(in) :: simulated_time_s, z_high_m
    real(dp), intent(in) :: potential_mean_v, potential_std_v, potential_min_v, potential_max_v

    write (unit_id, '(i0,2(a,es24.16),a,i0,4(a,es24.16))') &
      batch_idx, ',', simulated_time_s, ',', z_high_m, ',', sample_n, &
      ',', potential_mean_v, ',', potential_std_v, ',', potential_min_v, ',', potential_max_v
  end subroutine write_top_reference_history_snapshot

  subroutine open_matching_plane_history_writer(app, resumed, history_opened, history_unit)
    type(app_config), intent(in) :: app
    logical, intent(in) :: resumed
    logical, intent(out) :: history_opened
    integer, intent(out) :: history_unit
    character(len=1024) :: path
    integer :: ios
    logical :: file_exists

    history_opened = .false.
    history_unit = -1
    if (.not. app%write_output) return
    path = trim(app%output_dir)//'/matching_plane_history.csv'
    if (app%history_stride <= 0_i32 .or. &
        trim(lower_ascii(app%surface_current%model)) /= 'matching_plane_quasistatic') then
      if (.not. resumed) call delete_matching_plane_history_if_exists(path)
      return
    end if
    call ensure_output_dir(app%output_dir)
    inquire (file=trim(path), exist=file_exists)
    if (resumed) then
      open (newunit=history_unit, file=trim(path), status='unknown', position='append', action='write', iostat=ios)
    else
      open (newunit=history_unit, file=trim(path), status='replace', action='write', iostat=ios)
    end if
    if (ios /= 0) error stop 'Failed to open matching-plane history file.'
    if (.not. resumed .or. .not. file_exists) then
      write (history_unit, '(a)') &
        'batch,simulated_time_s,D_H_C_m2,phi_H_V,electron_inward_flux_m2_s,ion_inward_flux_m2_s,'// &
        'electron_access_potential_V,ion_access_potential_V,photoelectron_barrier_potential_V,'// &
        'photoelectron_outward_flux_m2_s,'// &
        'photoelectron_mean_normal_energy_eV,electron_outward_flux_m2_s,ion_outward_flux_m2_s,'// &
        'photoelectron_return_flux_m2_s,photoelectron_escape_flux_m2_s,iterations,residual'
    end if
    history_opened = .true.
  end subroutine open_matching_plane_history_writer

  !> fresh runで生成条件を満たさない旧matching historyを残さない。
  subroutine delete_matching_plane_history_if_exists(path)
    character(len=*), intent(in) :: path
    integer :: unit_id, ios
    logical :: file_exists

    inquire (file=trim(path), exist=file_exists)
    if (.not. file_exists) return
    open (newunit=unit_id, file=trim(path), status='old', iostat=ios)
    if (ios /= 0) error stop 'Failed to open stale matching-plane history file for deletion.'
    close (unit_id, status='delete', iostat=ios)
    if (ios /= 0) error stop 'Failed to delete stale matching-plane history file.'
  end subroutine delete_matching_plane_history_if_exists

  subroutine write_matching_plane_history_snapshot(unit_id, batch_idx, simulated_time_s, stats)
    integer, intent(in) :: unit_id
    integer(i32), intent(in) :: batch_idx
    real(dp), intent(in) :: simulated_time_s
    type(sim_stats), intent(in) :: stats

    if (.not. stats%matching_plane_state_valid) return
    write (unit_id, '(i0,14(a,es24.16),a,i0,a,es24.16)') &
      batch_idx, ',', simulated_time_s, ',', stats%matching_plane_displacement_c_m2, &
      ',', stats%matching_plane_phi_v, ',', stats%matching_plane_response(2), &
      ',', stats%matching_plane_response(3), ',', stats%matching_plane_response(4), &
      ',', stats%matching_plane_response(5), ',', stats%matching_plane_response(6), &
      ',', stats%matching_plane_feedback(1), &
      ',', stats%matching_plane_feedback(2), ',', stats%matching_plane_feedback(3), &
      ',', stats%matching_plane_feedback(4), ',', stats%matching_plane_photoelectron_return_flux_m2_s, &
      ',', stats%matching_plane_photoelectron_escape_flux_m2_s, ',', stats%matching_plane_iterations, &
      ',', stats%matching_plane_residual
  end subroutine write_matching_plane_history_snapshot

  !> 実行結果の主要統計を標準出力へ表示する。
  !! @param[in] mesh 実行後のメッシュ情報。
  !! @param[in] stats 実行後の統計値。
  subroutine print_run_summary(mesh, stats)
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    integer(i32) :: dielectric_count

    print '(a,i0)', 'mesh nelem=', mesh%nelem
    print '(a,i0)', 'processed_particles=', stats%processed_particles
    print '(a,i0)', 'absorbed=', stats%absorbed
    print '(a,i0)', 'escaped=', stats%escaped
    print '(a,i0)', 'batches=', stats%batches
    print '(a,i0)', 'escaped_boundary=', stats%escaped_boundary
    print '(a,i0)', 'survived_max_step=', stats%survived_max_step
    print '(a,i0)', 'multiple_box_events_retry_attempted=', stats%multiple_box_events_retry_attempted
    print '(a,i0)', 'multiple_box_events_retry_resolved=', stats%multiple_box_events_retry_resolved
    print '(a,i0)', 'multiple_box_events_soft_discarded=', stats%multiple_box_events_soft_discarded
    print '(a,es12.4)', 'multiple_box_events_soft_discard_fraction=', soft_discard_fraction(stats)
    print '(a,es12.4)', 'multiple_box_events_soft_discarded_abs_charge_C=', &
      stats%multiple_box_events_soft_discarded_abs_charge
    print '(a,es12.4)', 'last_rel_change=', stats%last_rel_change
    print '(a,es12.4)', 'simulated_time_s=', stats%simulated_time
    print '(a,i0)', 'adaptive_nonzero_mode_rejected_trials=', &
      stats%adaptive_nonzero_mode_rejected_trials
    print '(a,es12.4)', 'adaptive_nonzero_mode_last_batch_duration_s=', &
      stats%adaptive_nonzero_mode_last_batch_duration
    print '(a,es12.4)', 'adaptive_nonzero_mode_last_potential_step_V=', &
      stats%adaptive_nonzero_mode_last_potential_step
    print '(a,i0)', 'adaptive_nonzero_mode_omp_threads=', &
      stats%adaptive_nonzero_mode_omp_threads
    print '(a,l1)', 'matching_plane_state_valid=', stats%matching_plane_state_valid
    if (stats%matching_plane_state_valid) then
      print '(a,es12.4)', 'matching_plane_displacement_C_m2=', stats%matching_plane_displacement_c_m2
      print '(a,es12.4)', 'matching_plane_phi_V=', stats%matching_plane_phi_v
      print '(a,6(1x,es12.4))', 'matching_plane_response=', stats%matching_plane_response
      print '(a,4(1x,es12.4))', 'matching_plane_feedback=', stats%matching_plane_feedback
      print '(a,es12.4)', 'matching_plane_photoelectron_return_flux_m2_s=', &
        stats%matching_plane_photoelectron_return_flux_m2_s
      print '(a,es12.4)', 'matching_plane_photoelectron_escape_flux_m2_s=', &
        stats%matching_plane_photoelectron_escape_flux_m2_s
      print '(a,i0)', 'matching_plane_iterations=', stats%matching_plane_iterations
      print '(a,es12.4)', 'matching_plane_residual=', stats%matching_plane_residual
    end if
    print '(a,*(es12.4,1x))', 'mesh charges=', mesh%q_elem
    dielectric_count = count_dielectric_surfaces(mesh)
    if (dielectric_count > 0_i32) then
      print '(a,i0)', 'surface_model_dielectric_elem_count=', dielectric_count
      print '(a)', 'surface_model_note=dielectric surface models are metadata-only in this version.'
    end if
  end subroutine print_run_summary

  !> 解析結果を `summary.txt` / `charges.csv` / `mesh_triangles.csv` などとして保存する。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 書き出し対象のメッシュ。
  !! @param[in] stats 書き出し対象の統計値。
  !! @param[in] cfg 出力設定を含むアプリ設定。
  subroutine write_result_files( &
    out_dir, mesh, stats, cfg, mpi_world_size, mesh_potential_v, charge_ledger, electrostatic_diagnostics &
    )
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    type(app_config), intent(in) :: cfg
    integer(i32), intent(in), optional :: mpi_world_size
    real(dp), intent(in), optional :: mesh_potential_v(:)
    type(charge_ledger_type), intent(in), optional :: charge_ledger
    type(electrostatic_diagnostics_type), intent(in), optional :: electrostatic_diagnostics
    call ensure_output_dir(out_dir)
    call begin_checkpoint_publish(out_dir)
    call write_summary_file( &
      out_dir, mesh, stats, cfg, mpi_world_size=mpi_world_size, charge_ledger=charge_ledger, &
      electrostatic_diagnostics=electrostatic_diagnostics &
      )
    call write_charges_file(out_dir, mesh)
    if (cfg%write_mesh_potential) then
      if (.not. present(mesh_potential_v)) then
        error stop 'write_result_files: mesh_potential_v is required when write_mesh_potential is enabled.'
      end if
      call write_mesh_potential_file(out_dir, mesh, mesh_potential_v)
    end if
    call write_mesh_file(out_dir, mesh)
    call write_mesh_sources_file(out_dir, mesh, cfg)
    if (present(charge_ledger)) call write_charge_ledger_file(out_dir, charge_ledger)
  end subroutine write_result_files

  !> 再開に必要な root-rank 共通状態だけをチェックポイントへ保存する。
  subroutine write_checkpoint_state_files(out_dir, mesh, stats, cfg, mpi_world_size, charge_ledger)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    type(app_config), intent(in) :: cfg
    integer(i32), intent(in), optional :: mpi_world_size
    type(charge_ledger_type), intent(in), optional :: charge_ledger

    call ensure_output_dir(out_dir)
    call begin_checkpoint_publish(out_dir)
    call write_summary_file( &
      out_dir, mesh, stats, cfg, mpi_world_size=mpi_world_size, charge_ledger=charge_ledger &
      )
    call write_charges_file(out_dir, mesh)
    if (present(charge_ledger)) call write_charge_ledger_file(out_dir, charge_ledger)
  end subroutine write_checkpoint_state_files

  !> 出力ディレクトリを作成する。
  !! @param[in] out_dir 作成対象ディレクトリのパス。
  subroutine ensure_output_dir(out_dir)
    character(len=*), intent(in) :: out_dir
    integer :: status

    call create_directories(out_dir, status)
    select case (status)
    case (filesystem_success)
      return
    case (filesystem_empty_path)
      error stop 'Failed to create output directory: output path is empty.'
    case (filesystem_not_directory)
      error stop 'Failed to create output directory: a path component is not an accessible directory.'
    case (filesystem_os_error)
      error stop 'Failed to create output directory: operating-system directory creation failed.'
    case default
      error stop 'Failed to create output directory: unexpected filesystem status.'
    end select
  end subroutine ensure_output_dir

  !> 実行統計を `summary.txt` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh メッシュ情報（要素数を書き出す）。
  !! @param[in] stats 実行統計。
  subroutine write_summary_file( &
    out_dir, mesh, stats, cfg, mpi_world_size, charge_ledger, electrostatic_diagnostics &
    )
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    type(app_config), intent(in) :: cfg
    integer(i32), intent(in), optional :: mpi_world_size
    type(charge_ledger_type), intent(in), optional :: charge_ledger
    type(electrostatic_diagnostics_type), intent(in), optional :: electrostatic_diagnostics
    type(external_boundary_contract_type) :: resolved_boundary
    type(surface_current_model_result_type) :: current_model
    type(field_physics_config) :: field_config
    type(panel_kernel_config) :: panel_config
    character(len=1024) :: summary_path
    character(len=512) :: matching_response_message
    character(len=256) :: boundary_message
    character(len=16) :: resolved_field_solver, matching_response_fingerprint
    integer :: u, ios
    integer(i32) :: world_size, boundary_status, resolved_tree_leaf_max, matching_response_status
    real(dp) :: resolved_tree_theta

    call resolve_external_boundary_contract( &
      cfg%sim%reservoir_potential_model, cfg%sim%open_boundary_model, &
      resolved_boundary, boundary_status, boundary_message &
      )
    if (boundary_status /= external_boundary_ok) then
      error stop 'write_summary_file: invalid local boundary contract: '//trim(boundary_message)
    end if
    call evaluate_surface_current_model(cfg, current_model)
    call derive_field_panel_config(cfg%sim, field_config, panel_config)
    call resolve_field_solver_tree_params( &
      mesh%nelem, cfg%sim, resolved_tree_theta, resolved_tree_leaf_max &
      )
    resolved_field_solver = resolve_field_solver_mode(mesh%nelem, cfg%sim)
    summary_path = trim(out_dir)//'/summary.txt'
    open (newunit=u, file=trim(summary_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open summary file.'
    world_size = 1_i32
    if (present(mpi_world_size)) world_size = max(1_i32, mpi_world_size)
    write (u, '(a,i0)') 'checkpoint_schema_version=', checkpoint_schema_version_current
    write (u, '(a)') 'build_info_schema_version=1'
    write (u, '(a,a)') 'build_version=', beach_version
    write (u, '(a,a)') 'build_version_mode=', beach_version_mode
    write (u, '(a,a)') 'build_source_commit=', beach_source_commit
    write (u, '(a,a)') 'build_id=', beach_build_id
    write (u, '(a,a)') 'model_fingerprint=', model_fingerprint(cfg)
    write (u, '(a,a)') 'mesh_fingerprint=', mesh_fingerprint(mesh)
    write (u, '(a,a)') 'species_fingerprint=', species_fingerprint(cfg)
    write (u, '(a,i0)') 'mesh_nelem=', mesh%nelem
    write (u, '(a,i0)') 'mesh_count=', max(1_i32, maxval(mesh%elem_mesh_id))
    write (u, '(a,i0)') 'mpi_world_size=', world_size
    write (u, '(a,i0)') 'processed_particles=', stats%processed_particles
    write (u, '(a,i0)') 'absorbed=', stats%absorbed
    write (u, '(a,i0)') 'escaped=', stats%escaped
    write (u, '(a,i0)') 'batches=', stats%batches
    write (u, '(a,i0)') 'checkpoint_stride=', cfg%checkpoint_stride
    write (u, '(a,i0)') 'escaped_boundary=', stats%escaped_boundary
    write (u, '(a,i0)') 'survived_max_step=', stats%survived_max_step
    write (u, '(a,i0)') 'multiple_box_events_retry_attempted=', stats%multiple_box_events_retry_attempted
    write (u, '(a,i0)') 'multiple_box_events_retry_resolved=', stats%multiple_box_events_retry_resolved
    write (u, '(a,i0)') 'multiple_box_events_soft_discarded=', stats%multiple_box_events_soft_discarded
    write (u, '(a,es24.16)') 'multiple_box_events_soft_discard_fraction=', soft_discard_fraction(stats)
    write (u, '(a,es24.16)') 'multiple_box_events_soft_discarded_abs_charge_C=', &
      stats%multiple_box_events_soft_discarded_abs_charge
    write (u, '(a,es24.16)') 'last_rel_change=', stats%last_rel_change
    write (u, '(a,es24.16)') 'simulated_time_s=', stats%simulated_time
    write (u, '(a,i0)') 'adaptive_nonzero_mode_rejected_trials=', &
      stats%adaptive_nonzero_mode_rejected_trials
    write (u, '(a,es24.16)') 'adaptive_nonzero_mode_last_batch_duration_s=', &
      stats%adaptive_nonzero_mode_last_batch_duration
    write (u, '(a,es24.16)') 'adaptive_nonzero_mode_last_potential_step_V=', &
      stats%adaptive_nonzero_mode_last_potential_step
    write (u, '(a,i0)') 'adaptive_nonzero_mode_omp_threads=', &
      stats%adaptive_nonzero_mode_omp_threads
    write (u, '(a,l1)') 'matching_plane_state_valid=', stats%matching_plane_state_valid
    write (u, '(a,es24.16)') 'matching_plane_displacement_C_m2=', stats%matching_plane_displacement_c_m2
    write (u, '(a,es24.16)') 'matching_plane_phi_V=', stats%matching_plane_phi_v
    write (u, '(a,es24.16)') 'matching_plane_electron_inward_flux_m2_s=', stats%matching_plane_response(2)
    write (u, '(a,es24.16)') 'matching_plane_ion_inward_flux_m2_s=', stats%matching_plane_response(3)
    write (u, '(a,es24.16)') 'matching_plane_electron_access_potential_V=', stats%matching_plane_response(4)
    write (u, '(a,es24.16)') 'matching_plane_ion_access_potential_V=', stats%matching_plane_response(5)
    write (u, '(a,es24.16)') 'matching_plane_photoelectron_barrier_potential_V=', &
      stats%matching_plane_response(6)
    write (u, '(a,es24.16)') 'matching_plane_photoelectron_outward_flux_m2_s=', &
      stats%matching_plane_feedback(1)
    write (u, '(a,es24.16)') 'matching_plane_photoelectron_mean_normal_energy_eV=', &
      stats%matching_plane_feedback(2)
    write (u, '(a,es24.16)') 'matching_plane_electron_outward_flux_m2_s=', &
      stats%matching_plane_feedback(3)
    write (u, '(a,es24.16)') 'matching_plane_ion_outward_flux_m2_s=', &
      stats%matching_plane_feedback(4)
    write (u, '(a,es24.16)') 'matching_plane_photoelectron_return_flux_m2_s=', &
      stats%matching_plane_photoelectron_return_flux_m2_s
    write (u, '(a,es24.16)') 'matching_plane_photoelectron_escape_flux_m2_s=', &
      stats%matching_plane_photoelectron_escape_flux_m2_s
    write (u, '(a,i0)') 'matching_plane_iterations=', stats%matching_plane_iterations
    write (u, '(a,es24.16)') 'matching_plane_residual=', stats%matching_plane_residual
    write (u, '(a)') 'particle_time_centering=same_time_midpoint_boris'
    write (u, '(a,a)') 'field_backend=', trim(field_config%backend)
    write (u, '(a,a)') 'field_normalization=', trim(field_config%normalization)
    write (u, '(a)') 'field_source_model=triangle_p0'
    write (u, '(a,a)') 'field_kernel_id=', trim(panel_config%kernel_id)
    write (u, '(a)') 'field_reconstruction_schema_version=2'
    write (u, '(a,a)') 'field_reconstruction_resolved_field_solver=', trim(resolved_field_solver)
    write (u, '(a,i0)') 'field_reconstruction_fmm_expansion_order=', field_solver_fmm_expansion_order
    write (u, '(a,a)') 'field_reconstruction_field_bc_mode=', trim(cfg%sim%field_bc_mode)
    write (u, '(a,es24.16)') 'field_reconstruction_tree_theta=', resolved_tree_theta
    write (u, '(a,i0)') 'field_reconstruction_tree_leaf_max=', resolved_tree_leaf_max
    write (u, '(a,3(1x,es24.16))') 'field_reconstruction_e0_V_m=', cfg%sim%e0
    write (u, '(a,l1)') 'field_reconstruction_use_box=', cfg%sim%use_box
    write (u, '(a,3(1x,es24.16))') 'field_reconstruction_box_min_m=', cfg%sim%box_min
    write (u, '(a,3(1x,es24.16))') 'field_reconstruction_box_max_m=', cfg%sim%box_max
    write (u, '(a,3(1x,i0))') 'field_reconstruction_boundary_low=', cfg%sim%bc_low
    write (u, '(a,3(1x,i0))') 'field_reconstruction_boundary_high=', cfg%sim%bc_high
    write (u, '(a,i0)') 'field_reconstruction_periodic_image_layers=', cfg%sim%field_periodic_image_layers
    write (u, '(a,a)') 'field_reconstruction_periodic_far_correction=', &
      trim(cfg%sim%field_periodic_far_correction)
    write (u, '(a,a)') 'field_reconstruction_periodic_nonzero_mode_backend=', &
      trim(cfg%periodic2%nonzero_mode_backend)
    write (u, '(a,a)') 'field_reconstruction_periodic_zero_mode_policy=', &
      trim(cfg%periodic2%zero_mode_policy)
    write (u, '(a,a)') 'field_reconstruction_periodic_lower_boundary_model=', &
      trim(cfg%periodic2%lower_boundary_model)
    write (u, '(a,i0)') 'field_reconstruction_periodic_reference_mode_layers=', &
      cfg%periodic2%reference_mode_layers
    write (u, '(a,i0)') 'field_reconstruction_periodic_panel_quadrature_order=', &
      cfg%periodic2%panel_quadrature_order
    write (u, '(a,es24.16)') 'field_reconstruction_periodic_ewald_alpha=', cfg%sim%field_periodic_ewald_alpha
    write (u, '(a,i0)') 'field_reconstruction_periodic_ewald_layers=', cfg%sim%field_periodic_ewald_layers
    write (u, '(a,a)') 'field_reconstruction_periodic_cache_dir=', trim(cfg%sim%field_periodic_cache_dir)
    write (u, '(a,es24.16)') 'field_reconstruction_periodic_generation_tolerance=', &
      cfg%sim%field_periodic_generation_tolerance
    write (u, '(a,a)') 'periodic2_nonzero_mode_backend=', trim(cfg%periodic2%nonzero_mode_backend)
    write (u, '(a,a)') 'periodic2_zero_mode_policy=', trim(cfg%periodic2%zero_mode_policy)
    write (u, '(a,a)') 'periodic2_lower_boundary_model=', trim(cfg%periodic2%lower_boundary_model)
    write (u, '(a,es24.16)') 'periodic2_max_nonzero_mode_potential_step_V=', &
      cfg%periodic2%max_nonzero_mode_potential_step
    write (u, '(a,a)') 'periodic2_cache_dir=', trim(cfg%sim%field_periodic_cache_dir)
    write (u, '(a,es24.16)') 'periodic2_generation_tolerance=', &
      cfg%sim%field_periodic_generation_tolerance
    write (u, '(a,a)') 'reservoir_inflow_map=', trim(external_inflow_map_name(resolved_boundary%inflow_map))
    write (u, '(a,a)') 'particle_ordinary_open_model=', &
      trim(external_open_model_name(resolved_boundary%ordinary_open_model))
    write (u, '(a,a)') 'surface_current_model=', trim(current_model%model)
    if (current_model%active) then
      write (u, '(a,a)') 'surface_current_model_kinetic_contract=', trim(current_model%kinetic_contract)
      if (trim(current_model%model) == 'matching_plane_quasistatic') then
        write (u, '(a,a)') 'surface_current_model_response_backend=', &
          trim(lower_ascii(cfg%surface_current%response_backend))
        write (u, '(a,l1)') 'surface_current_model_implicit_zero_mode=', &
          cfg%surface_current%implicit_zero_mode
        select case (trim(lower_ascii(cfg%surface_current%response_backend)))
        case ('table')
          call get_matching_plane_response_content_fingerprint( &
            trim(cfg%surface_current%response_table_path), matching_response_fingerprint, &
            matching_response_status, matching_response_message &
            )
          if (matching_response_status /= matching_plane_response_ok) then
            error stop 'write_summary_file: matching-plane response fingerprint failed: '// &
              trim(matching_response_message)
          end if
          write (u, '(a,a)') 'surface_current_model_response_table_path=', &
            trim(cfg%surface_current%response_table_path)
          write (u, '(a,a)') 'surface_current_model_response_content_fingerprint=', &
            matching_response_fingerprint
        case ('zhao_online')
          write (u, '(a)') 'surface_current_model_response_contract=matching_plane_zhao_online_v1'
          write (u, '(a,a)') 'surface_current_model_zhao_branch=', &
            trim(lower_ascii(cfg%surface_current%zhao_branch))
          write (u, '(a,a)') 'surface_current_model_zhao_root_selection=', &
            trim(lower_ascii(cfg%surface_current%zhao_root_selection))
          write (u, '(a)') &
            'surface_current_model_outer_solver=charge_driven_finite_h_sagdeev'
          write (u, '(a)') &
            'surface_current_model_photoelectron_closure=moment_matched_half_maxwellian'
          write (u, '(a)') &
            'surface_current_model_ambient_outward_feedback=transparent'
          if (trim(lower_ascii(cfg%surface_current%zhao_root_selection)) == 'continuation') then
            write (u, '(a)') &
              'surface_current_model_outer_solver_state=accepted_endpoint_continuation_v1'
          else
            write (u, '(a)') 'surface_current_model_outer_solver_state=stateless'
          end if
        case default
          error stop 'write_summary_file: unknown matching-plane response backend.'
        end select
        write (u, '(a,es24.16)') 'surface_current_model_matching_plane_z_m=', cfg%sim%box_max(3)
        write (u, '(a,a)') 'surface_current_model_electron_species=', &
          trim(cfg%surface_current%electron_species)
        write (u, '(a,a)') 'surface_current_model_ion_species=', trim(cfg%surface_current%ion_species)
        write (u, '(a,a)') 'surface_current_model_photoelectron_species=', &
          trim(cfg%surface_current%photoelectron_species)
        write (u, '(a,es24.16)') 'surface_current_model_coupling_rtol=', cfg%surface_current%coupling_rtol
        write (u, '(a,4(1x,es24.16))') 'surface_current_model_coupling_atol=', &
          cfg%surface_current%coupling_atol
        write (u, '(a,i0)') 'surface_current_model_coupling_max_iterations=', &
          cfg%surface_current%coupling_max_iterations
        write (u, '(a,es24.16)') 'surface_current_model_coupling_relaxation=', &
          cfg%surface_current%coupling_relaxation
        write (u, '(a)') 'surface_current_model_dynamic_state_source=accepted_batch_fixed_point'
      else
        write (u, '(a,a)') 'surface_current_model_zhao_branch=', current_model%zhao_branch
        write (u, '(a,l1)') 'surface_current_model_photoelectron_active=', current_model%photoelectron_active
        write (u, '(a,es24.16)') 'surface_current_model_reference_area_m2=', current_model%reference_area_m2
        write (u, '(a,es24.16)') 'surface_current_model_phi0_V=', current_model%phi0_v
        write (u, '(a,es24.16)') 'surface_current_model_phi_m_V=', current_model%phi_m_v
        write (u, '(a,es24.16)') 'surface_current_model_ambient_electron_density_m3=', &
          current_model%ambient_electron_density_m3
        write (u, '(a,es24.16)') 'surface_current_model_electron_current_density_A_m2=', &
          current_model%electron_current_density_a_m2
        write (u, '(a,es24.16)') 'surface_current_model_ion_current_density_A_m2=', &
          current_model%ion_current_density_a_m2
        write (u, '(a,es24.16)') 'surface_current_model_pe_emission_current_density_A_m2=', &
          current_model%photoelectron_emission_current_density_a_m2
        write (u, '(a,es24.16)') 'surface_current_model_pe_escape_current_density_A_m2=', &
          current_model%photoelectron_escape_current_density_a_m2
        write (u, '(a,es24.16)') 'surface_current_model_pe_return_current_density_A_m2=', &
          current_model%photoelectron_return_current_density_a_m2
        write (u, '(a,es24.16)') 'surface_current_model_net_current_density_A_m2=', &
          current_model%net_current_density_a_m2
        write (u, '(a)') 'surface_current_model_current_budget_contract=surface_targets_plus_external_escape'
        write (u, '(a,es24.16)') 'surface_current_model_pe_budget_residual_current_density_A_m2=', &
          current_model%photoelectron_budget_residual_current_density_a_m2
        write (u, '(a,es24.16)') 'surface_current_model_surface_budget_residual_current_density_A_m2=', &
          current_model%surface_budget_residual_current_density_a_m2
        if (current_model%photoelectron_active) then
          write (u, '(a,es24.16)') 'surface_current_model_pe_escape_particle_current_A=', &
            current_model%escaped_particle_current_a(current_model%photoelectron_species_idx)
        else
          write (u, '(a,es24.16)') 'surface_current_model_pe_escape_particle_current_A=', 0.0_dp
        end if
        write (u, '(a,es24.16)') 'surface_current_model_electron_inflow_reservoir_potential_V=', &
          current_model%inflow_reservoir_potential_v(current_model%electron_species_idx)
        write (u, '(a,es24.16)') 'surface_current_model_electron_inflow_access_potential_V=', &
          current_model%inflow_access_potential_v(current_model%electron_species_idx)
        write (u, '(a,i0)') 'surface_current_model_electron_inflow_face=', &
          current_model%inflow_kinetic_face(current_model%electron_species_idx)
        if (current_model%photoelectron_active) then
          write (u, '(a,es24.16)') 'surface_current_model_pe_outflow_barrier_potential_V=', &
            current_model%outflow_barrier_potential_v(current_model%photoelectron_species_idx)
          write (u, '(a,i0)') 'surface_current_model_pe_outflow_barrier_face=', &
            current_model%outflow_barrier_face(current_model%photoelectron_species_idx)
        else
          write (u, '(a,es24.16)') 'surface_current_model_pe_outflow_barrier_potential_V=', 0.0_dp
          write (u, '(a,i0)') 'surface_current_model_pe_outflow_barrier_face=', 0_i32
        end if
      end if
    end if
    if (present(electrostatic_diagnostics)) then
      write (u, '(a,l1)') 'top_reference_available=', electrostatic_diagnostics%top_reference_available
      if (electrostatic_diagnostics%top_reference_available) then
        write (u, '(a)') 'top_reference_definition=box_z_high_plane_mean'
        write (u, '(a,i0)') 'top_reference_last_batch=', electrostatic_diagnostics%top_reference_last_batch
        write (u, '(a,es24.16)') 'top_reference_simulated_time_s=', electrostatic_diagnostics%top_reference_simulated_time
        write (u, '(a,es24.16)') 'top_reference_z_high_m=', electrostatic_diagnostics%top_reference_z_high
        write (u, '(a,i0)') 'top_reference_sample_n=', electrostatic_diagnostics%top_reference_sample_n
        write (u, '(a,es24.16)') 'top_reference_potential_mean_V=', electrostatic_diagnostics%top_reference_potential_mean
        write (u, '(a,es24.16)') 'top_reference_potential_std_V=', electrostatic_diagnostics%top_reference_potential_std
        write (u, '(a,es24.16)') 'top_reference_potential_min_V=', electrostatic_diagnostics%top_reference_potential_min
        write (u, '(a,es24.16)') 'top_reference_potential_max_V=', electrostatic_diagnostics%top_reference_potential_max
      end if
      write (u, '(a,l1)') 'electrostatic_split_periodic_active=', electrostatic_diagnostics%split_periodic_active
      write (u, '(a,a)') 'electrostatic_status=', trim(electrostatic_diagnostics%status)
      write (u, '(a,es24.16)') 'gauss_residual_C=', electrostatic_diagnostics%gauss_residual
      write (u, '(a,l1)') 'periodic2_cache_hit=', electrostatic_diagnostics%periodic_cache_hit
      write (u, '(a,i0)') 'periodic2_operator_build_count=', electrostatic_diagnostics%periodic_operator_build_count
      write (u, '(a,a)') 'periodic2_cache_fingerprint=', trim(electrostatic_diagnostics%periodic_cache_fingerprint)
      write (u, '(a,a)') 'periodic2_cache_path=', trim(electrostatic_diagnostics%periodic_cache_path)
    end if
    if (present(charge_ledger)) then
      write (u, '(a,i0)') 'charge_ledger_nspecies=', charge_ledger%nspecies
      write (u, '(a,i0)') 'charge_ledger_batch_count=', charge_ledger%batch_count
      write (u, '(a,es24.16)') 'charge_ledger_surface_charge_before_C=', charge_ledger%surface_charge_before
      write (u, '(a,es24.16)') 'charge_ledger_surface_charge_after_C=', charge_ledger%surface_charge_after
      write (u, '(a,es24.16)') 'charge_ledger_local_flight_charge_before_C=', &
        charge_ledger%local_flight_charge_before
      write (u, '(a,es24.16)') 'charge_ledger_local_flight_charge_after_C=', &
        charge_ledger%local_flight_charge_after
      write (u, '(a,es24.16)') 'charge_ledger_unresolved_stock_before_C=', &
        charge_ledger%unresolved_stock_before
      write (u, '(a,es24.16)') 'charge_ledger_unresolved_stock_after_C=', &
        charge_ledger%unresolved_stock_after
      write (u, '(a,es24.16)') 'charge_ledger_residual_C=', charge_ledger%residual()
      write (u, '(a,es24.16)') 'charge_ledger_discarded_unresolved_abs_C=', &
        charge_ledger%discarded_unresolved_abs()
      write (u, '(a,es24.16)') 'charge_ledger_neutral_return_correction_C=', &
        finite_charge_sum(charge_ledger%neutral_return_correction, 'summary neutral-return correction')
      write (u, '(a,es24.16)') 'charge_ledger_fixed_current_correction_C=', &
        finite_charge_sum(charge_ledger%fixed_current_correction, 'summary fixed-current correction')
      write (u, '(a,es24.16)') 'charge_ledger_fixed_absorbed_applied_charge_C=', &
        finite_charge_sum(charge_ledger%fixed_absorbed_target_charge, 'summary fixed absorbed applied charge')
      write (u, '(a,es24.16)') 'charge_ledger_fixed_emission_applied_charge_C=', &
        finite_charge_sum(charge_ledger%fixed_emission_target_charge, 'summary fixed emission applied charge')
      write (u, '(a,es24.16)') 'charge_ledger_raw_escape_charge_C=', &
        finite_charge_sum(charge_ledger%escaped_to_infinity, 'summary raw escape charge')
      write (u, '(a,es24.16)') 'charge_ledger_fixed_escape_target_charge_C=', &
        finite_charge_sum(charge_ledger%fixed_escape_target_charge, 'summary fixed escape target charge')
      write (u, '(a,es24.16)') 'charge_ledger_fixed_escape_applied_charge_C=', &
        finite_charge_sum(charge_ledger%fixed_escape_target_charge, 'summary fixed escape applied charge')
      write (u, '(a,es24.16)') 'charge_ledger_fixed_escape_correction_C=', &
        finite_charge_sum(charge_ledger%fixed_escape_correction, 'summary fixed escape correction')
      write (u, '(a,es24.16)') 'charge_ledger_fixed_applied_surface_net_charge_C=', &
        finite_charge_sum( &
        [charge_ledger%fixed_absorbed_target_charge, charge_ledger%fixed_emission_target_charge], &
        'summary fixed applied surface net charge' &
        )
      if (current_model%active .and. current_model%photoelectron_active) then
        write (u, '(a,es24.16)') 'charge_ledger_fixed_pe_continuity_residual_C=', &
          finite_charge_sum( &
          [ &
          charge_ledger%fixed_absorbed_target_charge(current_model%photoelectron_species_idx), &
          charge_ledger%fixed_emission_target_charge(current_model%photoelectron_species_idx), &
          charge_ledger%fixed_escape_target_charge(current_model%photoelectron_species_idx) &
          ], &
          'summary fixed photoelectron continuity residual' &
          )
      end if
    end if
    if (count_dielectric_surfaces(mesh) > 0_i32) then
      write (u, '(a,i0)') 'surface_model_dielectric_elem_count=', count_dielectric_surfaces(mesh)
      write (u, '(a)') 'surface_model_note=metadata_only_dielectric_present'
    end if
    close (u)
  end subroutine write_summary_file

  !> 解決済みの流入モデルを summary 用の安定した語彙へ変換する。
  function external_inflow_map_name(inflow_map) result(name)
    integer(i32), intent(in) :: inflow_map
    character(len=32) :: name

    select case (inflow_map)
    case (external_inflow_none)
      name = 'source_vdf'
    case (external_inflow_scalar_barrier)
      name = 'infinity_barrier'
    case default
      error stop 'write_summary_file: unknown resolved external inflow map.'
    end select
  end function external_inflow_map_name

  !> 解決済みの通常 open 境界モデルを summary 用の安定した語彙へ変換する。
  function external_open_model_name(open_model) result(name)
    integer(i32), intent(in) :: open_model
    character(len=32) :: name

    select case (open_model)
    case (external_open_escape)
      name = 'escape'
    case (external_open_potential_barrier)
      name = 'potential_barrier'
    case default
      error stop 'write_summary_file: unknown resolved ordinary open model.'
    end select
  end function external_open_model_name

  !> species 別の signed charge flux と粒子数を `charge_ledger.csv` に保存する。
  subroutine write_charge_ledger_file(out_dir, ledger)
    character(len=*), intent(in) :: out_dir
    type(charge_ledger_type), intent(in) :: ledger
    character(len=1024) :: path
    integer :: u, ios, species_idx

    if (ledger%nspecies < 1_i32 .or. .not. allocated(ledger%injected_from_remote)) then
      error stop 'write_charge_ledger_file requires an initialized ledger.'
    end if
    path = trim(out_dir)//'/charge_ledger.csv'
    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open charge_ledger.csv.'
    write (u, '(a)') &
      'batch,species_idx,injected_from_remote_C,emitted_from_surface_C,absorbed_on_surface_C,'// &
      'escaped_to_infinity_C,discarded_unresolved_C,'// &
      'neutral_return_correction_C,neutral_return_weight_scale,neutral_return_unresolved_fraction,'// &
      'fixed_absorbed_target_charge_C,fixed_absorbed_weight_scale,'// &
      'fixed_emission_target_charge_C,fixed_emission_weight_scale,fixed_current_correction_C,'// &
      'fixed_absorbed_applied_charge_C,fixed_emission_applied_charge_C,'// &
      'fixed_escape_target_charge_C,fixed_escape_applied_charge_C,fixed_escape_correction_C,'// &
      'injected_count,emitted_count,absorbed_count,escaped_count,discarded_unresolved_count'
    ! Keep the three applied-charge columns as file-format aliases for existing
    ! readers; the ledger stores only the identical target values.
    do species_idx = 1, ledger%nspecies
      write (u, '(i0,a,i0,18(a,es24.16),5(a,i0))') &
        ledger%batch_count, ',', species_idx, &
        ',', ledger%injected_from_remote(species_idx), &
        ',', ledger%emitted_from_surface(species_idx), &
        ',', ledger%absorbed_on_surface(species_idx), &
        ',', ledger%escaped_to_infinity(species_idx), &
        ',', ledger%discarded_unresolved(species_idx), &
        ',', ledger%neutral_return_correction(species_idx), &
        ',', ledger%neutral_return_weight_scale(species_idx), &
        ',', ledger%neutral_return_unresolved_fraction(species_idx), &
        ',', ledger%fixed_absorbed_target_charge(species_idx), &
        ',', ledger%fixed_absorbed_weight_scale(species_idx), &
        ',', ledger%fixed_emission_target_charge(species_idx), &
        ',', ledger%fixed_emission_weight_scale(species_idx), &
        ',', ledger%fixed_current_correction(species_idx), &
        ',', ledger%fixed_absorbed_target_charge(species_idx), &
        ',', ledger%fixed_emission_target_charge(species_idx), &
        ',', ledger%fixed_escape_target_charge(species_idx), &
        ',', ledger%fixed_escape_target_charge(species_idx), &
        ',', ledger%fixed_escape_correction(species_idx), &
        ',', ledger%injected_count(species_idx), &
        ',', ledger%emitted_count(species_idx), &
        ',', ledger%absorbed_count(species_idx), &
        ',', ledger%escaped_count(species_idx), &
        ',', ledger%discarded_unresolved_count(species_idx)
    end do
    close (u)
  end subroutine write_charge_ledger_file

  !> 要素電荷を `charges.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素電荷を含むメッシュ情報。
  subroutine write_charges_file(out_dir, mesh)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    character(len=1024) :: charges_path
    integer :: u, ios, i

    charges_path = trim(out_dir)//'/charges.csv'
    open (newunit=u, file=trim(charges_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open charges file.'
    write (u, '(a)') 'elem_idx,charge_C'
    do i = 1, mesh%nelem
      write (u, '(i0,a,es24.16)') i, ',', mesh%q_elem(i)
    end do
    close (u)
  end subroutine write_charges_file

  !> 事前計算済み電位を `mesh_potential.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素情報（要素数の検証用）。
  !! @param[in] potential_v 各要素重心での電位 [V]。
  subroutine write_mesh_potential_file(out_dir, mesh, potential_v)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: potential_v(:)
    character(len=1024) :: potential_path
    integer :: u, ios, i

    if (size(potential_v) /= mesh%nelem) error stop 'precomputed mesh potential size mismatch.'

    potential_path = trim(out_dir)//'/mesh_potential.csv'
    open (newunit=u, file=trim(potential_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open mesh_potential.csv.'
    write (u, '(a)') 'elem_idx,potential_V'
    do i = 1, mesh%nelem
      write (u, '(i0,a,es24.16)') i, ',', potential_v(i)
    end do
    close (u)
  end subroutine write_mesh_potential_file

  !> 三角形メッシュを `mesh_triangles.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 頂点座標と要素電荷を含むメッシュ情報。
  subroutine write_mesh_file(out_dir, mesh)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    character(len=1024) :: mesh_path
    integer :: u, ios, i

    mesh_path = trim(out_dir)//'/mesh_triangles.csv'
    open (newunit=u, file=trim(mesh_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open mesh file.'
    write (u, '(a)') 'elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id'
    do i = 1, mesh%nelem
      write (u, '(i0,10(a,es24.16),a,i0)') i, ',', mesh%v0(1, i), ',', mesh%v0(2, i), ',', mesh%v0(3, i), &
        ',', mesh%v1(1, i), ',', mesh%v1(2, i), ',', mesh%v1(3, i), &
        ',', mesh%v2(1, i), ',', mesh%v2(2, i), ',', mesh%v2(3, i), ',', &
        mesh%q_elem(i), ',', mesh%elem_mesh_id(i)
    end do
    close (u)
  end subroutine write_mesh_file

  !> メッシュ識別情報を `mesh_sources.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素ごとの `mesh_id` を含むメッシュ情報。
  !! @param[in] cfg 元の入力設定。
  subroutine write_mesh_sources_file(out_dir, mesh, cfg)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(app_config), intent(in) :: cfg
    character(len=1024) :: path
    character(len=16) :: mode_key, source_kind, template_kind
    character(len=16) :: surface_model
    real(dp) :: epsilon_r
    logical :: has_obj
    integer :: u, ios, i, mesh_id
    integer(i32) :: elem_count

    path = trim(out_dir)//'/mesh_sources.csv'
    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open mesh_sources.csv.'
    write (u, '(a)') 'mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count'

    mode_key = trim(lower_ascii(cfg%mesh_mode))
    if (mode_key == 'obj') then
      source_kind = 'obj'
    else if (mode_key == 'template') then
      source_kind = 'template'
    else
      inquire (file=trim(cfg%obj_path), exist=has_obj)
      if (has_obj) then
        source_kind = 'obj'
      else
        source_kind = 'template'
      end if
    end if

    if (source_kind == 'obj') then
      surface_model = mesh_surface_model_name(mesh, 1)
      epsilon_r = mesh_epsilon_r(mesh, 1)
      write (u, '(i0,a,a,a,a,a,a,a,es16.8,a,i0)') &
        1, ',', 'obj', ',', 'obj', ',', trim(surface_model), ',', epsilon_r, ',', mesh%nelem
      close (u)
      return
    end if

    mesh_id = 0
    do i = 1, size(cfg%templates)
      if (.not. cfg%templates(i)%enabled) cycle
      mesh_id = mesh_id + 1
      template_kind = trim(lower_ascii(cfg%templates(i)%kind))
      surface_model = mesh_surface_model_name(mesh, mesh_id)
      epsilon_r = mesh_epsilon_r(mesh, mesh_id)
      elem_count = int(count(mesh%elem_mesh_id == mesh_id), kind=i32)
      write (u, '(i0,a,a,a,a,a,a,a,es16.8,a,i0)') &
        mesh_id, ',', 'template', ',', trim(template_kind), ',', trim(surface_model), ',', epsilon_r, ',', elem_count
    end do
    close (u)
  end subroutine write_mesh_sources_file

  !> mesh_id に対応する表面モデル名を返す。
  !! 複数モデルが混在している場合は最初の要素のモデル名を代表値として返す。
  function mesh_surface_model_name(mesh, mesh_id) result(name)
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: mesh_id
    character(len=16) :: name
    integer :: i

    name = 'insulator'
    if (.not. allocated(mesh%elem_surface_model)) return
    do i = 1, mesh%nelem
      if (mesh%elem_mesh_id(i) /= mesh_id) cycle
      select case (mesh%elem_surface_model(i))
      case (surface_model_insulator)
        name = 'insulator'
      case (surface_model_conductor)
        name = 'conductor'
      case (surface_model_dielectric)
        name = 'dielectric'
      case default
        name = 'unknown'
      end select
      return
    end do
  end function mesh_surface_model_name

  !> mesh_id に対応する相対誘電率を返す。
  function mesh_epsilon_r(mesh, mesh_id) result(epsilon_r)
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: mesh_id
    real(dp) :: epsilon_r
    integer :: i

    epsilon_r = 1.0d0
    if (.not. allocated(mesh%elem_epsilon_r)) return
    do i = 1, mesh%nelem
      if (mesh%elem_mesh_id(i) /= mesh_id) cycle
      epsilon_r = mesh%elem_epsilon_r(i)
      return
    end do
  end function mesh_epsilon_r

  !> 処理済みmacro particleに対する累積soft-discard率を返す。
  pure real(dp) function soft_discard_fraction(stats) result(fraction)
    type(sim_stats), intent(in) :: stats

    fraction = 0.0_dp
    if (stats%processed_particles <= 0) return
    fraction = real(stats%multiple_box_events_soft_discarded, dp)/real(stats%processed_particles, dp)
  end function soft_discard_fraction

  !> 現行物理モデルで未分岐の dielectric 要素数を数える。
  integer(i32) function count_dielectric_surfaces(mesh) result(n)
    type(mesh_type), intent(in) :: mesh

    n = 0_i32
    if (.not. allocated(mesh%elem_surface_model)) return
    n = int(count(mesh%elem_surface_model == surface_model_dielectric), kind=i32)
  end function count_dielectric_surfaces

end module bem_output_writer
