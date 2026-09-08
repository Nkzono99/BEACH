!> 実行サマリ・最終CSV・履歴CSVの出力を担当するモジュール。
module bem_output_writer
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_stats
  use bem_app_config_types, only: app_config
  use bem_charge_ledger, only: charge_ledger_type
  use bem_checkpoint_contract, only: begin_checkpoint_publish
  use bem_electrostatic_snapshot, only: electrostatic_diagnostics_type
  use bem_filesystem, only: &
    create_directories, &
    filesystem_empty_path, &
    filesystem_not_directory, &
    filesystem_os_error, &
    filesystem_success
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

  interface
    module subroutine open_history_writer(app, resumed, history_opened, history_unit)
      type(app_config), intent(in) :: app
      logical, intent(in) :: resumed
      logical, intent(out) :: history_opened
      integer, intent(out) :: history_unit
    end subroutine open_history_writer

    module subroutine open_potential_history_writer(app, resumed, potential_history_opened, potential_history_unit)
      type(app_config), intent(in) :: app
      logical, intent(in) :: resumed
      logical, intent(out) :: potential_history_opened
      integer, intent(out) :: potential_history_unit
    end subroutine open_potential_history_writer

    module subroutine open_top_reference_history_writer(app, resumed, history_opened, history_unit)
      type(app_config), intent(in) :: app
      logical, intent(in) :: resumed
      logical, intent(out) :: history_opened
      integer, intent(out) :: history_unit
    end subroutine open_top_reference_history_writer

    module subroutine write_top_reference_history_snapshot( &
      unit_id, batch_idx, simulated_time_s, z_high_m, sample_n, &
      potential_mean_v, potential_std_v, potential_min_v, potential_max_v &
      )
      integer, intent(in) :: unit_id
      integer(i32), intent(in) :: batch_idx, sample_n
      real(dp), intent(in) :: simulated_time_s, z_high_m
      real(dp), intent(in) :: potential_mean_v, potential_std_v, potential_min_v, potential_max_v
    end subroutine write_top_reference_history_snapshot

    module subroutine open_matching_plane_history_writer(app, resumed, history_opened, history_unit)
      type(app_config), intent(in) :: app
      logical, intent(in) :: resumed
      logical, intent(out) :: history_opened
      integer, intent(out) :: history_unit
    end subroutine open_matching_plane_history_writer

    module subroutine write_matching_plane_history_snapshot(unit_id, batch_idx, simulated_time_s, stats)
      integer, intent(in) :: unit_id
      integer(i32), intent(in) :: batch_idx
      real(dp), intent(in) :: simulated_time_s
      type(sim_stats), intent(in) :: stats
    end subroutine write_matching_plane_history_snapshot

    module subroutine write_summary_file( &
      out_dir, mesh, stats, cfg, mpi_world_size, charge_ledger, electrostatic_diagnostics &
      )
      character(len=*), intent(in) :: out_dir
      type(mesh_type), intent(in) :: mesh
      type(sim_stats), intent(in) :: stats
      type(app_config), intent(in) :: cfg
      integer(i32), intent(in), optional :: mpi_world_size
      type(charge_ledger_type), intent(in), optional :: charge_ledger
      type(electrostatic_diagnostics_type), intent(in), optional :: electrostatic_diagnostics
    end subroutine write_summary_file

    module pure real(dp) function soft_discard_fraction(stats) result(fraction)
      type(sim_stats), intent(in) :: stats
    end function soft_discard_fraction

    module integer(i32) function count_dielectric_surfaces(mesh) result(n)
      type(mesh_type), intent(in) :: mesh
    end function count_dielectric_surfaces

    module subroutine write_charge_ledger_file(out_dir, ledger)
      character(len=*), intent(in) :: out_dir
      type(charge_ledger_type), intent(in) :: ledger
    end subroutine write_charge_ledger_file

    module subroutine write_charges_file(out_dir, mesh)
      character(len=*), intent(in) :: out_dir
      type(mesh_type), intent(in) :: mesh
    end subroutine write_charges_file

    module subroutine write_mesh_potential_file(out_dir, mesh, potential_v)
      character(len=*), intent(in) :: out_dir
      type(mesh_type), intent(in) :: mesh
      real(dp), intent(in) :: potential_v(:)
    end subroutine write_mesh_potential_file

    module subroutine write_mesh_file(out_dir, mesh)
      character(len=*), intent(in) :: out_dir
      type(mesh_type), intent(in) :: mesh
    end subroutine write_mesh_file

    module subroutine write_mesh_sources_file(out_dir, mesh, cfg)
      character(len=*), intent(in) :: out_dir
      type(mesh_type), intent(in) :: mesh
      type(app_config), intent(in) :: cfg
    end subroutine write_mesh_sources_file
  end interface

contains

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

end module bem_output_writer
