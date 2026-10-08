!> 履歴ファイルの生成・再開とスナップショットの出力。
submodule(bem_output_writer) bem_output_writer_history
  use bem_string_utils, only: lower_ascii
  implicit none
contains

  !> 履歴 CSV のオープンとヘッダ初期化を行う。
  !! @param[in] app 出力設定を含むアプリ設定。
  !! @param[in] resumed 再開実行かどうか。
  !! @param[out] history_opened 履歴ファイルを開けた場合に `.true.`。
  !! @param[out] history_unit 履歴CSVの出力ユニット番号（未使用時は `-1`）。
  module procedure open_history_writer
  character(len=1024) :: history_path

  history_opened = .false.
  history_unit = -1
  if (.not. app%write_output) return
  if (app%history_stride <= 0) return

  call ensure_output_dir(app%output_dir)

  history_path = trim(app%output_dir)//'/charge_history.csv'
  call open_history_file( &
    history_path, resumed, &
    'batch,processed_particles,rel_change,elem_idx,charge_C', &
    history_unit, 'Failed to open charge history file.' &
    )
  history_opened = .true.
  end procedure open_history_writer

  !> 電位履歴 CSV のオープンとヘッダ初期化を行う。
  !! @param[in] app 出力設定を含むアプリ設定。
  !! @param[in] resumed 再開実行かどうか。
  !! @param[out] potential_history_opened ファイルを開けた場合に `.true.`。
  !! @param[out] potential_history_unit 電位履歴CSVの出力ユニット番号（未使用時は `-1`）。
  module procedure open_potential_history_writer
  character(len=1024) :: path

  potential_history_opened = .false.
  potential_history_unit = -1
  if (.not. app%write_output) return
  if (.not. app%write_potential_history) return
  if (app%history_stride <= 0) return

  call ensure_output_dir(app%output_dir)

  path = trim(app%output_dir)//'/potential_history.csv'
  call open_history_file( &
    path, resumed, &
    'batch,elem_idx,potential_V', &
    potential_history_unit, 'Failed to open potential history file.' &
    )
  potential_history_opened = .true.
  end procedure open_potential_history_writer

  !> z-high 面平均電位基準の履歴 CSV をオープンし、必要ならヘッダを初期化する。
  !!
  !! 要素電位履歴の companion file として `output.write_potential_history`
  !! および `output.history_stride` に連動する。
  module procedure open_top_reference_history_writer
  character(len=1024) :: path

  history_opened = .false.
  history_unit = -1
  if (.not. app%write_output) return
  if (.not. app%write_potential_history) return
  if (app%history_stride <= 0) return
  if (.not. app%sim%use_box) return

  call ensure_output_dir(app%output_dir)

  path = trim(app%output_dir)//'/top_reference_history.csv'
  call open_history_file( &
    path, resumed, &
    'batch,simulated_time_s,z_high_m,sample_n,potential_mean_V,potential_std_V,'// &
    'potential_min_V,potential_max_V', &
    history_unit, 'Failed to open top-reference potential history file.' &
    )
  history_opened = .true.
  end procedure open_top_reference_history_writer

  !> z-high 面平均電位基準を履歴 CSV へ1行書き出す。
  module procedure write_top_reference_history_snapshot

  write (unit_id, '(i0,2(a,es24.16),a,i0,4(a,es24.16))') &
    batch_idx, ',', simulated_time_s, ',', z_high_m, ',', sample_n, &
    ',', potential_mean_v, ',', potential_std_v, ',', potential_min_v, ',', potential_max_v
  end procedure write_top_reference_history_snapshot

  module procedure open_matching_plane_history_writer
  character(len=1024) :: path

  history_opened = .false.
  history_unit = -1
  if (.not. app%write_output) return
  path = trim(app%output_dir)//'/matching_plane_history.csv'
  if (app%history_stride <= 0_i32 .or. &
      trim(lower_ascii(app%surface_current%model)) /= 'zhao_stationary' .or. &
      app%surface_current%outflow_refresh_batches <= 0_i32) then
    if (.not. resumed) call delete_history_if_exists(path)
    return
  end if
  call ensure_output_dir(app%output_dir)
  call open_history_file( &
    path, resumed, &
    'batch,simulated_time_s,D_H_C_m2,phi_H_V,electron_inward_flux_m2_s,ion_inward_flux_m2_s,'// &
    'electron_access_potential_V,ion_access_potential_V,photoelectron_barrier_potential_V,'// &
    'photoelectron_outward_flux_m2_s,'// &
    'photoelectron_mean_normal_energy_eV,electron_outward_flux_m2_s,ion_outward_flux_m2_s,'// &
    'photoelectron_return_flux_m2_s,photoelectron_escape_flux_m2_s,iterations,residual', &
    history_unit, 'Failed to open matching-plane history file.' &
    )
  history_opened = .true.
  end procedure open_matching_plane_history_writer

  !> fresh runで生成条件を満たさない旧条件付きhistoryを残さない。
  subroutine delete_history_if_exists(path)
    character(len=*), intent(in) :: path
    integer :: unit_id, ios
    logical :: file_exists

    inquire (file=trim(path), exist=file_exists)
    if (.not. file_exists) return
    open (newunit=unit_id, file=trim(path), status='old', iostat=ios)
    if (ios /= 0) error stop 'Failed to open stale history file for deletion: '//trim(path)
    close (unit_id, status='delete', iostat=ios)
    if (ios /= 0) error stop 'Failed to delete stale history file: '//trim(path)
  end subroutine delete_history_if_exists

  module procedure write_matching_plane_history_snapshot

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
  end procedure write_matching_plane_history_snapshot

  !> 固定電流closureの追跡電荷と目標電荷を、species別channel別に窓合計で記録する。
  module procedure open_fixed_current_history_writer
  character(len=1024) :: path
  integer(i32) :: species_idx
  logical :: fixed_current_used

  history_opened = .false.
  history_unit = -1
  if (.not. app%write_output) return
  path = trim(app%output_dir)//'/fixed_current_history.csv'
  fixed_current_used = trim(lower_ascii(app%surface_current%model)) == 'zhao_stationary'
  do species_idx = 1_i32, app%n_particle_species
    if (.not. app%particle_species(species_idx)%enabled) cycle
    fixed_current_used = fixed_current_used .or. &
                         trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) == 'fixed_current'
  end do
  if (app%history_stride <= 0_i32 .or. .not. fixed_current_used) then
    if (.not. resumed) call delete_history_if_exists(path)
    return
  end if
  call ensure_output_dir(app%output_dir)
  call open_history_file( &
    path, resumed, &
    'batch,simulated_time_s,window_batches,species_idx,channel,tracked_charge_C,target_charge_C,target_over_tracked', &
    history_unit, 'Failed to open fixed-current history file.' &
    )
  history_opened = .true.
  end procedure open_fixed_current_history_writer

  !> 窓内でtargetを適用したchannelだけを書く。追跡電荷が0なら比はnanとする。
  module procedure write_fixed_current_history_rows
  character(len=*), parameter :: channel_names(3) = [character(len=8) :: 'absorbed', 'emission', 'escape']
  character(len=32) :: ratio_text
  integer(i32) :: species_idx, channel

  do species_idx = 1_i32, int(size(channel_active, 2), i32)
    do channel = 1_i32, 3_i32
      if (.not. channel_active(channel, species_idx)) cycle
      if (tracked_charge_c(channel, species_idx) /= 0.0_dp) then
        write (ratio_text, '(es24.16)') target_charge_c(channel, species_idx)/tracked_charge_c(channel, species_idx)
      else
        ratio_text = 'nan'
      end if
      write (unit_id, '(i0,a,es24.16,2(a,i0),a,a,2(a,es24.16),a,a)') &
        batch_idx, ',', simulated_time_s, ',', window_batches, ',', species_idx, ',', trim(channel_names(channel)), &
        ',', tracked_charge_c(channel, species_idx), ',', target_charge_c(channel, species_idx), &
        ',', trim(adjustl(ratio_text))
    end do
  end do
  end procedure write_fixed_current_history_rows

  !> Append or replace a history file; an existing resumed file keeps its original header.
  subroutine open_history_file(path, resumed, header, unit_id, failure_message)
    character(len=*), intent(in) :: path, header, failure_message
    logical, intent(in) :: resumed
    integer, intent(out) :: unit_id
    integer :: ios
    logical :: exists

    inquire (file=trim(path), exist=exists)
    if (resumed) then
      open (newunit=unit_id, file=trim(path), status='unknown', position='append', action='write', iostat=ios)
    else
      open (newunit=unit_id, file=trim(path), status='replace', action='write', iostat=ios)
    end if
    if (ios /= 0) error stop failure_message
    if (.not. resumed .or. .not. exists) write (unit_id, '(a)') header
  end subroutine open_history_file

end submodule bem_output_writer_history
