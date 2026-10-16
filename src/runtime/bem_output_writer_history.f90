!> 履歴ファイルの生成・再開とスナップショットの出力。
submodule(bem_output_writer) bem_output_writer_history
  use bem_string_utils, only: lower_ascii
  use bem_filesystem, only: atomic_rename
  use bem_outflow_spectrum, only: outflow_spectrum_bins, outflow_spectrum_edges_ev
  implicit none
  character(len=*), parameter :: matching_history_legacy_header = &
                                 'batch,simulated_time_s,D_H_C_m2,phi_H_V,electron_inward_flux_m2_s,ion_inward_flux_m2_s,'// &
                                 'electron_access_potential_V,ion_access_potential_V,photoelectron_barrier_potential_V,'// &
                                 'photoelectron_outward_flux_m2_s,'// &
                                 'photoelectron_mean_normal_energy_eV,electron_outward_flux_m2_s,ion_outward_flux_m2_s,'// &
                                 'photoelectron_return_flux_m2_s,photoelectron_escape_flux_m2_s,iterations,residual'
  character(len=*), parameter :: matching_history_v3_header = &
                                 matching_history_legacy_header//',secondary_share,secondary_budget_residual_a_m2'
  character(len=*), parameter :: matching_history_header = &
                                 matching_history_v3_header//',E_top_sheath_V_m,E_top_cell_mean_V_m'
  character(len=*), parameter :: top_history_v3_header = &
                                 'batch,simulated_time_s,z_high_m,sample_n,potential_mean_V,potential_std_V,'// &
                                 'potential_min_V,potential_max_V'
  character(len=*), parameter :: top_history_header = &
                                 top_history_v3_header//',switch_phi_rms_mismatch_V,switch_Ez_rms_mismatch_V_m'
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
  if (resumed) call extend_history_columns(path, top_history_v3_header, top_history_header, ',0,0')
  call open_history_file( &
    path, resumed, &
    top_history_header, &
    history_unit, 'Failed to open top-reference potential history file.' &
    )
  history_opened = .true.
  end procedure open_top_reference_history_writer

  !> z-high 面平均電位基準を履歴 CSV へ1行書き出す。
  module procedure write_top_reference_history_snapshot

  real(dp) :: phi_rms, ez_rms

  phi_rms = 0.0_dp
  ez_rms = 0.0_dp
  if (present(switch_phi_rms)) phi_rms = switch_phi_rms
  if (present(switch_ez_rms)) ez_rms = switch_ez_rms
  write (unit_id, '(i0,2(a,es24.16),a,i0,6(a,es24.16))') &
    batch_idx, ',', simulated_time_s, ',', z_high_m, ',', sample_n, &
    ',', potential_mean_v, ',', potential_std_v, ',', potential_min_v, ',', potential_max_v, ',', phi_rms, ',', ez_rms
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
  if (resumed) then
    call extend_history_columns(path, matching_history_legacy_header, matching_history_header, ',0,0,NaN,NaN')
    call extend_history_columns(path, matching_history_v3_header, matching_history_header, ',NaN,NaN')
  end if
  call open_history_file( &
    path, resumed, matching_history_header, &
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
  write (unit_id, '(i0,14(a,es24.16),a,i0,5(a,es24.16))') &
    batch_idx, ',', simulated_time_s, ',', stats%matching_plane_displacement_c_m2, &
    ',', stats%matching_plane_phi_v, ',', stats%matching_plane_response(2), &
    ',', stats%matching_plane_response(3), ',', stats%matching_plane_response(4), &
    ',', stats%matching_plane_response(5), ',', stats%matching_plane_response(6), &
    ',', stats%matching_plane_feedback(1), &
    ',', stats%matching_plane_feedback(2), ',', stats%matching_plane_feedback(3), &
    ',', stats%matching_plane_feedback(4), ',', stats%matching_plane_photoelectron_return_flux_m2_s, &
    ',', stats%matching_plane_photoelectron_escape_flux_m2_s, ',', stats%matching_plane_iterations, &
    ',', stats%matching_plane_residual, ',', stats%matching_plane_secondary_share, &
    ',', stats%secondary_budget_residual_a_m2, ',', stats%surface_current_top_field_v_m, &
    ',', stats%matching_plane_top_cell_mean_v_m
  end procedure write_matching_plane_history_snapshot

  !> 外部根の更新が使う H 外向きスペクトルを、matching-plane history と同じ条件と間隔で書く。
  module procedure open_matching_plane_spectrum_history_writer
  character(len=1024) :: path

  history_opened = .false.
  history_unit = -1
  if (.not. app%write_output) return
  path = trim(app%output_dir)//'/matching_plane_spectrum_history.csv'
  if (app%history_stride <= 0_i32 .or. &
      trim(lower_ascii(app%surface_current%model)) /= 'zhao_stationary' .or. &
      app%surface_current%outflow_refresh_batches <= 0_i32) then
    if (.not. resumed) call delete_history_if_exists(path)
    return
  end if
  call ensure_output_dir(app%output_dir)
  call open_history_file( &
    path, resumed, 'batch,simulated_time_s,bin,energy_low_eV,energy_high_eV,photoelectron_outward_flux_m2_s', &
    history_unit, 'Failed to open matching-plane spectrum history file.' &
    )
  history_opened = .true.
  end procedure open_matching_plane_spectrum_history_writer

  !> 現在の外部源のうち流束が正の bin だけを書く。表面放出の半 Maxwell 源（初期根）では何も書かない。
  module procedure write_matching_plane_spectrum_history_rows
  real(dp) :: edges(outflow_spectrum_bins + 1)
  integer(i32) :: j

  if (.not. stats%matching_plane_state_valid) return
  edges = outflow_spectrum_edges_ev()
  do j = 1_i32, outflow_spectrum_bins
    if (stats%matching_plane_photoelectron_spectrum_m2_s(j) <= 0.0_dp) cycle
    write (unit_id, '(i0,a,es24.16,a,i0,3(a,es24.16))') batch_idx, ',', simulated_time_s, ',', j, &
      ',', edges(j), ',', edges(j + 1), ',', stats%matching_plane_photoelectron_spectrum_m2_s(j)
  end do
  end procedure write_matching_plane_spectrum_history_rows

  module procedure open_secondary_emission_history_writer
  character(len=1024) :: path
  logical :: active
  integer :: s

  history_opened = .false.
  history_unit = -1
  if (.not. app%write_output) return
  path = trim(app%output_dir)//'/secondary_emission_history.csv'
  active = .false.
  do s = 1, app%n_particle_species
    active = active .or. (app%particle_species(s)%enabled .and. app%particle_species(s)%secondary_emission_enabled)
  end do
  if (.not. active) then
    if (.not. resumed) call delete_history_if_exists(path)
    return
  end if
  call ensure_output_dir(app%output_dir)
  call open_history_file(path, resumed, &
                         'batch,incident_species,hits,emitted_count,emitted_charge_c,mean_yield,generation_capped_count', &
                         history_unit, 'Failed to open secondary-emission history file.')
  history_opened = .true.
  end procedure open_secondary_emission_history_writer

  module procedure write_secondary_emission_history_rows
  integer :: s, j
  real(dp) :: mean_yield
  character(len=:), allocatable :: key

  do s = 1, app%n_particle_species
    if (.not. app%particle_species(s)%enabled .or. .not. app%particle_species(s)%secondary_emission_enabled) cycle
    mean_yield = 0.0_dp
    if (counts(1, s) > 0_i64) mean_yield = values(2, s)/real(counts(1, s), dp)
    key = '"'
    do j = 1, len_trim(app%particle_species(s)%species_key)
      key = key//app%particle_species(s)%species_key(j:j)
      if (app%particle_species(s)%species_key(j:j) == '"') key = key//'"'
    end do
    key = key//'"'
    write (unit_id, '(i0,a,a,2(a,i0),2(a,es24.16),a,i0)') batch_idx, ',', key, &
      ',', counts(1, s), ',', counts(2, s), ',', values(1, s), ',', mean_yield, ',', counts(3, s)
  end do
  end procedure write_secondary_emission_history_rows

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

  !> Preserve released 3.0.0 CSV rows on resume (supported while these histories can be resumed).
  !! Missing historical fields are NaN; histories predating the switch have zero mismatch.
  subroutine extend_history_columns(path, old_header, new_header, missing_values)
    character(len=*), intent(in) :: path, old_header, new_header, missing_values
    character(len=8192) :: line
    character(len=:), allocatable :: temporary
    logical :: exists
    integer :: input, output, ios, status

    inquire (file=trim(path), exist=exists)
    if (.not. exists) return
    open (newunit=input, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Failed to read history for resume.'
    read (input, '(a)', iostat=ios) line
    if (ios /= 0 .or. trim(line) /= old_header) then
      close (input)
      return
    end if
    temporary = trim(path)//'.columns.tmp'
    open (newunit=output, file=temporary, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to extend history columns.'
    write (output, '(a)') new_header
    do
      read (input, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) error stop 'Failed to read legacy history row.'
      write (output, '(a)') trim(line)//missing_values
    end do
    close (input)
    close (output)
    call atomic_rename(temporary, trim(path), status)
    if (status /= filesystem_success) error stop 'Failed to replace legacy history.'
  end subroutine extend_history_columns

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
