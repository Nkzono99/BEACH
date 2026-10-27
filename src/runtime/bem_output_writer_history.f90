!> 履歴ファイルの生成・再開とスナップショットの出力。
submodule(bem_output_writer) bem_output_writer_history
  use bem_string_utils, only: lower_ascii
  use bem_filesystem, only: atomic_rename, remove_file_if_exists
  use bem_performance_profile, only: performance_profile_file
  use bem_outflow_spectrum, only: outflow_spectrum_bins, outflow_spectrum_edges_ev
  implicit none
  ! Released 2.0.x/3.0.x headers (17/19-column matching, 8-column top) are upgraded on resume;
  ! remove them when schema-10 checkpoints stop loading.
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

  history_path = trim(app%output_dir)//'/'//charge_history_file
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

  path = trim(app%output_dir)//'/'//potential_history_file
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

  path = trim(app%output_dir)//'/'//top_reference_history_file
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
  path = trim(app%output_dir)//'/'//matching_plane_history_file
  if (app%history_stride <= 0_i32 .or. &
      trim(lower_ascii(app%surface_current%model)) /= 'zhao_stationary' .or. &
      app%surface_current%outflow_refresh_batches <= 0_i32) return
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

  module procedure write_matching_plane_history_snapshot

  if (.not. stats%matching_plane%state_valid) return
  write (unit_id, '(i0,14(a,es24.16),a,i0,5(a,es24.16))') &
    batch_idx, ',', simulated_time_s, ',', stats%matching_plane%displacement_c_m2, &
    ',', stats%matching_plane%phi_v, ',', stats%matching_plane%response(2), &
    ',', stats%matching_plane%response(3), ',', stats%matching_plane%response(4), &
    ',', stats%matching_plane%response(5), ',', stats%matching_plane%response(6), &
    ',', stats%matching_plane%feedback(1), &
    ',', stats%matching_plane%feedback(2), ',', stats%matching_plane%feedback(3), &
    ',', stats%matching_plane%feedback(4), ',', stats%matching_plane%photoelectron_return_flux_m2_s, &
    ',', stats%matching_plane%photoelectron_escape_flux_m2_s, ',', stats%matching_plane%iterations, &
    ',', stats%matching_plane%residual, ',', stats%matching_plane%secondary_share, &
    ',', stats%matching_plane%secondary_budget_residual_a_m2, ',', stats%matching_plane%top_field_v_m, &
    ',', stats%matching_plane%top_cell_mean_v_m
  end procedure write_matching_plane_history_snapshot

  !> 外部根の更新が使う H 外向きスペクトルを、matching-plane history と同じ条件と間隔で書く。
  module procedure open_matching_plane_spectrum_history_writer
  character(len=1024) :: path

  history_opened = .false.
  history_unit = -1
  if (.not. app%write_output) return
  path = trim(app%output_dir)//'/'//matching_plane_spectrum_history_file
  if (app%history_stride <= 0_i32 .or. &
      trim(lower_ascii(app%surface_current%model)) /= 'zhao_stationary' .or. &
      app%surface_current%outflow_refresh_batches <= 0_i32) return
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

  if (.not. stats%matching_plane%state_valid) return
  edges = outflow_spectrum_edges_ev()
  do j = 1_i32, outflow_spectrum_bins
    if (stats%matching_plane%photoelectron_spectrum_m2_s(j) <= 0.0_dp) cycle
    write (unit_id, '(i0,a,es24.16,a,i0,3(a,es24.16))') batch_idx, ',', simulated_time_s, ',', j, &
      ',', edges(j), ',', edges(j + 1), ',', stats%matching_plane%photoelectron_spectrum_m2_s(j)
  end do
  end procedure write_matching_plane_spectrum_history_rows

  module procedure open_secondary_emission_history_writer
  character(len=1024) :: path
  logical :: active
  integer :: s

  history_opened = .false.
  history_unit = -1
  if (.not. app%write_output) return
  path = trim(app%output_dir)//'/'//secondary_emission_history_file
  active = .false.
  do s = 1, app%n_particle_species
    active = active .or. (app%particle_species(s)%enabled .and. app%particle_species(s)%secondary_emission_enabled)
  end do
  if (.not. active) return
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
  path = trim(app%output_dir)//'/'//fixed_current_history_file
  fixed_current_used = trim(lower_ascii(app%surface_current%model)) == 'zhao_stationary'
  do species_idx = 1_i32, app%n_particle_species
    if (.not. app%particle_species(species_idx)%enabled) cycle
    fixed_current_used = fixed_current_used .or. &
                         trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) == 'fixed_current'
  end do
  if (app%history_stride <= 0_i32 .or. .not. fixed_current_used) return
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

  !> fresh run の開始時に、前の run が書き今回は書かないかもしれない任意出力を消す。
  !! 既知の名前だけを対象にし、出力先にある利用者のファイルには触れない。
  module procedure remove_previous_run_outputs
  integer :: k

  do k = 1, size(history_files)
    call remove_output_file(trim(out_dir)//'/'//trim(history_files(k)))
  end do
  call remove_output_file(trim(out_dir)//'/'//mesh_potential_file)
  call remove_output_file(trim(out_dir)//'/'//performance_profile_file)
  end procedure remove_previous_run_outputs

  subroutine remove_output_file(path)
    character(len=*), intent(in) :: path
    integer :: status

    call remove_file_if_exists(path, status)
    if (status /= filesystem_success) error stop 'Failed to remove a previous-run output file: '//trim(path)
  end subroutine remove_output_file

  !> 再開した checkpoint の batch より後の行と、末尾の改行のない不完全な行を各 history から除く。
  !! 全 history は 1 列目が batch で、1 行目を header とする。
  !! checkpoint batch 自身の snapshot は完全である。batch loop は history 行を書いてから checkpoint を書き、
  !! checkpoint を書く batch では history unit を flush してから完了 manifest を公開するので、
  !! 再開できる checkpoint の batch 以下の行は、checkpoint が見える時点で全て改行まで書かれている。
  module procedure trim_history_files
  integer :: k

  do k = 1, size(history_files)
    call trim_history_file(trim(out_dir)//'/'//trim(history_files(k)), last_batch)
  end do
  end procedure trim_history_files

  !> 除く行があるときだけ、残す行を一時ファイルへ写して原子的に置き換える。
  !! 0 byte のファイルと存在しないファイルはそのまま残し、header は writer が開くときに書く。
  subroutine trim_history_file(path, last_batch)
    character(len=*), intent(in) :: path
    integer(i32), intent(in) :: last_batch

    character(len=:), allocatable :: temporary
    integer(i64) :: file_size
    integer :: output, ios, status
    logical :: exists, changed

    inquire (file=path, exist=exists, size=file_size)
    if (.not. exists .or. file_size <= 0_i64) return
    call filter_history_rows(path, file_size, last_batch, -1, changed)
    if (.not. changed) return
    temporary = path//'.trim.tmp'
    open (newunit=output, file=temporary, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to create a trimmed history file: '//path
    call filter_history_rows(path, file_size, last_batch, output, changed)
    close (output, iostat=ios)
    if (ios /= 0) error stop 'Failed to close a trimmed history file: '//path
    call atomic_rename(temporary, path, status)
    if (status /= filesystem_success) error stop 'Failed to replace a trimmed history file: '//path
  end subroutine trim_history_file

  !> history を stream で読み、header と batch <= last_batch の完全な行だけを残す。
  !! output が -1 なら何も書かず、除く行があるかだけを返す。
  subroutine filter_history_rows(path, file_size, last_batch, output, changed)
    character(len=*), intent(in) :: path
    integer(i64), intent(in) :: file_size
    integer(i32), intent(in) :: last_batch
    integer, intent(in) :: output
    logical, intent(out) :: changed

    integer, parameter :: chunk_len = 65536
    character(len=chunk_len) :: chunk
    character(len=:), allocatable :: carry
    integer(i64) :: next_pos
    integer :: input, ios, n, start, newline
    logical :: header_seen

    changed = .false.
    header_seen = .false.
    carry = ''
    open (newunit=input, file=path, access='stream', form='unformatted', status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Failed to read a history file for resume: '//path
    next_pos = 1_i64
    do while (next_pos <= file_size)
      n = int(min(int(chunk_len, i64), file_size - next_pos + 1_i64))
      read (input, pos=next_pos, iostat=ios) chunk(1:n)
      if (ios /= 0) error stop 'Failed to read a history file for resume: '//path
      next_pos = next_pos + int(n, i64)
      start = 1
      do
        newline = index(chunk(start:n), achar(10))
        if (newline == 0) exit
        if (len(carry) == 0) then
          call filter_line(chunk(start:start + newline - 2))
        else
          call filter_line(carry//chunk(start:start + newline - 2))
          carry = ''
        end if
        start = start + newline
      end do
      if (start <= n) carry = carry//chunk(start:n)
    end do
    close (input)
    ! 改行で終わらない末尾は書きかけの行なので捨てる。
    if (len(carry) > 0) changed = .true.

  contains

    subroutine filter_line(line)
      character(len=*), intent(in) :: line

      if (header_seen) then
        if (.not. row_batch_at_most(line, last_batch)) then
          changed = .true.
          return
        end if
      end if
      header_seen = .true.
      if (output /= -1) write (output, '(a)') line
    end subroutine filter_line

  end subroutine filter_history_rows

  !> 行頭の batch（数字列に ',' が続く）が last_batch 以下なら真。解釈できない行は残さない。
  pure logical function row_batch_at_most(line, last_batch) result(keep)
    character(len=*), intent(in) :: line
    integer(i32), intent(in) :: last_batch

    integer(i64) :: batch
    integer :: j, digit

    keep = .false.
    batch = 0_i64
    do j = 1, len(line)
      if (line(j:j) == ',') exit
      digit = iachar(line(j:j)) - iachar('0')
      if (digit < 0 .or. digit > 9 .or. j > 18) return
      batch = 10_i64*batch + int(digit, i64)
    end do
    if (j == 1 .or. j > len(line)) return
    keep = batch <= int(last_batch, i64)
  end function row_batch_at_most

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

  !> Append or replace a history file; a resumed non-empty file keeps its original header.
  !! A resumed file that is missing or empty (the previous run stopped before writing its header) gets one.
  subroutine open_history_file(path, resumed, header, unit_id, failure_message)
    character(len=*), intent(in) :: path, header, failure_message
    logical, intent(in) :: resumed
    integer, intent(out) :: unit_id
    integer :: ios
    integer(i64) :: file_size
    logical :: exists

    inquire (file=trim(path), exist=exists, size=file_size)
    if (resumed) then
      open (newunit=unit_id, file=trim(path), status='unknown', position='append', action='write', iostat=ios)
    else
      open (newunit=unit_id, file=trim(path), status='replace', action='write', iostat=ios)
    end if
    if (ios /= 0) error stop failure_message
    if (.not. resumed .or. .not. exists .or. file_size == 0_i64) write (unit_id, '(a)') header
  end subroutine open_history_file

end submodule bem_output_writer_history
