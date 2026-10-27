!> 主ループへ渡す出力先（history の unit と間隔、終了時のメッシュ電位）を 1 つの型にまとめる。
!! unit が -1 の history は書かない。CLI は root rank で writer を開いてこの型へ詰め、終了時に閉じる。
module bem_run_outputs
  use bem_kinds, only: dp, i32
  use bem_app_config_types, only: app_config
  use bem_output_writer, only: open_history_writer, open_potential_history_writer, open_top_reference_history_writer, &
                               open_matching_plane_history_writer, open_fixed_current_history_writer, &
                               open_matching_plane_spectrum_history_writer, open_secondary_emission_history_writer, &
                               remove_previous_run_outputs, trim_history_files
  use bem_periodic_checkpoint, only: invalidate_output_checkpoints
  implicit none
  private

  public :: open_run_outputs
  public :: flush_run_outputs
  public :: close_run_outputs

  type, public :: run_outputs_type
    !> 電荷 history（output.history.stride_batches ごと）。
    integer :: history_unit = -1
    !> 要素電位 history。間隔は history_stride。
    integer :: potential_history_unit = -1
    !> z-high 参照電位 history。要素電位 history と同じ batch で書くので、potential_history_unit が有効なときだけ使う。
    integer :: top_reference_history_unit = -1
    integer :: matching_plane_history_unit = -1
    !> 外部源のスペクトル history。外部根を更新した batch で、間隔によらず書く。
    integer :: matching_plane_spectrum_history_unit = -1
    integer :: fixed_current_history_unit = -1
    !> SEE history。毎 batch 書く。
    integer :: secondary_emission_history_unit = -1
    !> history を書く accepted batch の間隔。0 以下なら matching-plane・spectrum・固定電流 history を書かず、
    !! 電荷・要素電位 history は毎 batch の扱いになる（CLI の writer はこのとき開かれない）。
    integer(i32) :: history_stride = 1_i32
    !> 終了時に要素重心の電位を mesh_potential_v へ求める（root rank だけ）。
    logical :: mesh_potential = .false.
    real(dp), allocatable :: mesh_potential_v(:)
  end type run_outputs_type

contains

  !> root rank で出力先をこの run のものにしてから、設定どおりに history writer を開く。root 以外は全て無効（unit=-1）。
  !! @param[in] resume_batch 再開した checkpoint の batch 数（`resumed` のときだけ使う）。
  subroutine open_run_outputs(app, resumed, resume_batch, is_root, outputs)
    type(app_config), intent(in) :: app
    logical, intent(in) :: resumed, is_root
    integer(i32), intent(in) :: resume_batch
    type(run_outputs_type), intent(out) :: outputs
    logical :: opened
    integer :: unit_id

    outputs%history_stride = app%history_stride
    outputs%mesh_potential = app%write_mesh_potential
    if (.not. is_root) return
    if (app%write_output) call claim_output_dir(app, resumed, resume_batch)
    call open_history_writer(app, resumed, opened, unit_id)
    if (opened) outputs%history_unit = unit_id
    call open_potential_history_writer(app, resumed, opened, unit_id)
    if (opened) outputs%potential_history_unit = unit_id
    call open_top_reference_history_writer(app, resumed, opened, unit_id)
    if (opened) outputs%top_reference_history_unit = unit_id
    call open_matching_plane_history_writer(app, resumed, opened, unit_id)
    if (opened) outputs%matching_plane_history_unit = unit_id
    call open_fixed_current_history_writer(app, resumed, opened, unit_id)
    if (opened) outputs%fixed_current_history_unit = unit_id
    call open_matching_plane_spectrum_history_writer(app, resumed, opened, unit_id)
    if (opened) outputs%matching_plane_spectrum_history_unit = unit_id
    call open_secondary_emission_history_writer(app, resumed, opened, unit_id)
    if (opened) outputs%secondary_emission_history_unit = unit_id
  end subroutine open_run_outputs

  !> 出力先に残る BEACH の成果物を、これから続ける run と同じ世代にそろえる。root rank だけが呼ぶ。
  !! 再開では、checkpoint より後の history 行（巻き戻した試行の行）を除き、checkpoint のファイルは変更しない。
  !! `run.restart.from` が別の directory なら、その中は何も変更しない。
  !! fresh run では、まず前の run の checkpoint を再開候補から外し、次に今回書かないかもしれない任意出力を消す。
  !! 出力先を読み書きするのは root だけで（checkpoint の slot 選択も root が行って broadcast する）、
  !! fresh run の他の rank は出力先を読まないので、追加の barrier は要らない。
  subroutine claim_output_dir(app, resumed, resume_batch)
    type(app_config), intent(in) :: app
    logical, intent(in) :: resumed
    integer(i32), intent(in) :: resume_batch

    if (resumed) then
      call trim_history_files(trim(app%output_dir), resume_batch)
    else
      call invalidate_output_checkpoints(trim(app%output_dir))
      call remove_previous_run_outputs(trim(app%output_dir))
    end if
  end subroutine claim_output_dir

  !> 開いている history の行を書き出す。checkpoint の完了 manifest を公開する前に呼び、
  !! 再開できる batch 以下の history 行がすべて改行まで書かれていることを保証する。
  subroutine flush_run_outputs(outputs)
    type(run_outputs_type), intent(in) :: outputs

    call flush_unit(outputs%history_unit)
    call flush_unit(outputs%potential_history_unit)
    call flush_unit(outputs%top_reference_history_unit)
    call flush_unit(outputs%matching_plane_history_unit)
    call flush_unit(outputs%fixed_current_history_unit)
    call flush_unit(outputs%secondary_emission_history_unit)
    call flush_unit(outputs%matching_plane_spectrum_history_unit)
  end subroutine flush_run_outputs

  subroutine flush_unit(unit_id)
    integer, intent(in) :: unit_id

    if (unit_id /= -1) flush (unit_id)
  end subroutine flush_unit

  !> 開いている history writer を閉じる。
  subroutine close_run_outputs(outputs)
    type(run_outputs_type), intent(inout) :: outputs

    call close_unit(outputs%history_unit)
    call close_unit(outputs%potential_history_unit)
    call close_unit(outputs%top_reference_history_unit)
    call close_unit(outputs%matching_plane_history_unit)
    call close_unit(outputs%fixed_current_history_unit)
    call close_unit(outputs%secondary_emission_history_unit)
    call close_unit(outputs%matching_plane_spectrum_history_unit)
  end subroutine close_run_outputs

  subroutine close_unit(unit_id)
    integer, intent(inout) :: unit_id

    if (unit_id == -1) return
    close (unit_id)
    unit_id = -1
  end subroutine close_unit

end module bem_run_outputs
