!> 実行サマリ・最終CSV・履歴CSVの出力を担当するモジュール。
module bem_output_writer
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, sim_stats, sim_config, bc_periodic
  use bem_app_config_types, only: app_config
  use bem_coulomb_fmm_core, only: eval_potential_points
  use bem_coulomb_fmm_types, only: fmm_plan_type, fmm_state_type, reset_fmm_plan, initialize_fmm_state, reset_fmm_state
  use bem_coulomb_fmm_periodic, only: use_periodic2_m2l_root_oracle
  use bem_coulomb_fmm_periodic_ewald, only: precompute_periodic2_ewald_data, &
                                            add_periodic2_exact_ewald_potential_correction_all_sources
  implicit none
  private

  real(dp), parameter :: pi_dp = acos(-1.0d0)

  public :: open_history_writer
  public :: print_run_summary
  public :: write_result_files
  public :: ensure_output_dir
  public :: compute_mesh_potential_fmm_core

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

  !> 実行結果の主要統計を標準出力へ表示する。
  !! @param[in] mesh 実行後のメッシュ情報。
  !! @param[in] stats 実行後の統計値。
  subroutine print_run_summary(mesh, stats)
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats

    print '(a,i0)', 'mesh nelem=', mesh%nelem
    print '(a,i0)', 'processed_particles=', stats%processed_particles
    print '(a,i0)', 'absorbed=', stats%absorbed
    print '(a,i0)', 'escaped=', stats%escaped
    print '(a,i0)', 'batches=', stats%batches
    print '(a,i0)', 'escaped_boundary=', stats%escaped_boundary
    print '(a,i0)', 'survived_max_step=', stats%survived_max_step
    print '(a,es12.4)', 'last_rel_change=', stats%last_rel_change
    print '(a,*(es12.4,1x))', 'mesh charges=', mesh%q_elem
  end subroutine print_run_summary

  !> 解析結果を `summary.txt` / `charges.csv` / `mesh_triangles.csv` などとして保存する。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 書き出し対象のメッシュ。
  !! @param[in] stats 書き出し対象の統計値。
  !! @param[in] cfg 出力設定を含むアプリ設定。
  subroutine write_result_files(out_dir, mesh, stats, cfg, mpi_world_size, mesh_potential_v)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    type(app_config), intent(in) :: cfg
    integer(i32), intent(in), optional :: mpi_world_size
    real(dp), intent(in), optional :: mesh_potential_v(:)

    call ensure_output_dir(out_dir)
    call write_summary_file(out_dir, mesh, stats, mpi_world_size=mpi_world_size)
    call write_charges_file(out_dir, mesh)
    if (cfg%write_mesh_potential) then
      if (present(mesh_potential_v)) then
        call write_mesh_potential_file(out_dir, mesh, cfg%sim, precomputed_potential_v=mesh_potential_v)
      else
        call write_mesh_potential_file(out_dir, mesh, cfg%sim)
      end if
    end if
    call write_mesh_file(out_dir, mesh)
    call write_mesh_sources_file(out_dir, mesh, cfg)
  end subroutine write_result_files

  !> 出力ディレクトリを作成する。
  !! @param[in] out_dir 作成対象ディレクトリのパス。
  subroutine ensure_output_dir(out_dir)
    character(len=*), intent(in) :: out_dir
    character(len=1024) :: cmd
    integer :: ios

    cmd = 'mkdir -p "'//trim(out_dir)//'"'
    call execute_command_line(trim(cmd), wait=.true., exitstat=ios)
    if (ios /= 0) error stop 'Failed to create output directory.'
  end subroutine ensure_output_dir

  !> 実行統計を `summary.txt` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh メッシュ情報（要素数を書き出す）。
  !! @param[in] stats 実行統計。
  subroutine write_summary_file(out_dir, mesh, stats, mpi_world_size)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    integer(i32), intent(in), optional :: mpi_world_size
    character(len=1024) :: summary_path
    integer :: u, ios
    integer(i32) :: world_size

    summary_path = trim(out_dir)//'/summary.txt'
    open (newunit=u, file=trim(summary_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open summary file.'
    world_size = 1_i32
    if (present(mpi_world_size)) world_size = max(1_i32, mpi_world_size)
    write (u, '(a,i0)') 'mesh_nelem=', mesh%nelem
    write (u, '(a,i0)') 'mesh_count=', max(1_i32, maxval(mesh%elem_mesh_id))
    write (u, '(a,i0)') 'mpi_world_size=', world_size
    write (u, '(a,i0)') 'processed_particles=', stats%processed_particles
    write (u, '(a,i0)') 'absorbed=', stats%absorbed
    write (u, '(a,i0)') 'escaped=', stats%escaped
    write (u, '(a,i0)') 'batches=', stats%batches
    write (u, '(a,i0)') 'escaped_boundary=', stats%escaped_boundary
    write (u, '(a,i0)') 'survived_max_step=', stats%survived_max_step
    write (u, '(a,es24.16)') 'last_rel_change=', stats%last_rel_change
    close (u)
  end subroutine write_summary_file

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

  !> 要素重心で評価した電位を `mesh_potential.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素重心と電荷を含むメッシュ情報。
  !! @param[in] sim 場の境界条件を含む実行設定。
  subroutine write_mesh_potential_file(out_dir, mesh, sim, precomputed_potential_v)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    real(dp), intent(in), optional :: precomputed_potential_v(:)
    character(len=1024) :: potential_path
    real(dp), allocatable :: potential_v(:)
    integer :: u, ios, i

    allocate (potential_v(mesh%nelem))
    if (present(precomputed_potential_v)) then
      if (size(precomputed_potential_v) /= mesh%nelem) error stop 'precomputed mesh potential size mismatch.'
      potential_v = precomputed_potential_v
    else
      call compute_mesh_potential(mesh, sim, potential_v)
    end if

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
    logical :: has_obj
    integer :: u, ios, i, mesh_id
    integer(i32) :: elem_count

    path = trim(out_dir)//'/mesh_sources.csv'
    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open mesh_sources.csv.'
    write (u, '(a)') 'mesh_id,source_kind,template_kind,elem_count'

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
      write (u, '(i0,a,a,a,a,a,i0)') 1, ',', 'obj', ',', 'obj', ',', mesh%nelem
      close (u)
      return
    end if

    mesh_id = 0
    do i = 1, size(cfg%templates)
      if (.not. cfg%templates(i)%enabled) cycle
      mesh_id = mesh_id + 1
      template_kind = trim(lower_ascii(cfg%templates(i)%kind))
      elem_count = int(count(mesh%elem_mesh_id == mesh_id), kind=i32)
      write (u, '(i0,a,a,a,a,a,i0)') mesh_id, ',', 'template', ',', trim(template_kind), ',', elem_count
    end do
    close (u)
  end subroutine write_mesh_sources_file

  !> 構築済み Coulomb FMM core を使って要素重心電位を計算する。
  !! `periodic2 + m2l_root_oracle` では root local の定数 mode が 0 固定のため、
  !! 1 点の exact 電位で gauge を補正する。
  subroutine compute_mesh_potential_fmm_core(mesh, plan, state, potential_v)
    type(mesh_type), intent(in) :: mesh
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    real(dp), intent(out) :: potential_v(:)

    real(dp) :: phi_anchor_exact, phi_offset, self_coeff, softening
    integer(i32) :: i

    if (size(potential_v) /= mesh%nelem) error stop 'mesh potential output array size mismatch.'
    if (.not. plan%built .or. .not. state%ready) error stop 'FMM core is not ready for mesh potential output.'
    if (plan%nsrc /= mesh%nelem) error stop 'FMM source count does not match mesh size.'

    potential_v = 0.0d0
    if (mesh%nelem <= 0) return

    softening = plan%options%softening
    if (softening < 0.0d0) error stop 'mesh potential output requires non-negative FMM softening.'

    call eval_potential_points(plan, state, mesh%centers, potential_v)

    if (use_periodic2_m2l_root_oracle(plan)) then
      call compute_exact_potential_point_from_core(plan, state, mesh%centers(:, 1), phi_anchor_exact)
      phi_offset = phi_anchor_exact - potential_v(1)
      potential_v = potential_v + phi_offset
    end if

    if (softening <= 0.0d0) then
      do i = 1, mesh%nelem
        self_coeff = 2.0d0*sqrt(pi_dp)/max(mesh%h_elem(i), sqrt(tiny(1.0d0)))
        potential_v(i) = potential_v(i) + self_coeff*mesh%q_elem(i)
      end do
    end if

    potential_v = k_coulomb*potential_v
  end subroutine compute_mesh_potential_fmm_core

  !> FMM source/state を使って 1 点の exact 電位を O(N) で再評価する。
  !! periodic2 では explicit image shell と oracle residual を含む。
  subroutine compute_exact_potential_point_from_core(plan, state, target_in, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: target_in(3)
    real(dp), intent(out) :: phi

    integer(i32) :: j, img_i, img_j, axis1, axis2, nimg
    real(dp) :: target(3), shifted(3), dx, dy, dz, r2, soft2, min_dist2

    phi = 0.0d0
    if (.not. plan%built .or. .not. state%ready) return

    soft2 = plan%options%softening*plan%options%softening
    min_dist2 = tiny(1.0d0)
    target = target_in

    if (.not. plan%options%use_periodic2) then
      do j = 1, plan%nsrc
        if (abs(state%src_q(j)) <= tiny(1.0d0)) cycle
        dx = target(1) - plan%src_pos(1, j)
        dy = target(2) - plan%src_pos(2, j)
        dz = target(3) - plan%src_pos(3, j)
        r2 = dx*dx + dy*dy + dz*dz + soft2
        if (r2 <= min_dist2) cycle
        phi = phi + state%src_q(j)/sqrt(r2)
      end do
      return
    end if

    axis1 = plan%options%periodic_axes(1)
    axis2 = plan%options%periodic_axes(2)
    nimg = max(0_i32, plan%options%periodic_image_layers)
    target(axis1) = wrap_periodic_coord(target(axis1), plan%options%target_box_min(axis1), plan%options%periodic_len(1))
    target(axis2) = wrap_periodic_coord(target(axis2), plan%options%target_box_min(axis2), plan%options%periodic_len(2))

    do j = 1, plan%nsrc
      if (abs(state%src_q(j)) <= tiny(1.0d0)) cycle
      dx = target(1) - plan%src_pos(1, j)
      dy = target(2) - plan%src_pos(2, j)
      dz = target(3) - plan%src_pos(3, j)
      r2 = dx*dx + dy*dy + dz*dz + soft2
      if (r2 > min_dist2) phi = phi + state%src_q(j)/sqrt(r2)
    end do

    do img_i = -nimg, nimg
      do img_j = -nimg, nimg
        if (img_i == 0_i32 .and. img_j == 0_i32) cycle
        do j = 1, plan%nsrc
          if (abs(state%src_q(j)) <= tiny(1.0d0)) cycle
          shifted = plan%src_pos(:, j)
          shifted(axis1) = shifted(axis1) + real(img_i, dp)*plan%options%periodic_len(1)
          shifted(axis2) = shifted(axis2) + real(img_j, dp)*plan%options%periodic_len(2)
          dx = target(1) - shifted(1)
          dy = target(2) - shifted(2)
          dz = target(3) - shifted(3)
          r2 = dx*dx + dy*dy + dz*dz + soft2
          if (r2 <= min_dist2) cycle
          phi = phi + state%src_q(j)/sqrt(r2)
        end do
      end do
    end do

    if (use_periodic2_m2l_root_oracle(plan)) then
      call add_periodic2_exact_ewald_potential_correction_all_sources(plan, state, target, phi)
    end if
  end subroutine compute_exact_potential_point_from_core

  !> 要素重心での電位を計算する。
  !! `periodic2` では explicit image shell に加えて、
  !! `m2l_root_oracle` 相当の exact Ewald residual も反映する。
  !! @param[in] mesh 要素重心・代表長・電荷を含むメッシュ情報。
  !! @param[in] sim softening と境界条件を含む実行設定。
  !! @param[out] potential_v 各要素重心での電位 [V]。
  subroutine compute_mesh_potential(mesh, sim, potential_v)
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    real(dp), intent(out) :: potential_v(:)

    integer(i32) :: i, j, ix, iy, axis1, axis2, image_layers
    integer(i32) :: periodic_axes(2)
    real(dp) :: softening, soft2, min_dist2, self_coeff
    real(dp) :: periodic_len(2), periodic_origins(2)
    real(dp) :: target(3), shifted(3), dx, dy, dz, r2
    logical :: use_periodic2, use_exact_ewald_residual
    type(fmm_plan_type) :: ewald_plan
    type(fmm_state_type) :: ewald_state

    if (size(potential_v) /= mesh%nelem) error stop 'mesh potential output array size mismatch.'

    softening = sim%softening
    if (softening < 0.0d0) error stop 'sim.softening must be >= 0 for mesh potential output.'
    soft2 = softening*softening
    min_dist2 = tiny(1.0d0)
    potential_v = 0.0d0
    if (mesh%nelem <= 0) return

    call resolve_periodic2_output(sim, use_periodic2, periodic_axes, periodic_len, periodic_origins, image_layers)
    axis1 = periodic_axes(1)
    axis2 = periodic_axes(2)
    call initialize_fmm_state(ewald_state)
    call initialize_periodic2_ewald_context( &
      mesh, sim, use_periodic2, periodic_axes, periodic_len, image_layers, ewald_plan, ewald_state, use_exact_ewald_residual &
      )

    do i = 1, mesh%nelem
      if (softening > 0.0d0) then
        self_coeff = 1.0d0/softening
      else
        self_coeff = 2.0d0*sqrt(pi_dp)/max(mesh%h_elem(i), sqrt(tiny(1.0d0)))
      end if
      potential_v(i) = self_coeff*mesh%q_elem(i)
    end do

    if (.not. use_periodic2) then
      do i = 1, mesh%nelem
        do j = 1, mesh%nelem
          if (j == i) cycle
          dx = mesh%center_x(i) - mesh%center_x(j)
          dy = mesh%center_y(i) - mesh%center_y(j)
          dz = mesh%center_z(i) - mesh%center_z(j)
          r2 = dx*dx + dy*dy + dz*dz + soft2
          potential_v(i) = potential_v(i) + mesh%q_elem(j)/sqrt(max(r2, min_dist2))
        end do
      end do
      potential_v = k_coulomb*potential_v
      return
    end if

    do i = 1, mesh%nelem
      target = mesh%centers(:, i)
      target(axis1) = wrap_periodic_coord(target(axis1), periodic_origins(1), periodic_len(1))
      target(axis2) = wrap_periodic_coord(target(axis2), periodic_origins(2), periodic_len(2))

      do j = 1, mesh%nelem
        if (j == i) cycle
        dx = target(1) - mesh%centers(1, j)
        dy = target(2) - mesh%centers(2, j)
        dz = target(3) - mesh%centers(3, j)
        r2 = dx*dx + dy*dy + dz*dz + soft2
        potential_v(i) = potential_v(i) + mesh%q_elem(j)/sqrt(max(r2, min_dist2))
      end do

      do ix = -image_layers, image_layers
        do iy = -image_layers, image_layers
          if (ix == 0_i32 .and. iy == 0_i32) cycle
          do j = 1, mesh%nelem
            shifted = mesh%centers(:, j)
            shifted(axis1) = shifted(axis1) + real(ix, dp)*periodic_len(1)
            shifted(axis2) = shifted(axis2) + real(iy, dp)*periodic_len(2)
            dx = target(1) - shifted(1)
            dy = target(2) - shifted(2)
            dz = target(3) - shifted(3)
            r2 = dx*dx + dy*dy + dz*dz + soft2
            potential_v(i) = potential_v(i) + mesh%q_elem(j)/sqrt(max(r2, min_dist2))
          end do
        end do
      end do

      if (use_exact_ewald_residual) then
        call add_periodic2_exact_ewald_potential_correction_all_sources(ewald_plan, ewald_state, target, potential_v(i))
      end if
    end do

    if (use_exact_ewald_residual) then
      call reset_fmm_state(ewald_state)
      call reset_fmm_plan(ewald_plan)
    end if
    potential_v = k_coulomb*potential_v
  end subroutine compute_mesh_potential

  !> exact Ewald residual 用の最小 periodic2 context を作る。
  subroutine initialize_periodic2_ewald_context( &
    mesh, sim, use_periodic2, periodic_axes, periodic_len, image_layers, plan, state, enabled &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    logical, intent(in) :: use_periodic2
    integer(i32), intent(in) :: periodic_axes(2), image_layers
    real(dp), intent(in) :: periodic_len(2)
    type(fmm_plan_type), intent(inout) :: plan
    type(fmm_state_type), intent(inout) :: state
    logical, intent(out) :: enabled

    character(len=16) :: far_correction
    integer(i32) :: ewald_layers

    enabled = .false.
    call reset_fmm_plan(plan)
    call reset_fmm_state(state)
    if (.not. use_periodic2) return

    far_correction = trim(lower_ascii(sim%field_periodic_far_correction))
    ewald_layers = max(0_i32, sim%field_periodic_ewald_layers)
    select case (trim(far_correction))
    case ('auto', 'none')
      far_correction = 'm2l_root_oracle'
      ewald_layers = max(1_i32, ewald_layers)
    case ('m2l_root_oracle')
      ewald_layers = max(1_i32, ewald_layers)
    case default
      return
    end select

    plan%options%use_periodic2 = .true.
    plan%options%periodic_far_correction = far_correction
    plan%options%periodic_axes = periodic_axes
    plan%options%periodic_len = periodic_len
    plan%options%periodic_image_layers = image_layers
    plan%options%periodic_ewald_alpha = sim%field_periodic_ewald_alpha
    plan%options%periodic_ewald_layers = ewald_layers
    plan%options%softening = sim%softening
    plan%options%target_box_min = sim%box_min
    plan%options%target_box_max = sim%box_max
    plan%nsrc = mesh%nelem
    allocate (plan%src_pos(3, mesh%nelem))
    plan%src_pos = mesh%centers
    allocate (state%src_q(mesh%nelem))
    state%src_q = mesh%q_elem
    state%ready = .true.

    call precompute_periodic2_ewald_data(plan)
    enabled = plan%periodic_ewald%ready
    if (.not. enabled) then
      call reset_fmm_state(state)
      call reset_fmm_plan(plan)
    end if
  end subroutine initialize_periodic2_ewald_context

  !> `periodic2` 用の周期軸・周期長・原点・画像層数を解決する。
  subroutine resolve_periodic2_output(sim, use_periodic2, periodic_axes, periodic_len, periodic_origins, image_layers)
    type(sim_config), intent(in) :: sim
    logical, intent(out) :: use_periodic2
    integer(i32), intent(out) :: periodic_axes(2)
    real(dp), intent(out) :: periodic_len(2)
    real(dp), intent(out) :: periodic_origins(2)
    integer(i32), intent(out) :: image_layers

    character(len=16) :: field_bc_mode
    integer(i32) :: axis, n_periodic
    real(dp) :: span

    use_periodic2 = .false.
    periodic_axes = 0_i32
    periodic_len = 0.0d0
    periodic_origins = 0.0d0
    image_layers = max(0_i32, sim%field_periodic_image_layers)

    field_bc_mode = trim(lower_ascii(sim%field_bc_mode))
    if (field_bc_mode /= 'periodic2') return
    if (.not. sim%use_box) error stop 'mesh potential output for periodic2 requires sim.use_box=true.'

    n_periodic = 0_i32
    do axis = 1_i32, 3_i32
      if ((sim%bc_low(axis) == bc_periodic) .neqv. (sim%bc_high(axis) == bc_periodic)) then
        error stop 'mesh potential output requires bc_low(axis)=bc_high(axis)=periodic for periodic axes.'
      end if
      if (sim%bc_low(axis) == bc_periodic) then
        n_periodic = n_periodic + 1_i32
        if (n_periodic <= 2_i32) periodic_axes(n_periodic) = axis
      end if
    end do
    if (n_periodic /= 2_i32) then
      error stop 'mesh potential output with sim.field_bc_mode="periodic2" requires exactly two periodic axes.'
    end if

    do axis = 1_i32, 2_i32
      span = sim%box_max(periodic_axes(axis)) - sim%box_min(periodic_axes(axis))
      if (span <= 0.0d0) error stop 'mesh potential output requires positive box length on periodic axes.'
      periodic_len(axis) = span
      periodic_origins(axis) = sim%box_min(periodic_axes(axis))
    end do

    use_periodic2 = .true.
  end subroutine resolve_periodic2_output

  !> ASCII 英字を小文字化する。
  pure function lower_ascii(s) result(out)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: out
    integer :: i, code

    out = s
    do i = 1, len(s)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        out(i:i) = achar(code + 32)
      end if
    end do
  end function lower_ascii

  !> 周期軸上の座標を `[origin, origin + length)` に折り返す。
  pure real(dp) function wrap_periodic_coord(x, origin, periodic_len) result(wrapped)
    real(dp), intent(in) :: x, origin, periodic_len

    wrapped = origin + modulo(x - origin, periodic_len)
  end function wrap_periodic_coord

end module bem_output_writer
