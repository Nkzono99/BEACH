!> Coulomb FMM plan 構築と tree トポロジ前計算。
module bem_coulomb_fmm_plan_ops
  use, intrinsic :: ieee_arithmetic, only: ieee_class, ieee_negative_zero, operator(==)
  use bem_filesystem, only: default_periodic_cache_directory
  use bem_kinds, only: dp, sp, i32
  use bem_coulomb_fmm_types, only: fmm_options_type, fmm_plan_type, fmm_tree_type, reset_fmm_plan, &
                                   periodic_source_proxy_half_size, periodic_operator_source_radius_ratio, &
                                   periodic_operator_target_radius_ratio, m2l_lane_width
  use bem_coulomb_fmm_basis, only: initialize_basis_tables, build_axis_powers, compute_laplace_derivatives
  use bem_coulomb_fmm_periodic, only: has_valid_target_box, build_periodic_shift_values, distance_to_source_bbox, &
                                      distance_to_source_bbox_periodic, wrap_src_pos_to_primary_cell, &
                                      use_periodic2_cached_kneq0
  use bem_coulomb_fmm_periodic_ewald, only: precompute_periodic2_ewald_data
  use bem_coulomb_fmm_periodic_root_ops, only: precompute_periodic_root_operator
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_periodic_zero_mode_plan, only: build_periodic_zero_mode_height_plan, periodic_zero_mode_ok
  use bem_periodic_shell_mean, only: build_periodic_shell_mean_plan
  use bem_coulomb_fmm_tree_utils, only: octant_index, append_i32_buffer, nodes_well_separated
  implicit none
  private

  public :: core_build_plan_impl
  public :: core_build_panel_plan_impl
  public :: core_prepare_panel_sources_impl
  public :: core_complete_plan_impl
  public :: core_destroy_plan_impl

contains

  !> FMM 計画と木構造、転送演算子を構築する。
  !! @param[inout] plan 構築対象の FMM 計画。
  !! @param[in] src_pos ソース位置 `(3,n)` [m]。
  !! @param[in] options FMM 設定。
  subroutine core_build_plan_impl(plan, src_pos, options)
    type(fmm_plan_type), intent(inout) :: plan
    real(dp), intent(in) :: src_pos(:, :)
    type(fmm_options_type), intent(in) :: options
    call prepare_plan_sources(plan, src_pos, options)
    call core_complete_plan_impl(plan)
  end subroutine core_build_plan_impl

  subroutine core_build_panel_plan_impl(plan, v0, v1, v2, options)
    type(fmm_plan_type), intent(inout) :: plan
    real(dp), intent(in) :: v0(:, :), v1(:, :), v2(:, :)
    type(fmm_options_type), intent(in) :: options

    call core_prepare_panel_sources_impl(plan, v0, v1, v2, options)
    call core_complete_plan_impl(plan)
  end subroutine core_build_panel_plan_impl

  !> Panel sources only: options, centroids and panel geometry, as the first part of a plan build.
  !! The plan then supports exact-direct evaluation; `core_complete_plan_impl` adds the FMM part.
  subroutine core_prepare_panel_sources_impl(plan, v0, v1, v2, options)
    type(fmm_plan_type), intent(inout) :: plan
    real(dp), intent(in) :: v0(:, :), v1(:, :), v2(:, :)
    type(fmm_options_type), intent(in) :: options
    real(dp), allocatable :: src_pos(:, :)

    if (size(v0, 1) /= 3 .or. size(v1, 1) /= 3 .or. size(v2, 1) /= 3 .or. &
        size(v0, 2) /= size(v1, 2) .or. size(v0, 2) /= size(v2, 2)) then
      error stop 'FMM panel core expects v0/v1/v2(3,n).'
    end if
    allocate (src_pos(3, size(v0, 2)))
    src_pos = (v0 + v1 + v2)/3.0_dp
    call prepare_plan_sources(plan, src_pos, options, v0, v1, v2)
  end subroutine core_prepare_panel_sources_impl

  subroutine prepare_plan_sources(plan, src_pos, options, v0, v1, v2)
    type(fmm_plan_type), intent(inout) :: plan
    real(dp), intent(in) :: src_pos(:, :)
    type(fmm_options_type), intent(in) :: options
    real(dp), intent(in), optional :: v0(:, :), v1(:, :), v2(:, :)
    integer(i32) :: nsrc, idx, panel_status
    real(dp), allocatable :: original_centroid(:, :)
    real(dp) :: shift(3), shifted_vertex(3, 3)

    if (allocated(plan%alpha) .or. allocated(plan%src_pos) .or. allocated(plan%source_tree%elem_order)) then
      call core_destroy_plan_impl(plan)
    end if
    plan%options = options
    nsrc = int(size(src_pos, 2), i32)
    plan%nsrc = nsrc
    if (size(src_pos, 1) /= 3) error stop 'FMM core expects src_pos(3,n).'
    if (options%leaf_max <= 0_i32) error stop 'FMM leaf_max must be > 0.'
    if (options%order < 1_i32) error stop 'FMM order must be >= 1.'
    if (options%softening < 0.0d0) error stop 'FMM softening must be >= 0.'
    if (options%use_periodic2) then
      if (any(options%periodic_axes < 1_i32) .or. any(options%periodic_axes > 3_i32)) then
        error stop 'periodic2 requires periodic axes in 1..3.'
      end if
      if (options%periodic_axes(1) == options%periodic_axes(2)) then
        error stop 'periodic2 requires two distinct axes.'
      end if
      if (any(options%periodic_len <= 0.0d0)) then
        error stop 'periodic2 requires positive periodic lengths.'
      end if
      if (.not. has_valid_target_box(plan%options)) then
        error stop 'periodic2 requires a valid target box.'
      end if
      select case (trim(plan%options%periodic_far_correction))
      case ('auto')
        plan%options%periodic_far_correction = 'none'
      case ('none')
        continue
      case ('m2l_root_oracle')
        error stop 'periodic far correction "m2l_root_oracle" was removed; use "cached_kneq0".'
      case ('cached_kneq0')
        continue
      case default
        error stop 'Unsupported periodic far correction in FMM core.'
      end select
      if (trim(plan%options%periodic_far_correction) == 'cached_kneq0') then
        if (plan%options%periodic_ewald_layers < 1_i32) then
          error stop 'periodic2 root-operator far correction requires periodic_ewald_layers >= 1.'
        end if
      end if
      if (trim(plan%options%periodic_far_correction) == 'cached_kneq0') then
        if (plan%options%periodic_image_layers < 1_i32) error stop 'cached_kneq0 requires at least one image layer.'
        if (len_trim(plan%options%periodic_cache_dir) == 0) then
          plan%options%periodic_cache_dir = default_periodic_cache_directory()
        end if
        if (plan%options%periodic_generation_tolerance <= 0.0_dp) then
          error stop 'cached_kneq0 requires a positive generation tolerance.'
        end if
      end if
    end if

    allocate (plan%src_pos(3, max(0_i32, nsrc)))
    if (nsrc > 0_i32) plan%src_pos = src_pos

    plan%panel_source = present(v0)
    if (plan%panel_source) then
      if (options%softening /= 0.0_dp) error stop 'FMM panel source requires zero softening.'
      allocate (plan%panel_geometry(nsrc), original_centroid(3, nsrc))
      do idx = 1_i32, nsrc
        call init_panel_geometry(v0(:, idx), v1(:, idx), v2(:, idx), plan%panel_geometry(idx), panel_status)
        if (panel_status /= panel_geometry_ok) error stop 'FMM panel source contains an invalid triangle.'
        original_centroid(:, idx) = plan%panel_geometry(idx)%centroid
      end do
    end if

    ! A finite image sum must keep the same representatives when neutral sources are added or
    ! a subset is evaluated. Use the box-based primary cell for every periodic plan.
    if (options%use_periodic2 .and. nsrc > 0_i32) call wrap_src_pos_to_primary_cell(plan)
    if (plan%panel_source) then
      do idx = 1_i32, nsrc
        shift = plan%src_pos(:, idx) - original_centroid(:, idx)
        ! A panel that the periodic wrap did not move keeps the geometry built above. Adding the zero
        ! shift would only turn a -0.0 vertex coordinate into +0.0, so such panels are rebuilt as before.
        if (all(shift == 0.0_dp) .and. .not. has_negative_zero(plan%panel_geometry(idx)%vertex)) cycle
        shifted_vertex = plan%panel_geometry(idx)%vertex + spread(shift, 2, 3)
        call init_panel_geometry( &
          shifted_vertex(:, 1), shifted_vertex(:, 2), shifted_vertex(:, 3), plan%panel_geometry(idx), panel_status &
          )
        if (panel_status /= panel_geometry_ok) error stop 'Periodic wrapping produced an invalid FMM panel source.'
      end do
      deallocate (original_centroid)
    end if
    plan%sources_ready = .true.
  end subroutine prepare_plan_sources

  pure logical function has_negative_zero(values)
    real(dp), intent(in) :: values(:, :)

    has_negative_zero = any(ieee_class(values) == ieee_negative_zero)
  end function has_negative_zero

  !> Build the tree, interaction lists and operators of a plan whose sources are prepared.
  subroutine core_complete_plan_impl(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: nsrc, idx, zero_status, free_axis
    real(dp), allocatable :: source_heights(:, :)
    character(len=256) :: zero_message

    nsrc = plan%nsrc
    if (use_periodic2_cached_kneq0(plan)) then
      free_axis = 6_i32 - sum(plan%options%periodic_axes)
      allocate (source_heights(3, nsrc))
      do idx = 1_i32, nsrc
        if (plan%panel_source) then
          source_heights(:, idx) = plan%panel_geometry(idx)%vertex(free_axis, :)
        else
          source_heights(:, idx) = plan%src_pos(free_axis, idx)
        end if
      end do
      call build_periodic_zero_mode_height_plan( &
        source_heights, product(plan%options%periodic_len), plan%periodic_k0_plan, zero_status, zero_message &
        )
      if (zero_status /= periodic_zero_mode_ok) error stop 'cached periodic k0 plan: '//trim(zero_message)
      plan%periodic_k0_free_axis = free_axis
      plan%periodic_k0_ready = .true.
      deallocate (source_heights)
    end if

    call initialize_basis_tables(plan, plan%options%order)
    if (use_periodic2_cached_kneq0(plan)) call precompute_periodic_source_grid(plan)
    call build_source_tree(plan)
    call precompute_source_p2m_basis(plan)
    call build_target_topology(plan)
    call build_interactions(plan)
    call precompute_translation_operators(plan)
    call precompute_periodic_shell_mean(plan)
    call precompute_periodic2_ewald_data(plan)
    call precompute_periodic_root_operator(plan)
    call precompute_m2l_derivatives(plan)
    plan%built = .true.
  end subroutine core_complete_plan_impl

  !> cached_kneq0 で引く、有限画像殻の面平均の滑らかな残りを線源ごとの高さ表にする。
  !! 平均をとるセルは線源の木の根の中心に置き、遠方補正の当てはめと同じ定義を使う。
  subroutine precompute_periodic_shell_mean(plan)
    type(fmm_plan_type), intent(inout) :: plan
    real(dp), parameter :: spacing_per_period = 32.0_dp
    real(dp), allocatable :: vertices(:, :, :)
    integer(i32) :: idx, axes(2), free_axis

    if (.not. use_periodic2_cached_kneq0(plan) .or. .not. plan%periodic_k0_ready .or. plan%nsrc <= 0_i32) return
    axes = plan%options%periodic_axes
    free_axis = plan%periodic_k0_free_axis
    allocate (vertices(3, 3, plan%nsrc))
    do idx = 1_i32, plan%nsrc
      if (plan%panel_source) then
        vertices(:, :, idx) = plan%panel_geometry(idx)%vertex
      else
        vertices(:, :, idx) = spread(plan%src_pos(:, idx), 2, 3)
      end if
    end do
    call build_periodic_shell_mean_plan( &
      axes, free_axis, 0.5_dp*(plan%options%target_box_min(axes) + plan%options%target_box_max(axes)), &
      plan%options%periodic_len, plan%options%periodic_image_layers, &
      plan%options%softening, plan%options%target_box_min(free_axis), plan%options%target_box_max(free_axis), &
      minval(plan%options%periodic_len)/spacing_per_period, vertices, plan%periodic_shell_mean &
      )
  end subroutine precompute_periodic_shell_mean

  !> FMM 計画に確保した資源を解放する。
  !! @param[inout] plan 解放対象の FMM 計画。
  subroutine core_destroy_plan_impl(plan)
    type(fmm_plan_type), intent(inout) :: plan

    call reset_fmm_plan(plan)
  end subroutine core_destroy_plan_impl

  subroutine build_source_tree(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: i, max_node_guess

    if (plan%nsrc <= 0_i32) return
    max_node_guess = max(1_i32, 2_i32*plan%nsrc)
    call ensure_source_tree_capacity(plan, max_node_guess)
    allocate (plan%source_tree%elem_order(plan%nsrc))
    do i = 1_i32, plan%nsrc
      plan%source_tree%elem_order(i) = i
    end do
    plan%source_tree%nnode = 1_i32
    plan%source_tree%max_depth = 0_i32
    call build_source_node(plan, 1_i32, 1_i32, plan%nsrc, 0_i32)
    call rebuild_level_cache(plan%source_tree)
    call build_parent_index(plan%source_tree)
    call build_source_leaf_index(plan)
  end subroutine build_source_tree

  recursive subroutine build_source_node(plan, node_idx, start_idx, end_idx, depth)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32), intent(in) :: node_idx, start_idx, end_idx, depth
    integer(i32) :: count, p, idx, oct
    integer(i32) :: child_k, child_node, child_start, child_end
    integer(i32) :: counts(8), offsets(8), cursor(8)
    integer(i32), allocatable :: work(:)
    real(dp) :: bb_min(3), bb_max(3), span(3), center(3), split_eps

    count = end_idx - start_idx + 1_i32
    plan%source_tree%node_start(node_idx) = start_idx
    plan%source_tree%node_count(node_idx) = count
    plan%source_tree%child_count(node_idx) = 0_i32
    plan%source_tree%child_idx(:, node_idx) = 0_i32
    plan%source_tree%child_octant(:, node_idx) = 0_i32
    plan%source_tree%node_depth(node_idx) = depth
    plan%source_tree%max_depth = max(plan%source_tree%max_depth, depth)

    idx = plan%source_tree%elem_order(start_idx)
    if (plan%panel_source) then
      bb_min = minval(plan%panel_geometry(idx)%vertex, dim=2)
      bb_max = maxval(plan%panel_geometry(idx)%vertex, dim=2)
    else
      bb_min = plan%src_pos(:, idx)
      bb_max = bb_min
    end if
    do p = start_idx + 1_i32, end_idx
      idx = plan%source_tree%elem_order(p)
      if (plan%panel_source) then
        bb_min = min(bb_min, minval(plan%panel_geometry(idx)%vertex, dim=2))
        bb_max = max(bb_max, maxval(plan%panel_geometry(idx)%vertex, dim=2))
      else
        bb_min = min(bb_min, plan%src_pos(:, idx))
        bb_max = max(bb_max, plan%src_pos(:, idx))
      end if
    end do

    span = bb_max - bb_min
    center = 0.5d0*(bb_max + bb_min)
    plan%source_tree%node_center(:, node_idx) = center
    plan%source_tree%node_half_size(:, node_idx) = 0.5d0*span
    plan%source_tree%node_radius(node_idx) = sqrt(sum(plan%source_tree%node_half_size(:, node_idx)* &
                                                      plan%source_tree%node_half_size(:, node_idx)))
    if (count <= plan%options%leaf_max) return

    split_eps = 1.0d-12*max(1.0d0, maxval(abs(center)))
    if (maxval(span) <= split_eps) return

    counts = 0_i32
    do p = start_idx, end_idx
      idx = plan%source_tree%elem_order(p)
      oct = octant_index(plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), center)
      counts(oct) = counts(oct) + 1_i32
    end do
    if (maxval(counts) == count) return
    allocate (work(count))

    offsets(1) = 1_i32
    do oct = 2, 8
      offsets(oct) = offsets(oct - 1) + counts(oct - 1)
    end do
    cursor = offsets
    do p = start_idx, end_idx
      idx = plan%source_tree%elem_order(p)
      oct = octant_index(plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), center)
      work(cursor(oct)) = idx
      cursor(oct) = cursor(oct) + 1_i32
    end do
    plan%source_tree%elem_order(start_idx:end_idx) = work
    deallocate (work)

    child_k = 0_i32
    child_start = start_idx
    do oct = 1_i32, 8_i32
      if (counts(oct) <= 0_i32) cycle
      child_end = child_start + counts(oct) - 1_i32
      plan%source_tree%nnode = plan%source_tree%nnode + 1_i32
      if (plan%source_tree%nnode > plan%source_tree%max_node) error stop 'FMM source tree capacity exceeded.'
      child_node = plan%source_tree%nnode
      child_k = child_k + 1_i32
      plan%source_tree%child_idx(child_k, node_idx) = child_node
      plan%source_tree%child_octant(child_k, node_idx) = oct
      call build_source_node(plan, child_node, child_start, child_end, depth + 1_i32)
      child_start = child_end + 1_i32
    end do
    plan%source_tree%child_count(node_idx) = child_k
  end subroutine build_source_node

  subroutine ensure_source_tree_capacity(plan, max_node_needed)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32), intent(in) :: max_node_needed

    call allocate_tree_nodes(plan%source_tree, max_node_needed)
    allocate (plan%source_tree%node_start(max_node_needed), plan%source_tree%node_count(max_node_needed))
    plan%source_tree%node_start = 0_i32
    plan%source_tree%node_count = 0_i32
  end subroutine ensure_source_tree_capacity

  !> Allocate the per-node arrays shared by the source and target trees for max_node nodes, zeroed.
  subroutine allocate_tree_nodes(tree, max_node)
    type(fmm_tree_type), intent(inout) :: tree
    integer(i32), intent(in) :: max_node

    tree%max_node = max_node
    allocate (tree%child_count(max_node), tree%child_idx(8, max_node), tree%child_octant(8, max_node))
    allocate (tree%node_depth(max_node))
    allocate (tree%node_center(3, max_node), tree%node_half_size(3, max_node), tree%node_radius(max_node))
    tree%child_count = 0_i32
    tree%child_idx = 0_i32
    tree%child_octant = 0_i32
    tree%node_depth = 0_i32
    tree%node_center = 0.0d0
    tree%node_half_size = 0.0d0
    tree%node_radius = 0.0d0
  end subroutine allocate_tree_nodes

  !> 深さ順のノード表 level_start/level_nodes を作る。各深さの中ではノード番号の昇順に並ぶ。
  subroutine rebuild_level_cache(tree)
    type(fmm_tree_type), intent(inout) :: tree
    integer(i32) :: node_idx, depth
    integer(i32), allocatable :: depth_count(:), cursor(:)

    if (tree%nnode <= 0_i32) return
    allocate (depth_count(tree%max_depth + 1_i32), cursor(tree%max_depth + 1_i32))
    depth_count = 0_i32
    do node_idx = 1_i32, tree%nnode
      depth_count(tree%node_depth(node_idx) + 1_i32) = depth_count(tree%node_depth(node_idx) + 1_i32) + 1_i32
    end do
    allocate (tree%level_start(tree%max_depth + 2_i32), tree%level_nodes(tree%nnode))
    tree%level_start(1) = 1_i32
    do depth = 0_i32, tree%max_depth
      tree%level_start(depth + 2_i32) = tree%level_start(depth + 1_i32) + depth_count(depth + 1_i32)
    end do
    cursor = tree%level_start(1:tree%max_depth + 1_i32)
    do node_idx = 1_i32, tree%nnode
      depth = tree%node_depth(node_idx)
      tree%level_nodes(cursor(depth + 1_i32)) = node_idx
      cursor(depth + 1_i32) = cursor(depth + 1_i32) + 1_i32
    end do
    deallocate (depth_count, cursor)
  end subroutine rebuild_level_cache

  subroutine build_parent_index(tree)
    type(fmm_tree_type), intent(inout) :: tree
    integer(i32) :: parent_node, child_k, child_node

    if (allocated(tree%parent)) deallocate (tree%parent)
    allocate (tree%parent(max(1_i32, tree%nnode)))
    tree%parent = 0_i32
    do parent_node = 1_i32, tree%nnode
      do child_k = 1_i32, tree%child_count(parent_node)
        child_node = tree%child_idx(child_k, parent_node)
        tree%parent(child_node) = parent_node
      end do
    end do
  end subroutine build_parent_index

  subroutine build_source_leaf_index(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: node_idx, leaf_idx

    if (allocated(plan%source_leaf_nodes)) deallocate (plan%source_leaf_nodes)
    plan%nsource_leaf = 0_i32
    do node_idx = 1_i32, plan%source_tree%nnode
      if (plan%source_tree%child_count(node_idx) <= 0_i32) plan%nsource_leaf = plan%nsource_leaf + 1_i32
    end do
    if (plan%nsource_leaf <= 0_i32) return

    allocate (plan%source_leaf_nodes(plan%nsource_leaf))
    leaf_idx = 0_i32
    do node_idx = 1_i32, plan%source_tree%nnode
      if (plan%source_tree%child_count(node_idx) > 0_i32) cycle
      leaf_idx = leaf_idx + 1_i32
      plan%source_leaf_nodes(leaf_idx) = node_idx
    end do
  end subroutine build_source_leaf_index

  !> The far correction uses all box cells, including empty cells. The adaptive near tree is independent.
  !! Cells and anchors form box lattices with per-axis counts chosen so that every cell is nearly cubic
  !! and its half diagonal stays below the radius ratio times the shorter period; a tall box therefore
  !! adds cells along its long axis instead of refining the lateral subdivision.
  subroutine precompute_periodic_source_grid(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: cell, ix, iy, iz, source, index3(3), alpha_idx, axis, ncell(3), nanchor(3)
    real(dp) :: span(3), half(3), proxy_half(3), center(3), d(3), tolerance, period
    real(dp) :: xpow(0:plan%options%order), ypow(0:plan%options%order), zpow(0:plan%options%order)

    span = plan%options%target_box_max - plan%options%target_box_min
    period = minval(plan%options%periodic_len)
    nanchor = plan%options%periodic_operator_anchor_cells
    if (any(nanchor < 1_i32)) then
      nanchor = box_operator_cells(periodic_operator_target_radius_ratio)
      ! Two anchors across each period keep every anchor on at most one periodic seam per axis.
      nanchor(plan%options%periodic_axes) = max(nanchor(plan%options%periodic_axes), 2_i32)
    end if
    ncell = plan%options%periodic_operator_source_cells
    if (any(ncell < 1_i32)) then
      ! Source cells refine the anchor lattice by powers of two, so every anchor-cell offset is a multiple of
      ! half the source spacing and the generator can share Ewald residuals between translation classes.
      ncell = nanchor
      do while (norm2(0.5_dp*span/real(ncell, dp)) > periodic_operator_source_radius_ratio*period)
        axis = int(maxloc(span/real(ncell, dp), dim=1), i32)
        ncell(axis) = 2_i32*ncell(axis)
      end do
    end if
    plan%options%periodic_operator_anchor_cells = nanchor
    plan%options%periodic_operator_source_cells = ncell
    half = 0.5_dp*span/real(ncell, dp)
    plan%periodic_operator%source_half = half
    plan%periodic_operator%anchor_half = 0.5_dp*span/real(nanchor, dp)
    proxy_half = periodic_source_proxy_half_size(half)
    plan%periodic_operator%root_source_count = product(ncell)
    allocate (plan%periodic_operator%source_center(3, plan%periodic_operator%root_source_count), &
              plan%periodic_operator%source_cell(plan%nsrc))
    allocate (plan%periodic_operator%source_p2m_basis(plan%ncoef, plan%nsrc))
    do iz = 0_i32, ncell(3) - 1_i32
      do iy = 0_i32, ncell(2) - 1_i32
        do ix = 0_i32, ncell(1) - 1_i32
          cell = 1_i32 + ix + ncell(1)*(iy + ncell(2)*iz)
          plan%periodic_operator%source_center(:, cell) = plan%options%target_box_min + &
                                                          (real([ix, iy, iz], dp) + 0.5_dp)*2.0_dp*half
        end do
      end do
    end do
    plan%periodic_operator%anchor_count = product(nanchor)
    allocate (plan%periodic_operator%anchor_center(3, plan%periodic_operator%anchor_count))
    do iz = 0_i32, nanchor(3) - 1_i32
      do iy = 0_i32, nanchor(2) - 1_i32
        do ix = 0_i32, nanchor(1) - 1_i32
          cell = 1_i32 + ix + nanchor(1)*(iy + nanchor(2)*iz)
          plan%periodic_operator%anchor_center(:, cell) = &
            plan%options%target_box_min + (real([ix, iy, iz], dp) + 0.5_dp)*2.0_dp*plan%periodic_operator%anchor_half
        end do
      end do
    end do
    tolerance = 64.0_dp*epsilon(1.0_dp)*maxval(span)
    do source = 1_i32, plan%nsrc
      ! Primary-cell centroids on the upper free face belong to the last cell.
      index3 = min(ncell - 1_i32, max(0_i32, floor((plan%src_pos(:, source) - &
                                                    plan%options%target_box_min)/span*real(ncell, dp), kind=i32)))
      cell = 1_i32 + index3(1) + ncell(1)*(index3(2) + ncell(2)*index3(3))
      plan%periodic_operator%source_cell(source) = cell
      center = plan%periodic_operator%source_center(:, cell)
      if (plan%panel_source) then
        if (any(abs(plan%panel_geometry(source)%vertex - spread(center, 2, 3)) > &
                spread(proxy_half + tolerance, 2, 3))) then
          error stop 'cached_kneq0: panel exceeds its box-cell proxy shell (scale 1.5). '// &
            'Refine the mesh or decrease periodic_operator_source_cells; keep panels inside the source box.'
        end if
        do alpha_idx = 1_i32, plan%ncoef
          plan%periodic_operator%source_p2m_basis(alpha_idx, source) = &
            panel_monomial_average(plan%panel_geometry(source), center, plan%alpha(:, alpha_idx))/ &
            plan%alpha_factorial(alpha_idx)
        end do
      else
        d = plan%src_pos(:, source) - center
        if (any(abs(d) > proxy_half + tolerance)) then
          error stop 'cached_kneq0: point source is outside the box-cell proxy shell; enlarge the target box.'
        end if
        call build_axis_powers(d, plan%options%order, xpow, ypow, zpow)
        do alpha_idx = 1_i32, plan%ncoef
          plan%periodic_operator%source_p2m_basis(alpha_idx, source) = &
            xpow(plan%alpha(1, alpha_idx))*ypow(plan%alpha(2, alpha_idx))*zpow(plan%alpha(3, alpha_idx))/ &
            plan%alpha_factorial(alpha_idx)
        end do
      end if
    end do

  contains

    !> Smallest per-axis counts whose cell half diagonal is at most ratio times the shorter period,
    !! always splitting the axis with the widest cell so the cells stay nearly cubic.
    function box_operator_cells(ratio) result(cells)
      real(dp), intent(in) :: ratio
      integer(i32) :: cells(3), widest

      cells = 1_i32
      do while (norm2(0.5_dp*span/real(cells, dp)) > ratio*period)
        widest = int(maxloc(span/real(cells, dp), dim=1), i32)
        cells(widest) = cells(widest) + 1_i32
      end do
    end function box_operator_cells
  end subroutine precompute_periodic_source_grid

  subroutine precompute_source_p2m_basis(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: leaf_idx, node_idx, p, idx, p_end, alpha_idx
    real(dp) :: d(3)
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    if (allocated(plan%source_p2m_basis)) deallocate (plan%source_p2m_basis)
    if (plan%nsrc <= 0_i32) return

    allocate (plan%source_p2m_basis(plan%ncoef, plan%nsrc))
    plan%source_p2m_basis = 0.0d0
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(plan) private(leaf_idx, node_idx, p, idx, p_end, alpha_idx, d, xpow, ypow, zpow)
    do leaf_idx = 1_i32, plan%nsource_leaf
      node_idx = plan%source_leaf_nodes(leaf_idx)
      p_end = plan%source_tree%node_start(node_idx) + plan%source_tree%node_count(node_idx) - 1_i32
      do p = plan%source_tree%node_start(node_idx), p_end
        idx = plan%source_tree%elem_order(p)
        if (plan%panel_source) then
          do alpha_idx = 1_i32, plan%ncoef
            plan%source_p2m_basis(alpha_idx, p) = panel_monomial_average( &
                                                  plan%panel_geometry(idx), plan%source_tree%node_center(:, node_idx), &
                                                  plan%alpha(:, alpha_idx))/plan%alpha_factorial(alpha_idx)
          end do
        else
          d = plan%src_pos(:, idx) - plan%source_tree%node_center(:, node_idx)
          call build_axis_powers(d, plan%options%order, xpow, ypow, zpow)
          do alpha_idx = 1_i32, plan%ncoef
            plan%source_p2m_basis(alpha_idx, p) = xpow(plan%alpha(1, alpha_idx))*ypow(plan%alpha(2, alpha_idx)) &
                                                  *zpow(plan%alpha(3, alpha_idx))/plan%alpha_factorial(alpha_idx)
          end do
        end if
      end do
    end do
    !$omp end parallel do
  end subroutine precompute_source_p2m_basis

  pure real(dp) function panel_monomial_average(geometry, center, alpha) result(moment)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: center(3)
    integer(i32), intent(in) :: alpha(3)
    integer(i32) :: ix0, ix1, ix2, iy0, iy1, iy2, iz0, iz1, iz2
    integer(i32) :: lambda_power(3), degree
    real(dp) :: dv(3, 3), coefficient

    dv = geometry%vertex - spread(center, 2, 3)
    moment = 0.0_dp
    degree = sum(alpha)
    do ix0 = 0_i32, alpha(1)
      do ix1 = 0_i32, alpha(1) - ix0
        ix2 = alpha(1) - ix0 - ix1
        do iy0 = 0_i32, alpha(2)
          do iy1 = 0_i32, alpha(2) - iy0
            iy2 = alpha(2) - iy0 - iy1
            do iz0 = 0_i32, alpha(3)
              do iz1 = 0_i32, alpha(3) - iz0
                iz2 = alpha(3) - iz0 - iz1
                lambda_power = [ix0 + iy0 + iz0, ix1 + iy1 + iz1, ix2 + iy2 + iz2]
                coefficient = multinomial3(alpha(1), ix0, ix1, ix2)* &
                              multinomial3(alpha(2), iy0, iy1, iy2)* &
                              multinomial3(alpha(3), iz0, iz1, iz2)
                coefficient = coefficient*dv(1, 1)**ix0*dv(1, 2)**ix1*dv(1, 3)**ix2
                coefficient = coefficient*dv(2, 1)**iy0*dv(2, 2)**iy1*dv(2, 3)**iy2
                coefficient = coefficient*dv(3, 1)**iz0*dv(3, 2)**iz1*dv(3, 3)**iz2
                moment = moment + coefficient*2.0_dp*factorial_i(lambda_power(1))* &
                         factorial_i(lambda_power(2))*factorial_i(lambda_power(3))/factorial_i(degree + 2_i32)
              end do
            end do
          end do
        end do
      end do
    end do
  end function panel_monomial_average

  pure real(dp) function multinomial3(n, a, b, c) result(value)
    integer(i32), intent(in) :: n, a, b, c
    value = factorial_i(n)/(factorial_i(a)*factorial_i(b)*factorial_i(c))
  end function multinomial3

  pure real(dp) function factorial_i(n) result(value)
    integer(i32), intent(in) :: n
    integer(i32) :: k
    value = 1.0_dp
    do k = 2_i32, n
      value = value*real(k, dp)
    end do
  end function factorial_i

  subroutine build_target_topology(plan)
    type(fmm_plan_type), intent(inout) :: plan

    if (plan%source_tree%nnode <= 0_i32) then
      call clear_target_tree(plan)
      return
    end if

    if (has_valid_target_box(plan%options)) then
      call build_box_target_tree(plan, plan%options%use_periodic2)
    else
      call assign_source_leaves_as_targets(plan)
    end if
  end subroutine build_target_topology

  subroutine clear_target_tree(plan)
    type(fmm_plan_type), intent(inout) :: plan

    plan%target_tree = fmm_tree_type()
    if (allocated(plan%leaf_nodes)) deallocate (plan%leaf_nodes)
    if (allocated(plan%leaf_slot_of_node)) deallocate (plan%leaf_slot_of_node)
    plan%target_tree_ready = .false.
    plan%nleaf = 0_i32
  end subroutine clear_target_tree

  subroutine assign_source_leaves_as_targets(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: node_idx, leaf_idx

    call clear_target_tree(plan)
    do node_idx = 1_i32, plan%source_tree%nnode
      if (plan%source_tree%child_count(node_idx) <= 0_i32) plan%nleaf = plan%nleaf + 1_i32
    end do
    if (plan%nleaf <= 0_i32) return
    allocate (plan%leaf_nodes(plan%nleaf), plan%leaf_slot_of_node(plan%source_tree%max_node))
    plan%leaf_slot_of_node = 0_i32
    leaf_idx = 0_i32
    do node_idx = 1_i32, plan%source_tree%nnode
      if (plan%source_tree%child_count(node_idx) > 0_i32) cycle
      leaf_idx = leaf_idx + 1_i32
      plan%leaf_nodes(leaf_idx) = node_idx
      plan%leaf_slot_of_node(node_idx) = leaf_idx
    end do
  end subroutine assign_source_leaves_as_targets

  subroutine build_box_target_tree(plan, use_periodic_wrap)
    type(fmm_plan_type), intent(inout) :: plan
    logical, intent(in) :: use_periodic_wrap
    integer(i32) :: depth_min, depth_max
    integer(i32) :: leaf_idx, node_idx, target_cap, level_nodes, depth
    real(dp) :: root_min(3), root_max(3), src_center(3), src_half(3)

    if (plan%source_tree%nnode <= 0_i32) then
      call clear_target_tree(plan)
      return
    end if

    src_center = plan%source_tree%node_center(:, 1_i32)
    src_half = plan%source_tree%node_half_size(:, 1_i32)
    root_min = plan%options%target_box_min
    root_max = plan%options%target_box_max
    if (any(root_max <= root_min)) then
      root_min = src_center - src_half
      root_max = src_center + src_half
    end if

    depth_min = 2_i32
    depth_max = 4_i32
    target_cap = 0_i32
    level_nodes = 1_i32
    do depth = 0_i32, depth_max
      target_cap = target_cap + level_nodes
      if (depth < depth_max) then
        if (level_nodes > huge(level_nodes)/8_i32) error stop 'FMM target tree capacity overflow.'
        level_nodes = 8_i32*level_nodes
      end if
    end do

    call clear_target_tree(plan)
    ! depth_max bounds the complete target octree independently of the source tree.
    call init_box_target_root(plan%target_tree, target_cap, root_min, root_max)
    call build_target_box_node(plan, 1_i32, 0_i32, depth_min, depth_max, src_center, src_half, use_periodic_wrap)

    do node_idx = 1_i32, plan%target_tree%nnode
      if (plan%target_tree%child_count(node_idx) <= 0_i32) plan%nleaf = plan%nleaf + 1_i32
    end do
    if (plan%nleaf <= 0_i32) return

    allocate (plan%leaf_nodes(plan%nleaf), plan%leaf_slot_of_node(plan%target_tree%max_node))
    plan%leaf_slot_of_node = 0_i32
    leaf_idx = 0_i32
    do node_idx = 1_i32, plan%target_tree%nnode
      if (plan%target_tree%child_count(node_idx) > 0_i32) cycle
      leaf_idx = leaf_idx + 1_i32
      plan%leaf_nodes(leaf_idx) = node_idx
      plan%leaf_slot_of_node(node_idx) = leaf_idx
    end do

    call rebuild_level_cache(plan%target_tree)
    call build_parent_index(plan%target_tree)
    plan%target_tree_ready = .true.
  end subroutine build_box_target_tree

  !> Allocate a box target tree of capacity max_node with only its root, which covers [root_min, root_max].
  subroutine init_box_target_root(tree, max_node, root_min, root_max)
    type(fmm_tree_type), intent(inout) :: tree
    integer(i32), intent(in) :: max_node
    real(dp), intent(in) :: root_min(3), root_max(3)

    call allocate_tree_nodes(tree, max_node)
    tree%nnode = 1_i32
    tree%node_depth(1_i32) = 0_i32
    tree%node_center(:, 1_i32) = 0.5d0*(root_min + root_max)
    tree%node_half_size(:, 1_i32) = 0.5d0*(root_max - root_min)
    tree%node_radius(1_i32) = sqrt(sum(tree%node_half_size(:, 1_i32)*tree%node_half_size(:, 1_i32)))
  end subroutine init_box_target_root

  recursive subroutine build_target_box_node(plan, node_idx, depth, depth_min, depth_max, src_center, src_half, use_periodic_wrap)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32), intent(in) :: node_idx, depth, depth_min, depth_max
    real(dp), intent(in) :: src_center(3), src_half(3)
    logical, intent(in) :: use_periodic_wrap
    integer(i32) :: child_k, child_node, oct
    real(dp) :: child_half(3), offset(3)

    if (.not. target_should_split_box( &
        plan, node_idx, depth, depth_min, depth_max, src_center, src_half, use_periodic_wrap)) return
    child_half = 0.5d0*plan%target_tree%node_half_size(:, node_idx)
    if (maxval(child_half) <= tiny(1.0d0)) return

    child_k = 0_i32
    do oct = 1_i32, 8_i32
      plan%target_tree%nnode = plan%target_tree%nnode + 1_i32
      if (plan%target_tree%nnode > plan%target_tree%max_node) error stop 'FMM target tree capacity exceeded.'
      child_node = plan%target_tree%nnode
      child_k = child_k + 1_i32
      plan%target_tree%child_idx(child_k, node_idx) = child_node
      plan%target_tree%child_octant(child_k, node_idx) = oct
      plan%target_tree%node_depth(child_node) = depth + 1_i32
      plan%target_tree%max_depth = max(plan%target_tree%max_depth, depth + 1_i32)

      offset(1) = merge(-1.0d0, 1.0d0, iand(oct - 1_i32, 1_i32) == 0_i32)
      offset(2) = merge(-1.0d0, 1.0d0, iand(oct - 1_i32, 2_i32) == 0_i32)
      offset(3) = merge(-1.0d0, 1.0d0, iand(oct - 1_i32, 4_i32) == 0_i32)
      plan%target_tree%node_half_size(:, child_node) = child_half
      plan%target_tree%node_center(:, child_node) = plan%target_tree%node_center(:, node_idx) + offset*child_half
      plan%target_tree%node_radius(child_node) = sqrt(sum(child_half*child_half))
      call build_target_box_node( &
        plan, child_node, depth + 1_i32, depth_min, depth_max, src_center, src_half, use_periodic_wrap &
        )
    end do
    plan%target_tree%child_count(node_idx) = child_k
  end subroutine build_target_box_node

  logical function target_should_split_box(plan, node_idx, depth, depth_min, depth_max, src_center, src_half, use_periodic_wrap)
    type(fmm_plan_type), intent(in) :: plan
    integer(i32), intent(in) :: node_idx, depth, depth_min, depth_max
    real(dp), intent(in) :: src_center(3), src_half(3)
    logical, intent(in) :: use_periodic_wrap
    real(dp) :: dsrc, radius

    if (depth < depth_min) then
      target_should_split_box = .true.
      return
    end if
    if (depth >= depth_max) then
      target_should_split_box = .false.
      return
    end if

    radius = plan%target_tree%node_radius(node_idx)
    if (use_periodic_wrap) then
      dsrc = distance_to_source_bbox_periodic(plan, plan%target_tree%node_center(:, node_idx), src_center, src_half)
    else
      dsrc = distance_to_source_bbox(plan%target_tree%node_center(:, node_idx), src_center, src_half)
    end if
    target_should_split_box = (radius > 0.35d0*max(dsrc, 1.0d-12))
  end function target_should_split_box

  subroutine build_interactions(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: nshift

    nshift = 1_i32
    if (plan%options%use_periodic2) nshift = 2_i32*max(0_i32, plan%options%periodic_image_layers) + 1_i32
    allocate (plan%shift_axis1(nshift), plan%shift_axis2(nshift))
    call build_periodic_shift_values(plan, plan%shift_axis1, plan%shift_axis2, nshift)

    if (plan%nleaf <= 0_i32) return
    if (plan%target_tree_ready) then
      call build_interaction_lists(plan, plan%target_tree)
    else
      call build_interaction_lists(plan, plan%source_tree)
    end if
  end subroutine build_interactions

  !> Near lists, M2L pairs and near-field tables of every target leaf.
  !! @param[in] target_tree 評価用の木（plan%target_tree か plan%source_tree）。読むだけで、plan の他の成分を書く。
  subroutine build_interaction_lists(plan, target_tree)
    type(fmm_plan_type), intent(inout) :: plan
    type(fmm_tree_type), intent(in) :: target_tree
    integer(i32) :: leaf_idx, i, near_n, near_used, near_cap
    integer(i32) :: n_target_nodes
    integer(i32) :: pair_cap, pair_used
    integer(i32), allocatable :: near_work(:), near_buf(:)
    integer(i32), allocatable :: m2l_target_buf(:), m2l_source_buf(:), m2l_shift_idx1_buf(:), m2l_shift_idx2_buf(:)

    n_target_nodes = target_tree%nnode
    allocate (plan%near_start(plan%nleaf + 1_i32))
    allocate (near_work(plan%source_tree%nnode))
    near_cap = max(64_i32, plan%nleaf*16_i32)
    allocate (near_buf(near_cap))
    near_used = 0_i32
    do leaf_idx = 1_i32, plan%nleaf
      near_n = 0_i32
      call gather_leaf_near_sources(plan, target_tree, plan%leaf_nodes(leaf_idx), 1_i32, near_work, near_n)
      plan%near_start(leaf_idx) = near_used + 1_i32
      do i = 1_i32, near_n
        call append_i32_buffer(near_buf, near_used, near_cap, near_work(i))
      end do
    end do
    plan%near_start(plan%nleaf + 1_i32) = near_used + 1_i32
    allocate (plan%near_nodes(max(1_i32, near_used)))
    plan%near_nodes = 0_i32
    if (near_used > 0_i32) plan%near_nodes(1:near_used) = near_buf(1:near_used)

    pair_cap = max(64_i32, plan%nleaf*32_i32)
    allocate (m2l_target_buf(pair_cap), m2l_source_buf(pair_cap), m2l_shift_idx1_buf(pair_cap), m2l_shift_idx2_buf(pair_cap))
    m2l_target_buf = 0_i32
    m2l_source_buf = 0_i32
    m2l_shift_idx1_buf = 0_i32
    m2l_shift_idx2_buf = 0_i32
    pair_used = 0_i32
    call accumulate_cached_pairs( &
      plan, target_tree, 1_i32, 1_i32, m2l_target_buf, m2l_source_buf, m2l_shift_idx1_buf, m2l_shift_idx2_buf, &
      pair_used, pair_cap &
      )
    call build_near_source_index( &
      plan, target_tree, m2l_target_buf, m2l_source_buf, m2l_shift_idx1_buf, m2l_shift_idx2_buf, pair_used, pair_cap &
      )
    if (plan%panel_source) call build_near_panel_tables(plan, target_tree)
    allocate ( &
      plan%m2l_target_nodes(max(1_i32, pair_used)), plan%m2l_source_nodes(max(1_i32, pair_used)), &
      plan%m2l_shift_idx1(max(1_i32, pair_used)), plan%m2l_shift_idx2(max(1_i32, pair_used)) &
      )
    plan%m2l_target_nodes = 0_i32
    plan%m2l_source_nodes = 0_i32
    plan%m2l_shift_idx1 = 0_i32
    plan%m2l_shift_idx2 = 0_i32
    if (pair_used > 0_i32) then
      plan%m2l_target_nodes(1:pair_used) = m2l_target_buf(1:pair_used)
      plan%m2l_source_nodes(1:pair_used) = m2l_source_buf(1:pair_used)
      plan%m2l_shift_idx1(1:pair_used) = m2l_shift_idx1_buf(1:pair_used)
      plan%m2l_shift_idx2(1:pair_used) = m2l_shift_idx2_buf(1:pair_used)
    end if
    plan%m2l_pair_count = pair_used
    call sort_pairs_by_target(plan, n_target_nodes)
    call group_targets_for_lanes(plan, n_target_nodes)

    deallocate (near_work, near_buf, m2l_target_buf, m2l_source_buf, m2l_shift_idx1_buf, m2l_shift_idx2_buf)
  end subroutine build_interaction_lists

  subroutine build_near_source_index( &
    plan, target_tree, m2l_target_buf, m2l_source_buf, m2l_shift_idx1_buf, m2l_shift_idx2_buf, pair_used, pair_cap &
    )
    type(fmm_plan_type), intent(inout) :: plan
    type(fmm_tree_type), intent(in) :: target_tree
    integer(i32), allocatable, intent(inout) :: m2l_target_buf(:), m2l_source_buf(:), &
                                                m2l_shift_idx1_buf(:), m2l_shift_idx2_buf(:)
    integer(i32), intent(inout) :: pair_used, pair_cap
    integer(i32) :: leaf_idx, near_pos, near_node, p, p_end, near_source_used
    integer(i32) :: shift_idx1, shift_idx2, nshift, target_leaf
    real(dp) :: shift1, shift2

    if (allocated(plan%near_source_start)) deallocate (plan%near_source_start)
    if (allocated(plan%near_source_idx)) deallocate (plan%near_source_idx)
    if (allocated(plan%near_source_shift1)) deallocate (plan%near_source_shift1)
    if (allocated(plan%near_source_shift2)) deallocate (plan%near_source_shift2)
    if (plan%nleaf <= 0_i32) return

    allocate (plan%near_source_start(plan%nleaf + 1_i32))
    plan%near_source_start = 1_i32
    near_source_used = 0_i32
    nshift = max(1_i32, int(size(plan%shift_axis1), i32))
    do leaf_idx = 1_i32, plan%nleaf
      plan%near_source_start(leaf_idx) = near_source_used + 1_i32
      target_leaf = plan%leaf_nodes(leaf_idx)
      do near_pos = plan%near_start(leaf_idx), plan%near_start(leaf_idx + 1_i32) - 1_i32
        near_node = plan%near_nodes(near_pos)
        if (plan%options%use_periodic2) then
          do shift_idx1 = 1_i32, nshift
            shift1 = plan%shift_axis1(shift_idx1)
            do shift_idx2 = 1_i32, nshift
              shift2 = plan%shift_axis2(shift_idx2)
              if (.not. nodes_well_separated_shifted(plan, target_tree, target_leaf, near_node, shift1, shift2)) then
                near_source_used = near_source_used + plan%source_tree%node_count(near_node)
              end if
            end do
          end do
        else
          near_source_used = near_source_used + plan%source_tree%node_count(near_node)
        end if
      end do
    end do
    plan%near_source_start(plan%nleaf + 1_i32) = near_source_used + 1_i32

    allocate ( &
      plan%near_source_idx(max(1_i32, near_source_used)), &
      plan%near_source_shift1(max(1_i32, near_source_used)), &
      plan%near_source_shift2(max(1_i32, near_source_used)) &
      )
    plan%near_source_idx = 0_i32
    plan%near_source_shift1 = 0.0d0
    plan%near_source_shift2 = 0.0d0
    if (near_source_used <= 0_i32) return

    near_source_used = 0_i32
    do leaf_idx = 1_i32, plan%nleaf
      target_leaf = plan%leaf_nodes(leaf_idx)
      do near_pos = plan%near_start(leaf_idx), plan%near_start(leaf_idx + 1_i32) - 1_i32
        near_node = plan%near_nodes(near_pos)
        p_end = plan%source_tree%node_start(near_node) + plan%source_tree%node_count(near_node) - 1_i32
        if (plan%options%use_periodic2) then
          do shift_idx1 = 1_i32, nshift
            shift1 = plan%shift_axis1(shift_idx1)
            do shift_idx2 = 1_i32, nshift
              shift2 = plan%shift_axis2(shift_idx2)
              if (nodes_well_separated_shifted(plan, target_tree, target_leaf, near_node, shift1, shift2)) then
                call append_m2l_pair( &
                  m2l_target_buf, m2l_source_buf, m2l_shift_idx1_buf, m2l_shift_idx2_buf, pair_used, pair_cap, &
                  target_leaf, near_node, shift_idx1, shift_idx2 &
                  )
                cycle
              end if
              do p = plan%source_tree%node_start(near_node), p_end
                near_source_used = near_source_used + 1_i32
                plan%near_source_idx(near_source_used) = plan%source_tree%elem_order(p)
                plan%near_source_shift1(near_source_used) = shift1
                plan%near_source_shift2(near_source_used) = shift2
              end do
            end do
          end do
        else
          do p = plan%source_tree%node_start(near_node), p_end
            near_source_used = near_source_used + 1_i32
            plan%near_source_idx(near_source_used) = plan%source_tree%elem_order(p)
          end do
        end if
      end do
    end do
  end subroutine build_near_source_index

  !> Store every near-field panel of every leaf, translated by its image shift, as component arrays
  !! in `near_source_idx` order. The translation is the same `vertex + shift` the point kernel used
  !! to apply per evaluation, so the tables reproduce it bit for bit.
  subroutine build_near_panel_tables(plan, target_tree)
    type(fmm_plan_type), intent(inout) :: plan
    type(fmm_tree_type), intent(in) :: target_tree
    integer(i32) :: n, pos, idx, axes(2), leaf_idx
    real(dp) :: shift(3), center(3)

    n = 0_i32
    if (plan%nleaf > 0_i32) n = plan%near_source_start(plan%nleaf + 1_i32) - 1_i32
    n = max(0_i32, n)
    axes = 0_i32
    if (plan%options%use_periodic2) axes = plan%options%periodic_axes
    if (plan%options%near_field_single) then
      call build_near_panel_tables_single(plan, target_tree, n, axes)
      return
    end if
    allocate (plan%near_panels%vx(3, max(1_i32, n)), plan%near_panels%vy(3, max(1_i32, n)), &
              plan%near_panels%vz(3, max(1_i32, n)))
    allocate (plan%near_panels%nx(max(1_i32, n)), plan%near_panels%ny(max(1_i32, n)), plan%near_panels%nz(max(1_i32, n)))
    allocate (plan%near_panels%edge(3, max(1_i32, n)), plan%near_panels%ox(3, max(1_i32, n)), &
              plan%near_panels%oy(3, max(1_i32, n)))
    allocate (plan%near_panels%oz(3, max(1_i32, n)), plan%near_panels%area(max(1_i32, n)))
    plan%near_panels%vx = 0.0_dp
    plan%near_panels%vy = 0.0_dp
    plan%near_panels%vz = 0.0_dp
    plan%near_panels%nx = 0.0_dp
    plan%near_panels%ny = 0.0_dp
    plan%near_panels%nz = 0.0_dp
    plan%near_panels%edge = 0.0_dp
    plan%near_panels%ox = 0.0_dp
    plan%near_panels%oy = 0.0_dp
    plan%near_panels%oz = 0.0_dp
    plan%near_panels%area = 1.0_dp
    do pos = 1_i32, n
      idx = plan%near_source_idx(pos)
      shift = 0.0_dp
      if (axes(1) > 0_i32) shift(axes(1)) = shift(axes(1)) + plan%near_source_shift1(pos)
      if (axes(2) > 0_i32) shift(axes(2)) = shift(axes(2)) + plan%near_source_shift2(pos)
      plan%near_panels%vx(:, pos) = plan%panel_geometry(idx)%vertex(1, :) + shift(1)
      plan%near_panels%vy(:, pos) = plan%panel_geometry(idx)%vertex(2, :) + shift(2)
      plan%near_panels%vz(:, pos) = plan%panel_geometry(idx)%vertex(3, :) + shift(3)
      plan%near_panels%nx(pos) = plan%panel_geometry(idx)%normal(1)
      plan%near_panels%ny(pos) = plan%panel_geometry(idx)%normal(2)
      plan%near_panels%nz(pos) = plan%panel_geometry(idx)%normal(3)
      plan%near_panels%edge(:, pos) = plan%panel_geometry(idx)%edge_length
      plan%near_panels%ox(:, pos) = plan%panel_geometry(idx)%edge_outward(1, :)
      plan%near_panels%oy(:, pos) = plan%panel_geometry(idx)%edge_outward(2, :)
      plan%near_panels%oz(:, pos) = plan%panel_geometry(idx)%edge_outward(3, :)
      plan%near_panels%area(pos) = plan%panel_geometry(idx)%area
    end do
  contains

    !> Single-precision tables: each panel is stored relative to the center of the leaf whose near list
    !! holds it, so the rounding error stays at single-precision resolution of the leaf size rather than
    !! of the absolute coordinates. The leaf center is kept in double for shifting the target.
    subroutine build_near_panel_tables_single(plan, target_tree, n, axes)
      type(fmm_plan_type), intent(inout) :: plan
      type(fmm_tree_type), intent(in) :: target_tree
      integer(i32), intent(in) :: n, axes(2)

      allocate (plan%near_panels_sp%vx(3, max(1_i32, n)), plan%near_panels_sp%vy(3, max(1_i32, n)), &
                plan%near_panels_sp%vz(3, max(1_i32, n)))
      allocate (plan%near_panels_sp%nx(max(1_i32, n)), plan%near_panels_sp%ny(max(1_i32, n)), &
                plan%near_panels_sp%nz(max(1_i32, n)))
      allocate (plan%near_panels_sp%edge(3, max(1_i32, n)), plan%near_panels_sp%ox(3, max(1_i32, n)), &
                plan%near_panels_sp%oy(3, max(1_i32, n)))
      allocate (plan%near_panels_sp%oz(3, max(1_i32, n)), plan%near_panels_sp%area(max(1_i32, n)))
      allocate (plan%near_panels_sp%leaf_center(3, max(1_i32, plan%nleaf)))
      plan%near_panels_sp%vx = 0.0_sp
      plan%near_panels_sp%vy = 0.0_sp
      plan%near_panels_sp%vz = 0.0_sp
      plan%near_panels_sp%nx = 0.0_sp
      plan%near_panels_sp%ny = 0.0_sp
      plan%near_panels_sp%nz = 0.0_sp
      plan%near_panels_sp%edge = 0.0_sp
      plan%near_panels_sp%ox = 0.0_sp
      plan%near_panels_sp%oy = 0.0_sp
      plan%near_panels_sp%oz = 0.0_sp
      plan%near_panels_sp%area = 1.0_sp
      plan%near_panels_sp%leaf_center = 0.0_dp
      do leaf_idx = 1_i32, plan%nleaf
        center = target_tree%node_center(:, plan%leaf_nodes(leaf_idx))
        plan%near_panels_sp%leaf_center(:, leaf_idx) = center
        do pos = plan%near_source_start(leaf_idx), plan%near_source_start(leaf_idx + 1_i32) - 1_i32
          idx = plan%near_source_idx(pos)
          shift = 0.0_dp
          if (axes(1) > 0_i32) shift(axes(1)) = shift(axes(1)) + plan%near_source_shift1(pos)
          if (axes(2) > 0_i32) shift(axes(2)) = shift(axes(2)) + plan%near_source_shift2(pos)
          plan%near_panels_sp%vx(:, pos) = real(plan%panel_geometry(idx)%vertex(1, :) + shift(1) - center(1), sp)
          plan%near_panels_sp%vy(:, pos) = real(plan%panel_geometry(idx)%vertex(2, :) + shift(2) - center(2), sp)
          plan%near_panels_sp%vz(:, pos) = real(plan%panel_geometry(idx)%vertex(3, :) + shift(3) - center(3), sp)
          plan%near_panels_sp%nx(pos) = real(plan%panel_geometry(idx)%normal(1), sp)
          plan%near_panels_sp%ny(pos) = real(plan%panel_geometry(idx)%normal(2), sp)
          plan%near_panels_sp%nz(pos) = real(plan%panel_geometry(idx)%normal(3), sp)
          plan%near_panels_sp%edge(:, pos) = real(plan%panel_geometry(idx)%edge_length, sp)
          plan%near_panels_sp%ox(:, pos) = real(plan%panel_geometry(idx)%edge_outward(1, :), sp)
          plan%near_panels_sp%oy(:, pos) = real(plan%panel_geometry(idx)%edge_outward(2, :), sp)
          plan%near_panels_sp%oz(:, pos) = real(plan%panel_geometry(idx)%edge_outward(3, :), sp)
          plan%near_panels_sp%area(pos) = real(plan%panel_geometry(idx)%area, sp)
        end do
      end do
    end subroutine build_near_panel_tables_single
  end subroutine build_near_panel_tables

  !> Source leaves that are not well separated from target_leaf. Well-separated subtrees reach the
  !! leaf through the M2L pairs, so the traversal only stops there.
  recursive subroutine gather_leaf_near_sources(plan, target_tree, target_leaf, source_node, near_buf, near_n)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_tree_type), intent(in) :: target_tree
    integer(i32), intent(in) :: target_leaf, source_node
    integer(i32), intent(inout) :: near_buf(:), near_n
    integer(i32) :: child_k, child_node

    if (nodes_well_separated(plan, target_tree, target_leaf, source_node)) return
    if (plan%source_tree%child_count(source_node) <= 0_i32) then
      near_n = near_n + 1_i32
      near_buf(near_n) = source_node
      return
    end if
    do child_k = 1_i32, plan%source_tree%child_count(source_node)
      child_node = plan%source_tree%child_idx(child_k, source_node)
      call gather_leaf_near_sources(plan, target_tree, target_leaf, child_node, near_buf, near_n)
    end do
  end subroutine gather_leaf_near_sources

  recursive subroutine accumulate_cached_pairs( &
    plan, target_tree, t_node, s_node, target_buf, source_buf, shift_idx1_buf, shift_idx2_buf, pair_used, pair_cap &
    )
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_tree_type), intent(in) :: target_tree
    integer(i32), intent(in) :: t_node, s_node
    integer(i32), allocatable, intent(inout) :: target_buf(:), source_buf(:), shift_idx1_buf(:), shift_idx2_buf(:)
    integer(i32), intent(inout) :: pair_used, pair_cap
    integer(i32) :: tc, sc, t_child, s_child

    tc = target_tree%child_count(t_node)
    sc = plan%source_tree%child_count(s_node)
    if (nodes_well_separated(plan, target_tree, t_node, s_node)) then
      call append_m2l_pair( &
        target_buf, source_buf, shift_idx1_buf, shift_idx2_buf, pair_used, pair_cap, t_node, s_node, 0_i32, 0_i32 &
        )
      return
    end if
    if (tc <= 0_i32 .and. sc <= 0_i32) return
    if (sc > 0_i32 .and. (tc <= 0_i32 .or. plan%source_tree%node_radius(s_node) >= target_tree%node_radius(t_node))) then
      do s_child = 1_i32, sc
        call accumulate_cached_pairs( &
          plan, target_tree, t_node, plan%source_tree%child_idx(s_child, s_node), target_buf, source_buf, shift_idx1_buf, &
          shift_idx2_buf, pair_used, pair_cap &
          )
      end do
    else
      do t_child = 1_i32, tc
        call accumulate_cached_pairs( &
          plan, target_tree, target_tree%child_idx(t_child, t_node), s_node, &
          target_buf, source_buf, shift_idx1_buf, shift_idx2_buf, pair_used, pair_cap &
          )
      end do
    end if
  end subroutine accumulate_cached_pairs

  subroutine append_m2l_pair( &
    target_buf, source_buf, shift_idx1_buf, shift_idx2_buf, pair_used, pair_cap, t_node, s_node, shift_idx1, shift_idx2 &
    )
    integer(i32), allocatable, intent(inout) :: target_buf(:), source_buf(:), shift_idx1_buf(:), shift_idx2_buf(:)
    integer(i32), intent(inout) :: pair_used, pair_cap
    integer(i32), intent(in) :: t_node, s_node, shift_idx1, shift_idx2
    integer(i32), allocatable :: tmp_target(:), tmp_source(:), tmp_shift_idx1(:), tmp_shift_idx2(:)
    integer(i32) :: new_capacity

    if (pair_used >= pair_cap) then
      new_capacity = max(pair_cap*2_i32, pair_cap + 64_i32)
      allocate (tmp_target(new_capacity), tmp_source(new_capacity), &
                tmp_shift_idx1(new_capacity), tmp_shift_idx2(new_capacity))
      tmp_target = 0_i32
      tmp_source = 0_i32
      tmp_shift_idx1 = 0_i32
      tmp_shift_idx2 = 0_i32
      if (pair_used > 0_i32) then
        tmp_target(1:pair_used) = target_buf(1:pair_used)
        tmp_source(1:pair_used) = source_buf(1:pair_used)
        tmp_shift_idx1(1:pair_used) = shift_idx1_buf(1:pair_used)
        tmp_shift_idx2(1:pair_used) = shift_idx2_buf(1:pair_used)
      end if
      call move_alloc(tmp_target, target_buf)
      call move_alloc(tmp_source, source_buf)
      call move_alloc(tmp_shift_idx1, shift_idx1_buf)
      call move_alloc(tmp_shift_idx2, shift_idx2_buf)
      pair_cap = new_capacity
    end if
    pair_used = pair_used + 1_i32
    target_buf(pair_used) = t_node
    source_buf(pair_used) = s_node
    shift_idx1_buf(pair_used) = shift_idx1
    shift_idx2_buf(pair_used) = shift_idx2
  end subroutine append_m2l_pair

  !> M2L の組を target node 順に並べ替え、target ごとの先頭位置 `m2l_target_start` を作る。
  !! 同じ target の組は列挙順のまま残す（安定な並べ替え）。M2L はこの順に足すので、和の順番は列挙順と同じで、
  !! 1 つの target の組と微分表が連続した領域に並ぶ。
  subroutine sort_pairs_by_target(plan, n_target_nodes)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32), intent(in) :: n_target_nodes
    integer(i32) :: pair_idx, target_node
    integer(i32), allocatable :: counts(:), cursor(:), sorted_pos(:), buffer(:)

    allocate (plan%m2l_target_start(n_target_nodes + 1_i32))
    plan%m2l_target_start = 1_i32
    if (plan%m2l_pair_count <= 0_i32) return

    allocate (counts(n_target_nodes), cursor(n_target_nodes), sorted_pos(plan%m2l_pair_count))
    counts = 0_i32
    do pair_idx = 1_i32, plan%m2l_pair_count
      target_node = plan%m2l_target_nodes(pair_idx)
      counts(target_node) = counts(target_node) + 1_i32
    end do
    do target_node = 1_i32, n_target_nodes
      plan%m2l_target_start(target_node + 1_i32) = plan%m2l_target_start(target_node) + counts(target_node)
    end do
    cursor = plan%m2l_target_start(1:n_target_nodes)
    do pair_idx = 1_i32, plan%m2l_pair_count
      target_node = plan%m2l_target_nodes(pair_idx)
      sorted_pos(pair_idx) = cursor(target_node)
      cursor(target_node) = cursor(target_node) + 1_i32
    end do
    allocate (buffer(plan%m2l_pair_count))
    buffer = plan%m2l_target_nodes(1:plan%m2l_pair_count)
    plan%m2l_target_nodes(sorted_pos) = buffer
    buffer = plan%m2l_source_nodes(1:plan%m2l_pair_count)
    plan%m2l_source_nodes(sorted_pos) = buffer
    buffer = plan%m2l_shift_idx1(1:plan%m2l_pair_count)
    plan%m2l_shift_idx1(sorted_pos) = buffer
    buffer = plan%m2l_shift_idx2(1:plan%m2l_pair_count)
    plan%m2l_shift_idx2(sorted_pos) = buffer
  end subroutine sort_pairs_by_target

  !> target node を組の多い順（同数は node 番号順）に並べ、`m2l_lane_width` 個ずつ group にまとめる。
  !! group の step 数は最初の lane（組が最も多い）の組の数で、組の少ない lane の残りの step は空にする。
  subroutine group_targets_for_lanes(plan, n_target_nodes)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32), intent(in) :: n_target_nodes
    integer(i32) :: target_node, pair_count, max_count, position, group, lane, k
    integer(i32), allocatable :: count_start(:), order(:)

    plan%m2l_group_count = (n_target_nodes + m2l_lane_width - 1_i32)/m2l_lane_width
    allocate (plan%m2l_group_target(m2l_lane_width, max(1_i32, plan%m2l_group_count)))
    allocate (plan%m2l_group_start(plan%m2l_group_count + 1_i32))
    plan%m2l_group_target = 0_i32
    plan%m2l_group_start = 1_i32

    ! 組の数ごとの数え上げで、多い順・同数は node 番号順の並びを作る。
    max_count = 0_i32
    do target_node = 1_i32, n_target_nodes
      max_count = max(max_count, plan%m2l_target_start(target_node + 1_i32) - plan%m2l_target_start(target_node))
    end do
    allocate (count_start(0:max_count + 1_i32), order(max(1_i32, n_target_nodes)))
    count_start = 0_i32
    do target_node = 1_i32, n_target_nodes
      pair_count = plan%m2l_target_start(target_node + 1_i32) - plan%m2l_target_start(target_node)
      count_start(pair_count) = count_start(pair_count) + 1_i32
    end do
    position = 1_i32
    do pair_count = max_count, 0_i32, -1_i32
      k = count_start(pair_count)
      count_start(pair_count) = position
      position = position + k
    end do
    do target_node = 1_i32, n_target_nodes
      pair_count = plan%m2l_target_start(target_node + 1_i32) - plan%m2l_target_start(target_node)
      order(count_start(pair_count)) = target_node
      count_start(pair_count) = count_start(pair_count) + 1_i32
    end do

    do group = 1_i32, plan%m2l_group_count
      do lane = 1_i32, m2l_lane_width
        position = (group - 1_i32)*m2l_lane_width + lane
        if (position <= n_target_nodes) plan%m2l_group_target(lane, group) = order(position)
      end do
      target_node = plan%m2l_group_target(1, group)
      plan%m2l_group_start(group + 1_i32) = plan%m2l_group_start(group) + &
                                            plan%m2l_target_start(target_node + 1_i32) - plan%m2l_target_start(target_node)
    end do

    allocate (plan%m2l_step_source(m2l_lane_width, max(1_i32, plan%m2l_group_start(plan%m2l_group_count + 1_i32) - 1_i32)))
    plan%m2l_step_source = 0_i32
    do group = 1_i32, plan%m2l_group_count
      do lane = 1_i32, m2l_lane_width
        target_node = plan%m2l_group_target(lane, group)
        if (target_node <= 0_i32) cycle
        do k = 0_i32, plan%m2l_target_start(target_node + 1_i32) - plan%m2l_target_start(target_node) - 1_i32
          plan%m2l_step_source(lane, plan%m2l_group_start(group) + k) = &
            plan%m2l_source_nodes(plan%m2l_target_start(target_node) + k)
        end do
      end do
    end do
  end subroutine group_targets_for_lanes

  subroutine precompute_translation_operators(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: alpha_idx, beta_idx, gamma_idx, term_idx
    integer(i32) :: node_idx, parent_node
    integer(i32) :: delta(3)
    integer(i32) :: max_m2m_terms, max_l2l_terms
    real(dp) :: d(3)

    if (allocated(plan%m2m_term_count)) deallocate (plan%m2m_term_count)
    if (allocated(plan%m2m_alpha_list)) deallocate (plan%m2m_alpha_list)
    if (allocated(plan%m2m_delta_list)) deallocate (plan%m2m_delta_list)
    if (allocated(plan%l2l_term_count)) deallocate (plan%l2l_term_count)
    if (allocated(plan%l2l_gamma_list)) deallocate (plan%l2l_gamma_list)
    if (allocated(plan%l2l_delta_list)) deallocate (plan%l2l_delta_list)

    allocate (plan%m2m_term_count(plan%ncoef), plan%l2l_term_count(plan%ncoef))
    plan%m2m_term_count = 0_i32
    plan%l2l_term_count = 0_i32
    do beta_idx = 1_i32, plan%ncoef
      do alpha_idx = 1_i32, plan%ncoef
        if (any(plan%alpha(:, alpha_idx) > plan%alpha(:, beta_idx))) cycle
        plan%m2m_term_count(beta_idx) = plan%m2m_term_count(beta_idx) + 1_i32
      end do
    end do
    do alpha_idx = 1_i32, plan%ncoef
      do gamma_idx = 1_i32, plan%ncoef
        if (any(plan%alpha(:, gamma_idx) < plan%alpha(:, alpha_idx))) cycle
        plan%l2l_term_count(alpha_idx) = plan%l2l_term_count(alpha_idx) + 1_i32
      end do
    end do

    max_m2m_terms = max(1_i32, maxval(plan%m2m_term_count))
    max_l2l_terms = max(1_i32, maxval(plan%l2l_term_count))
    allocate (plan%m2m_alpha_list(max_m2m_terms, plan%ncoef), plan%m2m_delta_list(max_m2m_terms, plan%ncoef))
    allocate (plan%l2l_gamma_list(max_l2l_terms, plan%ncoef), plan%l2l_delta_list(max_l2l_terms, plan%ncoef))
    plan%m2m_alpha_list = 0_i32
    plan%m2m_delta_list = 0_i32
    plan%l2l_gamma_list = 0_i32
    plan%l2l_delta_list = 0_i32

    do beta_idx = 1_i32, plan%ncoef
      term_idx = 0_i32
      do alpha_idx = 1_i32, plan%ncoef
        if (any(plan%alpha(:, alpha_idx) > plan%alpha(:, beta_idx))) cycle
        delta = plan%alpha(:, beta_idx) - plan%alpha(:, alpha_idx)
        term_idx = term_idx + 1_i32
        plan%m2m_alpha_list(term_idx, beta_idx) = alpha_idx
        plan%m2m_delta_list(term_idx, beta_idx) = plan%alpha_map(delta(1), delta(2), delta(3))
      end do
    end do
    do alpha_idx = 1_i32, plan%ncoef
      term_idx = 0_i32
      do gamma_idx = 1_i32, plan%ncoef
        if (any(plan%alpha(:, gamma_idx) < plan%alpha(:, alpha_idx))) cycle
        delta = plan%alpha(:, gamma_idx) - plan%alpha(:, alpha_idx)
        term_idx = term_idx + 1_i32
        plan%l2l_gamma_list(term_idx, alpha_idx) = gamma_idx
        plan%l2l_delta_list(term_idx, alpha_idx) = plan%alpha_map(delta(1), delta(2), delta(3))
      end do
    end do

    if (allocated(plan%source_shift_monomial)) deallocate (plan%source_shift_monomial)
    allocate (plan%source_shift_monomial(plan%ncoef, max(1_i32, plan%source_tree%nnode)))
    plan%source_shift_monomial = 0.0d0
    if (plan%source_tree%nnode > 0_i32) plan%source_shift_monomial(1_i32, 1_i32) = 1.0d0
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(plan) private(node_idx, parent_node, d)
    do node_idx = 1_i32, plan%source_tree%nnode
      parent_node = plan%source_tree%parent(node_idx)
      if (parent_node <= 0_i32) cycle
      d = plan%source_tree%node_center(:, node_idx) - plan%source_tree%node_center(:, parent_node)
      call fill_shift_monomial(plan, d, plan%source_shift_monomial(:, node_idx))
    end do
    !$omp end parallel do

    if (plan%target_tree_ready) then
      call precompute_target_shift_monomials(plan, plan%target_tree)
    else
      call precompute_target_shift_monomials(plan, plan%source_tree)
    end if
  end subroutine precompute_translation_operators

  !> L2L monomials of the evaluation tree (plan%target_tree or plan%source_tree), which is only read.
  subroutine precompute_target_shift_monomials(plan, target_tree)
    type(fmm_plan_type), intent(inout) :: plan
    type(fmm_tree_type), intent(in) :: target_tree
    integer(i32) :: node_idx, n_target_nodes, parent_node
    real(dp) :: d(3)

    n_target_nodes = target_tree%nnode
    if (allocated(plan%target_shift_monomial)) deallocate (plan%target_shift_monomial)
    allocate (plan%target_shift_monomial(plan%ncoef, max(1_i32, n_target_nodes)))
    plan%target_shift_monomial = 0.0d0
    if (n_target_nodes > 0_i32) plan%target_shift_monomial(1_i32, 1_i32) = 1.0d0
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(plan, target_tree, n_target_nodes) private(node_idx, parent_node, d)
    do node_idx = 1_i32, n_target_nodes
      parent_node = target_tree%parent(node_idx)
      if (parent_node <= 0_i32) cycle
      d = target_tree%node_center(:, node_idx) - target_tree%node_center(:, parent_node)
      call fill_shift_monomial(plan, d, plan%target_shift_monomial(:, node_idx))
    end do
    !$omp end parallel do
  end subroutine precompute_target_shift_monomials

  !> Monomials `d^alpha/alpha!` of one parent-to-child offset, shared by the source (M2M) and target (L2L) trees.
  subroutine fill_shift_monomial(plan, d, monomial)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: d(3)
    real(dp), intent(out) :: monomial(:)
    integer(i32) :: delta_idx
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    call build_axis_powers(d, plan%options%order, xpow, ypow, zpow)
    do delta_idx = 1_i32, plan%ncoef
      monomial(delta_idx) = xpow(plan%alpha(1, delta_idx))*ypow(plan%alpha(2, delta_idx)) &
                            *zpow(plan%alpha(3, delta_idx))/plan%alpha_factorial(delta_idx)
    end do
  end subroutine fill_shift_monomial

  subroutine precompute_m2l_derivatives(plan)
    type(fmm_plan_type), intent(inout) :: plan

    if (plan%target_tree_ready) then
      call precompute_m2l_derivatives_for(plan, plan%target_tree)
    else
      call precompute_m2l_derivatives_for(plan, plan%source_tree)
    end if
  end subroutine precompute_m2l_derivatives

  !> M2L derivatives of every cached pair, stored at the pair's lane and step (see `m2l_group_target`);
  !! target nodes belong to target_tree (plan%target_tree or plan%source_tree).
  subroutine precompute_m2l_derivatives_for(plan, target_tree)
    type(fmm_plan_type), intent(inout) :: plan
    type(fmm_tree_type), intent(in) :: target_tree
    integer(i32) :: group, lane, k, pair_idx, img_i, img_j, nshift, step_count
    integer(i32) :: axis1, axis2, s_node, t_node, shift_idx1, shift_idx2
    real(dp) :: base(3), rimg(3)
    real(dp) :: deriv(plan%nderiv), deriv_tmp(plan%nderiv)

    if (plan%m2l_pair_count <= 0_i32) then
      allocate (plan%m2l_deriv(m2l_lane_width, plan%nderiv, 1))
      plan%m2l_deriv = 0.0d0
      return
    end if
    step_count = plan%m2l_group_start(plan%m2l_group_count + 1_i32) - 1_i32
    allocate (plan%m2l_deriv(m2l_lane_width, plan%nderiv, step_count))
    plan%m2l_deriv = 0.0d0

    axis1 = 0_i32
    axis2 = 0_i32
    if (plan%options%use_periodic2) then
      axis1 = plan%options%periodic_axes(1)
      axis2 = plan%options%periodic_axes(2)
    end if
    nshift = size(plan%shift_axis1)
    !$omp parallel do default(none) schedule(dynamic, 1) &
    !$omp shared(plan, target_tree, axis1, axis2, nshift) &
    !$omp private(group, lane, k, pair_idx, img_i, img_j, s_node, t_node, base, shift_idx1, shift_idx2) &
    !$omp private(rimg, deriv, deriv_tmp)
    do group = 1_i32, plan%m2l_group_count
      do lane = 1_i32, m2l_lane_width
        t_node = plan%m2l_group_target(lane, group)
        if (t_node <= 0_i32) cycle
        do k = 0_i32, plan%m2l_target_start(t_node + 1_i32) - plan%m2l_target_start(t_node) - 1_i32
          pair_idx = plan%m2l_target_start(t_node) + k
          s_node = plan%m2l_source_nodes(pair_idx)
          base = target_tree%node_center(:, t_node) - plan%source_tree%node_center(:, s_node)
          deriv = 0.0d0
          if (.not. plan%options%use_periodic2) then
            call compute_laplace_derivatives(plan, base, deriv_tmp)
            deriv = deriv + deriv_tmp
          else
            shift_idx1 = plan%m2l_shift_idx1(pair_idx)
            shift_idx2 = plan%m2l_shift_idx2(pair_idx)
            if (shift_idx1 > 0_i32 .and. shift_idx2 > 0_i32) then
              rimg = periodic_image_offset(base, axis1, axis2, plan%shift_axis1(shift_idx1), plan%shift_axis2(shift_idx2))
              call compute_laplace_derivatives(plan, rimg, deriv_tmp)
              deriv = deriv + deriv_tmp
            else
              do img_i = 1_i32, nshift
                do img_j = 1_i32, nshift
                  rimg = periodic_image_offset(base, axis1, axis2, plan%shift_axis1(img_i), plan%shift_axis2(img_j))
                  call compute_laplace_derivatives(plan, rimg, deriv_tmp)
                  deriv = deriv + deriv_tmp
                end do
              end do
            end if
          end if
          plan%m2l_deriv(lane, :, plan%m2l_group_start(group) + k) = deriv
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine precompute_m2l_derivatives_for

  !> Offset from the source image shifted by (shift1, shift2) along the two periodic axes to the target.
  pure function periodic_image_offset(base, axis1, axis2, shift1, shift2) result(rimg)
    real(dp), intent(in) :: base(3), shift1, shift2
    integer(i32), intent(in) :: axis1, axis2
    real(dp) :: rimg(3)

    rimg = base
    rimg(axis1) = rimg(axis1) - shift1
    rimg(axis2) = rimg(axis2) - shift2
  end function periodic_image_offset

  logical function nodes_well_separated_shifted(plan, target_tree, target_node, source_node, shift1, shift2)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_tree_type), intent(in) :: target_tree
    integer(i32), intent(in) :: target_node, source_node
    real(dp), intent(in) :: shift1, shift2
    real(dp) :: d(3), rs, rt, lhs, rhs, theta_eff, target_center(3), shifted_center(3)
    integer(i32) :: axis1, axis2

    target_center = target_tree%node_center(:, target_node)
    rt = target_tree%node_radius(target_node)

    shifted_center = plan%source_tree%node_center(:, source_node)
    axis1 = 0_i32
    axis2 = 0_i32
    if (plan%options%use_periodic2) then
      axis1 = plan%options%periodic_axes(1)
      axis2 = plan%options%periodic_axes(2)
    end if
    if (axis1 > 0_i32) shifted_center(axis1) = shifted_center(axis1) + shift1
    if (axis2 > 0_i32) shifted_center(axis2) = shifted_center(axis2) + shift2

    d = target_center - shifted_center
    if (sum(d*d) <= 0.0d0) then
      nodes_well_separated_shifted = .false.
      return
    end if

    rs = plan%source_tree%node_radius(source_node)
    theta_eff = plan%options%theta
    lhs = (rs + rt)*(rs + rt)
    rhs = (theta_eff*theta_eff)*sum(d*d)
    nodes_well_separated_shifted = (lhs < rhs)
  end function nodes_well_separated_shifted

end module bem_coulomb_fmm_plan_ops
