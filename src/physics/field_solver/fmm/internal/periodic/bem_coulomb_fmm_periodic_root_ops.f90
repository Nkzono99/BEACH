!> periodic2 root operator の前計算。
module bem_coulomb_fmm_periodic_root_ops
!$ use omp_lib, only: omp_get_num_threads
  use bem_kinds, only: dp, i32, i64
  use bem_mpi, only: mpi_context, mpi_is_root, mpi_get_rank_size, &
                     mpi_allreduce_sum_real_dp_array, mpi_bcast_i32_array, mpi_bcast_real_dp_array
  use bem_coulomb_fmm_types, only: fmm_plan_type, periodic_operator_proxy_shell_scale, periodic_source_proxy_half_size
  use bem_coulomb_fmm_basis, only: build_axis_powers, initialize_basis_tables
  use bem_coulomb_fmm_periodic, only: use_periodic2_cached_kneq0
  use bem_coulomb_fmm_periodic_ewald, only: add_periodic2_exact_ewald_field_potential_single_source
  use bem_periodic_shell_mean, only: periodic_shell_mean_point_excess
  use bem_coulomb_fmm_periodic_cache, only: periodic_cache_fingerprint
  use bem_coulomb_fmm_periodic_cell_cache, only: periodic_cell_cache_lock_type, lookup_cell_cache, &
                                                 publish_cell_cache_and_unlock
  use bem_regularized_qr, only: regularized_qr_type, prepare_regularized_qr, solve_regularized_qr
  implicit none
  private

  real(dp), parameter :: periodic_operator_proxy_multiplier = 8.0d0
  integer(i32), parameter :: periodic_operator_proxy_null_degrees = 2_i32
  real(dp), parameter :: periodic_operator_inner_proxy_scale = 0.85_dp
  real(dp), parameter :: periodic_operator_check_multiplier = 24.0d0
  real(dp), parameter :: periodic_operator_check_shell_scale = 0.92d0
  real(dp), parameter :: periodic_operator_shell_spacing = 1.0_dp/64.0_dp
  real(dp), parameter :: periodic_operator_lstsq_ridge = 1.0d-12
  real(dp), parameter :: periodic_operator_qr_tol = 1.0d-12
  real(dp), parameter :: periodic_operator_seam_weight_multiplier = 1.0d0

  type :: periodic_check_fit_type
    integer(i32), allocatable :: index(:)
    real(dp), allocatable :: weights(:), potential_weights(:), potential_mean_basis(:)
    type(regularized_qr_type) :: factorization
  end type periodic_check_fit_type

  type :: proxy_shell_table_type
    real(dp) :: h0 = 0.0_dp, dh = 0.0_dp
    real(dp), allocatable :: value(:, :, :), slope(:, :, :)
  end type proxy_shell_table_type

  public :: precompute_periodic_root_operator

contains

  !> periodic2 の root operator を前計算する。
  !! @param[inout] plan FMM 計画。
  subroutine precompute_periodic_root_operator(plan, reuse_translations)
    type(fmm_plan_type), intent(inout) :: plan
    ! Internal verification path: false computes every pair independently and bypasses disk cache.
    logical, intent(in), optional :: reuse_translations
    logical :: reuse

    if (allocated(plan%periodic_root_operator)) deallocate (plan%periodic_root_operator)
    if (allocated(plan%periodic_operator_cell_slot)) deallocate (plan%periodic_operator_cell_slot)
    if (allocated(plan%periodic_operator_cells)) deallocate (plan%periodic_operator_cells)
    plan%periodic_operator_built_count = 0_i32
    plan%periodic_operator_new_cell_count = 0_i32
    plan%periodic_root_operator_ready = .false.
    plan%periodic_root_target_count = 0_i32

    plan%periodic_cache_hit = .false.
    plan%periodic_cache_fingerprint = ''
    plan%periodic_cache_path = ''
    if (use_periodic2_cached_kneq0(plan)) then
      reuse = .true.
      if (present(reuse_translations)) reuse = reuse_translations
      ! MPI の初期化・終了は呼び出し側の責務。rank の分担は plan の options が持つ文脈に従う。
      call precompute_periodic_cached_operator_collective(plan, plan%options%mpi, reuse)
    end if
  end subroutine precompute_periodic_root_operator

  !> cached_kneq0 の root operator を組み立てる。cache にあるセルは MPI root が読み、無いセルは全 rank で分担して
  !! 作る。cache の読み書きとロックは root だけが bem_coulomb_fmm_periodic_cell_cache で行い、collective の順序は
  !! 全 rank で同じ。
  subroutine precompute_periodic_cached_operator_collective(plan, mpi, reuse_translations)
    type(fmm_plan_type), intent(inout) :: plan
    type(mpi_context), intent(in) :: mpi
    logical, intent(in) :: reuse_translations
    integer(i32) :: k, target_count, ncell, nreq, nmiss, m, slot, cell
    real(dp), allocatable :: flat_operator(:)
    type(periodic_cell_cache_lock_type) :: cache_lock
    logical, allocatable :: required(:), have(:)
    integer(i32), allocatable :: missing_flag(:), miss_cells(:), req_cells(:), cached_cells(:), cached_slot(:), &
                                 miss_slot(:)
    real(dp), allocatable :: cached_blocks(:, :, :), new_blocks(:, :, :)

    if (.not. plan%periodic_ewald%ready) return
    if (plan%ncoef <= 1_i32) return
    if (plan%periodic_operator_anchor_count <= 0_i32) return

    target_count = plan%periodic_operator_anchor_count
    ncell = plan%periodic_root_source_count
    ! The operator is mesh-independent per cell, but only cells that hold sources are needed now.
    allocate (required(ncell), have(ncell), missing_flag(ncell))
    required = .false.
    do k = 1_i32, plan%nsrc
      required(plan%periodic_source_cell(k)) = .true.
    end do
    nreq = count(required)
    req_cells = pack([(k, k=1_i32, ncell)], required)

    call set_cache_identity(plan)
    have = .false.
    if (reuse_translations .and. mpi_is_root(mpi)) then
      call lookup_cell_cache(plan, target_count, required, cached_cells, cached_blocks, have, cache_lock)
    end if
    missing_flag = merge(1_i32, 0_i32, required .and. .not. have)
    if (reuse_translations) call mpi_bcast_i32_array(mpi, missing_flag, 0_i32)
    nmiss = count(missing_flag == 1_i32)
    miss_cells = pack([(k, k=1_i32, ncell)], missing_flag == 1_i32)
    plan%periodic_operator_new_cell_count = nmiss
    if (nmiss > 0_i32) call build_missing_cell_operators(plan, mpi, reuse_translations, miss_cells, new_blocks)

    ! Compact operator for the cells this mesh uses: pair (slot-1)*anchors + anchor.
    allocate (plan%periodic_operator_cell_slot(ncell))
    plan%periodic_operator_cell_slot = 0_i32
    plan%periodic_operator_cell_slot(req_cells) = [(slot, slot=1_i32, nreq)]
    plan%periodic_operator_built_count = nreq
    plan%periodic_operator_cells = req_cells
    plan%periodic_root_target_count = target_count*nreq
    allocate (plan%periodic_root_operator(plan%ncoef, plan%ncoef, target_count*nreq))
    plan%periodic_root_operator = 0.0_dp
    allocate (miss_slot(ncell), cached_slot(ncell))
    miss_slot = 0_i32
    if (nmiss > 0_i32) miss_slot(miss_cells) = [(m, m=1_i32, nmiss)]
    cached_slot = 0_i32
    if (allocated(cached_cells)) then
      if (size(cached_cells) > 0) cached_slot(cached_cells) = [(m, m=1_i32, int(size(cached_cells), i32))]
    end if
    if (mpi_is_root(mpi) .or. .not. reuse_translations) then
      do slot = 1_i32, nreq
        cell = req_cells(slot)
        if (miss_slot(cell) > 0_i32) then
          m = miss_slot(cell)
          plan%periodic_root_operator(:, :, (slot - 1_i32)*target_count + 1_i32:slot*target_count) = &
            new_blocks(:, :, (m - 1_i32)*target_count + 1_i32:m*target_count)
        else
          m = cached_slot(cell)
          if (m <= 0_i32) error stop 'periodic operator cell is neither cached nor built.'
          plan%periodic_root_operator(:, :, (slot - 1_i32)*target_count + 1_i32:slot*target_count) = &
            cached_blocks(:, :, (m - 1_i32)*target_count + 1_i32:m*target_count)
        end if
      end do
    end if
    if (reuse_translations) then
      if (mpi_is_root(mpi) .and. cache_lock%held) then
        call publish_cell_cache_and_unlock(plan, target_count, cached_cells, cached_blocks, miss_cells, new_blocks, &
                                           cache_lock)
      end if
      ! Cached cells live only on the root; share the assembled compact operator.
      allocate (flat_operator(size(plan%periodic_root_operator)))
      if (mpi_is_root(mpi)) flat_operator = reshape(plan%periodic_root_operator, [size(flat_operator)])
      call mpi_bcast_real_dp_array(mpi, flat_operator, 0_i32)
      plan%periodic_root_operator = reshape(flat_operator, shape(plan%periodic_root_operator))
      deallocate (flat_operator)
    end if
    plan%periodic_cache_hit = reuse_translations .and. nmiss == 0_i32
    plan%periodic_root_operator_ready = .true.
  end subroutine precompute_periodic_cached_operator_collective

  !> cache に無いセル miss_cells の root operator ブロックを作る。proxy と check 点の fit を pair の translation class
  !! ごとに解き、class を rank で分担して最後に全 rank で足し合わせる。ファイル I/O はしない。
  !! new_blocks の pair (m-1)*anchors + anchor が miss_cells(m) のセルに対応する。
  subroutine build_missing_cell_operators(plan, mpi, reuse_translations, miss_cells, new_blocks)
    type(fmm_plan_type), intent(inout) :: plan
    type(mpi_context), intent(in) :: mpi
    logical, intent(in) :: reuse_translations
    integer(i32), intent(in) :: miss_cells(:)
    real(dp), allocatable, intent(out) :: new_blocks(:, :, :)
    integer(i32) :: nproxy, ncheck_bulk, ncheck_face, ncheck_max, ncheck, j, i, k, face, kind
    integer(i32) :: proxy_order, even_degree, even_count
    integer(i32) :: target_idx, target_count, rank, rank_count, source_idx, pair_count
    integer(i32) :: pair, class_idx, class_count, key(3), axis, mask, check_idx, proxy_idx, reflected_check
    integer(i32) :: reflection, shell_column, cells(3), cell_index(3), nmiss, m
    real(dp) :: source_half(3), proxy_half(3), target_center(3), target_half(3), delta(3), unit(3), point(3), target_z
    real(dp), allocatable :: proxy_points(:, :), check_points(:, :), selected_points(:, :)
    real(dp), allocatable :: proxy_to_multipole(:, :), proxy_to_local(:, :), proxy_pinv(:, :), full_proxy_pinv(:, :)
    real(dp), allocatable :: field_matrix(:, :), field_rhs(:), coeff(:), potential_rhs(:), flat_operator(:)
    real(dp), allocatable :: residual_e(:, :, :), residual_phi(:, :), class_delta(:, :)
    real(dp) :: e_res(3), phi_res, shell_excess, shell_excess_slope, signs(3)
    integer(i32), allocatable :: pair_class(:), class_key(:, :), target_kind(:), pair_reflection(:)
    type(fmm_plan_type) :: proxy_basis
    type(proxy_shell_table_type) :: shell_table
    type(periodic_check_fit_type) :: fit(0:15)
    logical :: used_kind(0:15)

    call mpi_get_rank_size(rank, rank_count, mpi)
    proxy_order = plan%options%order + periodic_operator_proxy_null_degrees
    even_degree = proxy_order/2_i32
    even_count = (even_degree + 1_i32)*(even_degree + 2_i32)*(even_degree + 3_i32)/6_i32
    nproxy = max(int(periodic_operator_proxy_multiplier*real(plan%ncoef, dp), i32), &
                 8_i32*((3_i32*even_count + 1_i32)/2_i32))
    ncheck_bulk = max(8_i32*plan%ncoef, int(periodic_operator_check_multiplier*real(plan%ncoef, dp), i32))
    ncheck_face = 4_i32*((max(8_i32, plan%ncoef) + 3_i32)/4_i32)
    ncheck_max = ncheck_bulk + 4_i32*ncheck_face
    target_count = plan%periodic_operator_anchor_count
    cells = plan%periodic_operator_source_cells
    nmiss = int(size(miss_cells), i32)
    pair_count = target_count*nmiss

    source_half = plan%periodic_source_half
    target_half = plan%periodic_anchor_half
    allocate (proxy_points(3, nproxy), check_points(3, ncheck_max))
    proxy_basis%options%order = proxy_order
    call initialize_basis_tables(proxy_basis, proxy_order)
    allocate (proxy_to_multipole(proxy_basis%ncoef, nproxy), proxy_to_local(plan%ncoef, nproxy))
    allocate (full_proxy_pinv(nproxy, proxy_basis%ncoef))
    allocate (proxy_pinv(nproxy, plan%ncoef), field_matrix(3_i32*ncheck_max, plan%ncoef - 1_i32))
    proxy_half = periodic_source_proxy_half_size(source_half)
    call build_root_surface_points([0.0_dp, 0.0_dp, 0.0_dp], proxy_half, nproxy, 0.13_dp, 1.0_dp, proxy_points)
    ! Two shells avoid the polynomial dependency of a single cuboid surface at degree six.
    do j = 1_i32, nproxy
      if (mod((j - 1_i32)/8_i32, 2_i32) == 0_i32) then
        proxy_points(:, j) = periodic_operator_inner_proxy_scale*proxy_points(:, j)
      end if
    end do
    call build_proxy_shell_table(plan, proxy_points, shell_table)
    call build_proxy_multipole_matrix(proxy_basis, [0.0_dp, 0.0_dp, 0.0_dp], proxy_points, proxy_to_multipole)
    call build_minimum_norm_pseudoinverse(proxy_to_multipole, full_proxy_pinv)
    ! Additional proxy moments are zero: the proxy must not invent the first omitted source moments.
    proxy_pinv = full_proxy_pinv(:, 1:plan%ncoef)
    call build_root_surface_points([0.0_dp, 0.0_dp, 0.0_dp], target_half, ncheck_bulk, 0.37_dp, &
                                   periodic_operator_check_shell_scale, check_points(:, 1:ncheck_bulk))
    ! Use the union of all four seam patterns for the expensive kernel. Each fit selects only its own faces.
    ncheck = ncheck_bulk
    do face = 0_i32, 3_i32
      axis = plan%options%periodic_axes(face/2_i32 + 1_i32)
      call append_face_check_points([0.0_dp, 0.0_dp, 0.0_dp], target_half, axis, &
                                    real(2_i32*mod(face, 2_i32) - 1_i32, dp)*target_half(axis), &
                                    ncheck_face, check_points, ncheck)
    end do
    allocate (target_kind(target_count), pair_reflection(pair_count), pair_class(pair_count), &
              class_key(3, pair_count), class_delta(3, pair_count))
    used_kind = .false.
    class_count = 0_i32
    ! Anchor spacings are power-of-two multiples of the source spacing, so offsets are multiples of unit.
    unit = plan%periodic_source_half
    do target_idx = 1_i32, target_count
      target_center = plan%periodic_anchor_center(:, target_idx)
      mask = 0_i32
      do face = 0_i32, 3_i32
        axis = plan%options%periodic_axes(face/2_i32 + 1_i32)
        if (mod(face, 2_i32) == 0_i32) then
          point(axis) = plan%options%target_box_min(axis)
        else
          point(axis) = plan%options%target_box_max(axis)
        end if
        if (abs(target_center(axis) + real(2_i32*mod(face, 2_i32) - 1_i32, dp)*target_half(axis) - &
                point(axis)) < 64.0_dp*epsilon(1.0_dp)*maxval(plan%options%target_box_max - &
                                                              plan%options%target_box_min)) then
          mask = ibset(mask, face)
        end if
      end do
      target_kind(target_idx) = mask
      used_kind(mask) = .true.
      do m = 1_i32, nmiss
        source_idx = miss_cells(m)
        pair = (m - 1_i32)*target_count + target_idx
        delta = target_center - plan%periodic_source_center(:, source_idx)
        key = nint(delta/unit, kind=i32)
        pair_reflection(pair) = 0_i32
        if (reuse_translations) then
          do axis = 1_i32, 3_i32
            if (key(axis) < 0_i32) pair_reflection(pair) = ibset(pair_reflection(pair), axis - 1_i32)
          end do
          key = abs(key)
          delta = abs(delta)
        end if
        class_idx = class_count + 1_i32
        if (reuse_translations) then
          do k = 1_i32, class_count
            if (.not. all(class_key(:, k) == key)) cycle
            class_idx = k
            exit
          end do
        end if
        if (class_idx > class_count) then
          class_count = class_idx
          class_key(:, class_idx) = key
          class_delta(:, class_idx) = delta
        end if
        pair_class(pair) = class_idx
      end do
    end do
    plan%periodic_operator_class_count = class_count
    ! QR is translation invariant; prepare it once per seam mask on each MPI rank.
    do kind = 0_i32, 15_i32
      if (.not. used_kind(kind)) cycle
      ncheck = ncheck_bulk + popcnt(kind)*ncheck_face
      allocate (fit(kind)%index(ncheck), fit(kind)%weights(ncheck), selected_points(3, ncheck))
      fit(kind)%index(1:ncheck_bulk) = [(i, i=1_i32, ncheck_bulk)]
      k = ncheck_bulk
      do face = 0_i32, 3_i32
        if (.not. btest(kind, face)) cycle
        fit(kind)%index(k + 1_i32:k + ncheck_face) = &
          [(ncheck_bulk + face*ncheck_face + i, i=1_i32, ncheck_face)]
        k = k + ncheck_face
      end do
      selected_points = check_points(:, fit(kind)%index)
      call set_periodic_check_weights(ncheck_bulk, ncheck, fit(kind)%weights)
      allocate (fit(kind)%potential_weights(ncheck), fit(kind)%potential_mean_basis(plan%ncoef))
      fit(kind)%potential_weights = fit(kind)%weights**2/sum(fit(kind)%weights**2)
      call prepare_potential_mean_basis(plan, selected_points, fit(kind)%potential_weights, fit(kind)%potential_mean_basis)
      call build_local_field_matrix(plan, [0.0_dp, 0.0_dp, 0.0_dp], selected_points, field_matrix(1:3_i32*ncheck, :))
      call apply_field_check_weights(field_matrix(1:3_i32*ncheck, :), fit(kind)%weights)
      call prepare_regularized_qr(fit(kind)%factorization, field_matrix(1:3_i32*ncheck, :), &
                                  periodic_operator_lstsq_ridge, periodic_operator_qr_tol)
      plan%periodic_qr_preparation_count = plan%periodic_qr_preparation_count + 1_i32
      deallocate (selected_points)
    end do
    allocate (new_blocks(plan%ncoef, plan%ncoef, pair_count))
    allocate (residual_e(3, ncheck_max, nproxy), residual_phi(ncheck_max, nproxy))
    new_blocks = 0.0_dp
    plan%periodic_operator_thread_count = max(plan%periodic_operator_thread_count, 1_i32)
    do class_idx = rank + 1_i32, class_count, rank_count
      delta = class_delta(:, class_idx)
!$omp parallel do default(shared) private(j, i, e_res, phi_res) schedule(static)
      do j = 1_i32, nproxy
        do i = 1_i32, ncheck_max
          e_res = 0.0_dp
          phi_res = 0.0_dp
          call add_periodic2_exact_ewald_field_potential_single_source( &
            plan, 1.0_dp, proxy_points(:, j), delta + check_points(:, i), e_res, phi_res)
          residual_e(:, i, j) = e_res
          residual_phi(i, j) = phi_res
        end do
      end do
!$omp end parallel do
      plan%periodic_operator_ewald_evaluation_count = plan%periodic_operator_ewald_evaluation_count + &
                                                      int(nproxy, i64)*int(ncheck_max, i64)
      do pair = 1_i32, pair_count
        if (pair_class(pair) /= class_idx) cycle
        target_idx = mod(pair - 1_i32, target_count) + 1_i32
        source_idx = miss_cells((pair - 1_i32)/target_count + 1_i32)
        cell_index = nint((plan%periodic_source_center(:, source_idx) - plan%options%target_box_min)/ &
                          (plan%options%target_box_max - plan%options%target_box_min)*real(cells, dp) - 0.5_dp, kind=i32)
        shell_column = 1_i32 + cell_index(plan%options%periodic_axes(1)) + &
                       cells(plan%options%periodic_axes(1))*cell_index(plan%options%periodic_axes(2))
        reflection = pair_reflection(pair)
        do axis = 1_i32, 3_i32
          signs(axis) = merge(-1.0_dp, 1.0_dp, btest(reflection, axis - 1_i32))
        end do
        kind = target_kind(target_idx)
        ncheck = size(fit(kind)%index)
        target_center = plan%periodic_anchor_center(:, target_idx)
        proxy_to_local = 0.0_dp
        plan%periodic_operator_local_target_count = plan%periodic_operator_local_target_count + 1_i32
!$omp parallel default(shared) private(j, i, check_idx, proxy_idx, reflected_check, field_rhs) &
!$omp private(potential_rhs, coeff, e_res, phi_res, shell_excess, shell_excess_slope, point, target_z)
        allocate (field_rhs(3_i32*ncheck), potential_rhs(ncheck), coeff(plan%ncoef - 1_i32))
!$omp single
!$      plan%periodic_operator_thread_count = max(plan%periodic_operator_thread_count, int(omp_get_num_threads(), i32))
!$omp end single
!$omp do schedule(static)
        do j = 1_i32, nproxy
          proxy_idx = 8_i32*((j - 1_i32)/8_i32) + ieor(mod(j - 1_i32, 8_i32), reflection) + 1_i32
          point = plan%periodic_source_center(:, source_idx) + proxy_points(:, j)
          do i = 1_i32, ncheck
            check_idx = fit(kind)%index(i)
            reflected_check = reflected_check_index(check_idx, ncheck_bulk, ncheck_face, reflection, &
                                                    plan%options%periodic_axes)
            e_res = signs*residual_e(:, reflected_check, proxy_idx)
            phi_res = residual_phi(reflected_check, proxy_idx)
            ! The shell mean uses absolute x/y and must be added for each physical pair.
            if (plan%periodic_shell_mean%ready) then
              target_z = target_center(plan%periodic_shell_mean%free_axis) + &
                         check_points(plan%periodic_shell_mean%free_axis, check_idx)
              call lookup_proxy_shell(shell_table, j, shell_column, target_z - point(plan%periodic_shell_mean%free_axis), &
                                      shell_excess, shell_excess_slope)
              phi_res = phi_res + shell_excess
              e_res(plan%periodic_shell_mean%free_axis) = e_res(plan%periodic_shell_mean%free_axis) - shell_excess_slope
            end if
            potential_rhs(i) = phi_res
            field_rhs(i) = fit(kind)%weights(i)*e_res(1)
            field_rhs(ncheck + i) = fit(kind)%weights(i)*e_res(2)
            field_rhs(2_i32*ncheck + i) = fit(kind)%weights(i)*e_res(3)
          end do
          call solve_regularized_qr(fit(kind)%factorization, field_rhs, coeff)
          proxy_to_local(2:plan%ncoef, j) = coeff
          proxy_to_local(1, j) = dot_product(fit(kind)%potential_weights, potential_rhs) - &
                                 dot_product(fit(kind)%potential_mean_basis(2:plan%ncoef), coeff)
        end do
!$omp end do
        deallocate (field_rhs, potential_rhs, coeff)
!$omp end parallel
        new_blocks(:, :, pair) = matmul(proxy_to_local, proxy_pinv)
      end do
    end do
    allocate (flat_operator(size(new_blocks)))
    flat_operator = reshape(new_blocks, [size(flat_operator)])
    call mpi_allreduce_sum_real_dp_array(mpi, flat_operator)
    new_blocks = reshape(flat_operator, shape(new_blocks))
    deallocate (flat_operator)
    if (mpi_is_root(mpi)) plan%periodic_operator_build_count = plan%periodic_operator_build_count + 1_i32
  end subroutine build_missing_cell_operators

  !> Only the source x/y and the height difference enter s. Reuse a table across all source/target z cells.
  subroutine build_proxy_shell_table(plan, proxy_points, table)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: proxy_points(:, :)
    type(proxy_shell_table_type), intent(out) :: table
    integer(i32) :: cells(3), nz, column, j, k, axes(2), free_axis
    real(dp) :: span(3), point(3), extent, h

    if (.not. plan%periodic_shell_mean%ready) return
    axes = plan%options%periodic_axes
    free_axis = plan%periodic_shell_mean%free_axis
    span = plan%options%target_box_max - plan%options%target_box_min
    cells = plan%periodic_operator_source_cells
    extent = span(free_axis) + maxval(abs(proxy_points(free_axis, :)))
    nz = ceiling(2.0_dp*extent/(periodic_operator_shell_spacing*minval(plan%options%periodic_len)), kind=i32) + 1_i32
    table%h0 = -extent
    table%dh = 2.0_dp*extent/real(nz - 1_i32, dp)
    allocate (table%value(nz, size(proxy_points, 2), cells(axes(1))*cells(axes(2))), &
              table%slope(nz, size(proxy_points, 2), cells(axes(1))*cells(axes(2))))
!$omp parallel do collapse(2) default(shared) private(column, j, k, point, h) schedule(static)
    do column = 1_i32, cells(axes(1))*cells(axes(2))
      do j = 1_i32, size(proxy_points, 2)
        point = 0.0_dp
        point(axes) = plan%options%target_box_min(axes) + span(axes)/real(cells(axes), dp)* &
                      (real([mod(column - 1_i32, cells(axes(1))), (column - 1_i32)/cells(axes(1))], dp) + 0.5_dp) + &
                      proxy_points(axes, j)
        do k = 1_i32, nz
          h = table%h0 + real(k - 1_i32, dp)*table%dh
          call periodic_shell_mean_point_excess(plan%periodic_shell_mean, point, h, &
                                                table%value(k, j, column), table%slope(k, j, column))
        end do
      end do
    end do
!$omp end parallel do
  end subroutine build_proxy_shell_table

  pure subroutine lookup_proxy_shell(table, proxy, column, h, value, slope)
    type(proxy_shell_table_type), intent(in) :: table
    integer(i32), intent(in) :: proxy, column
    real(dp), intent(in) :: h
    real(dp), intent(out) :: value, slope
    integer(i32) :: k
    real(dp) :: t, t2, t3, v0, v1, s0, s1

    t = (h - table%h0)/table%dh
    k = min(max(floor(t, kind=i32), 0_i32), size(table%value, 1) - 2_i32)
    t = t - real(k, dp)
    k = k + 1_i32
    t2 = t*t
    t3 = t2*t
    v0 = table%value(k, proxy, column)
    v1 = table%value(k + 1_i32, proxy, column)
    s0 = table%slope(k, proxy, column)
    s1 = table%slope(k + 1_i32, proxy, column)
    value = (2.0_dp*t3 - 3.0_dp*t2 + 1.0_dp)*v0 + (t3 - 2.0_dp*t2 + t)*table%dh*s0 + &
            (-2.0_dp*t3 + 3.0_dp*t2)*v1 + (t3 - t2)*table%dh*s1
    slope = ((6.0_dp*t2 - 6.0_dp*t)*v0 + (-6.0_dp*t2 + 6.0_dp*t)*v1)/table%dh + &
            (3.0_dp*t2 - 4.0_dp*t + 1.0_dp)*s0 + (3.0_dp*t2 - 2.0_dp*t)*s1
  end subroutine lookup_proxy_shell

  subroutine set_cache_identity(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: integer_values(13)
    real(dp) :: real_values(16)
    character(len=256) :: tag

    integer_values = [plan%options%order, plan%ncoef, plan%options%periodic_axes, &
                      plan%options%periodic_image_layers, plan%options%periodic_ewald_layers, &
                      plan%periodic_operator_source_cells, plan%periodic_operator_anchor_cells, &
                      periodic_operator_proxy_null_degrees]
    real_values = [plan%options%target_box_min, plan%options%target_box_max, plan%options%periodic_len, &
                   plan%options%softening, plan%periodic_ewald%alpha, plan%options%periodic_generation_tolerance, &
                   periodic_operator_proxy_shell_scale, periodic_operator_seam_weight_multiplier, &
                   periodic_operator_check_shell_scale, periodic_operator_shell_spacing, periodic_operator_inner_proxy_scale]
    tag = 'periodic-kneq0-generator-v12|box-lattice|sparse-cell-cache-v1|padded-box-proxy-v1|null-moments-v1|'// &
          'scaled-proxy-qr-v1|shell-mean-v1|real64|'
    if (plan%panel_source) then
      tag = trim(tag)//'kernel=triangle_p0|'
    else
      tag = trim(tag)//'kernel=point_softened|'
    end if
    tag = trim(tag)//'zero=state_subtracted|normalization=core'
    plan%periodic_cache_fingerprint = periodic_cache_fingerprint(tag, integer_values, real_values)
    plan%periodic_cache_path = trim(plan%options%periodic_cache_dir)//'/operator-'// &
                               trim(plan%periodic_cache_fingerprint)//'.bin'
  end subroutine set_cache_identity

  !> Complete reflection orbits make the sampled least-squares problem reflection invariant.
  subroutine build_root_surface_points(center, half_size, npoint, offset, scale, points)
    real(dp), intent(in) :: center(3), half_size(3), offset, scale
    integer(i32), intent(in) :: npoint
    real(dp), intent(out) :: points(3, npoint)
    integer(i32) :: idx, orbit, code, axis, normal
    real(dp) :: base(3), signs(3)
    real(dp), parameter :: g1 = 0.7548776662466927_dp, g2 = 0.5698402909980532_dp

    do idx = 1_i32, npoint
      orbit = (idx - 1_i32)/8_i32
      code = mod(idx - 1_i32, 8_i32)
      normal = mod(orbit, 3_i32) + 1_i32
      base = 0.0_dp
      base(normal) = 1.0_dp
      base(mod(normal, 3_i32) + 1_i32) = 0.05_dp + 0.9_dp*modulo(offset + real(orbit + 1_i32, dp)*g1, 1.0_dp)
      base(mod(normal + 1_i32, 3_i32) + 1_i32) = 0.05_dp + 0.9_dp*modulo(offset + real(orbit + 1_i32, dp)*g2, 1.0_dp)
      do axis = 1_i32, 3_i32
        signs(axis) = merge(-1.0_dp, 1.0_dp, btest(code, axis - 1_i32))
      end do
      points(:, idx) = center + scale*half_size*base*signs
    end do
  end subroutine build_root_surface_points

  subroutine append_face_check_points(center, half_size, normal_axis, face_position, npoint, points, point_count)
    real(dp), intent(in) :: center(3), half_size(3), face_position
    integer(i32), intent(in) :: normal_axis, npoint
    real(dp), intent(inout) :: points(:, :)
    integer(i32), intent(inout) :: point_count
    integer(i32) :: idx, axis, tangent_axis(2), tangent_count, orbit, code
    real(dp) :: u, v
    real(dp), parameter :: g1 = 0.7548776662466927_dp, g2 = 0.5698402909980532_dp

    tangent_count = 0_i32
    do axis = 1_i32, 3_i32
      if (axis == normal_axis) cycle
      tangent_count = tangent_count + 1_i32
      tangent_axis(tangent_count) = axis
    end do
    do idx = 1_i32, npoint
      orbit = (idx - 1_i32)/4_i32
      code = mod(idx - 1_i32, 4_i32)
      u = (0.025_dp + 0.95_dp*modulo(0.29_dp + real(orbit + 1_i32, dp)*g1, 1.0_dp))* &
          merge(-1.0_dp, 1.0_dp, btest(code, 0))
      v = (0.025_dp + 0.95_dp*modulo(0.41_dp + real(orbit + 1_i32, dp)*g2, 1.0_dp))* &
          merge(-1.0_dp, 1.0_dp, btest(code, 1))
      point_count = point_count + 1_i32
      points(:, point_count) = center
      points(normal_axis, point_count) = face_position
      points(tangent_axis(1), point_count) = center(tangent_axis(1)) + u*half_size(tangent_axis(1))
      points(tangent_axis(2), point_count) = center(tangent_axis(2)) + v*half_size(tangent_axis(2))
    end do
  end subroutine append_face_check_points

  pure integer(i32) function reflected_check_index(idx, bulk, per_face, reflection, axes) result(mapped)
    integer(i32), intent(in) :: idx, bulk, per_face, reflection, axes(2)
    integer(i32) :: face, local, code, bit, axis, normal

    if (idx <= bulk) then
      mapped = 8_i32*((idx - 1_i32)/8_i32) + ieor(mod(idx - 1_i32, 8_i32), reflection) + 1_i32
      return
    end if
    face = (idx - bulk - 1_i32)/per_face
    local = mod(idx - bulk - 1_i32, per_face)
    normal = axes(face/2_i32 + 1_i32)
    if (btest(reflection, normal - 1_i32)) face = ieor(face, 1_i32)
    code = mod(local, 4_i32)
    bit = 0_i32
    do axis = 1_i32, 3_i32
      if (axis == normal) cycle
      if (btest(reflection, axis - 1_i32)) code = ieor(code, 2_i32**bit)
      bit = bit + 1_i32
    end do
    mapped = bulk + face*per_face + 4_i32*(local/4_i32) + code + 1_i32
  end function reflected_check_index

  subroutine set_periodic_check_weights(bulk_count, point_count, weights)
    integer(i32), intent(in) :: bulk_count, point_count
    real(dp), intent(out) :: weights(:)
    real(dp) :: seam_weight

    weights(1:point_count) = 1.0_dp
    if (point_count <= bulk_count) return
    ! seam集合とbulk集合の総least-squares重みを揃え、境界だけで内部精度を支配しないようにする。
    seam_weight = periodic_operator_seam_weight_multiplier* &
                  sqrt(real(bulk_count, dp)/real(point_count - bulk_count, dp))
    weights(bulk_count + 1_i32:point_count) = seam_weight
  end subroutine set_periodic_check_weights

  subroutine apply_field_check_weights(matrix, weights)
    real(dp), intent(inout) :: matrix(:, :)
    real(dp), intent(in) :: weights(:)
    integer(i32) :: check_idx, ncheck

    ncheck = int(size(weights), i32)
    if (size(matrix, 1) /= 3_i32*ncheck) error stop 'field check-weight dimension mismatch.'
    do check_idx = 1_i32, ncheck
      matrix(check_idx, :) = weights(check_idx)*matrix(check_idx, :)
      matrix(ncheck + check_idx, :) = weights(check_idx)*matrix(ncheck + check_idx, :)
      matrix(2_i32*ncheck + check_idx, :) = weights(check_idx)*matrix(2_i32*ncheck + check_idx, :)
    end do
  end subroutine apply_field_check_weights

  subroutine build_proxy_multipole_matrix(plan, source_center, proxy_points, matrix)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: source_center(3), proxy_points(:, :)
    real(dp), intent(out) :: matrix(:, :)
    integer(i32) :: proxy_idx, beta_idx
    real(dp) :: d(3)
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    matrix = 0.0d0
    do proxy_idx = 1_i32, int(size(proxy_points, 2), i32)
      d = proxy_points(:, proxy_idx) - source_center
      call build_axis_powers(d, plan%options%order, xpow, ypow, zpow)
      do beta_idx = 1_i32, plan%ncoef
        matrix(beta_idx, proxy_idx) = xpow(plan%alpha(1, beta_idx))*ypow(plan%alpha(2, beta_idx)) &
                                      *zpow(plan%alpha(3, beta_idx))/plan%alpha_factorial(beta_idx)
      end do
    end do
  end subroutine build_proxy_multipole_matrix

  subroutine build_local_field_matrix(plan, target_center, check_points, matrix)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: target_center(3), check_points(:, :)
    real(dp), intent(out) :: matrix(:, :)
    integer(i32) :: check_idx, term_idx, coeff_idx, ncheck
    real(dp) :: d(3), monomial
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    matrix = 0.0d0
    ncheck = int(size(check_points, 2), i32)
    do check_idx = 1_i32, ncheck
      d = check_points(:, check_idx) - target_center
      call build_axis_powers(d, plan%options%order, xpow, ypow, zpow)
      do term_idx = 1_i32, plan%eval_term_count
        monomial = xpow(plan%eval_exp(1, term_idx))*ypow(plan%eval_exp(2, term_idx))*zpow(plan%eval_exp(3, term_idx))* &
                   plan%eval_inv_factorial(term_idx)
        coeff_idx = plan%eval_deriv_idx(1, term_idx)
        if (coeff_idx > 1_i32) matrix(check_idx, coeff_idx - 1_i32) = matrix(check_idx, coeff_idx - 1_i32) - monomial
        coeff_idx = plan%eval_deriv_idx(2, term_idx)
        if (coeff_idx > 1_i32) then
          matrix(ncheck + check_idx, coeff_idx - 1_i32) = matrix(ncheck + check_idx, coeff_idx - 1_i32) - monomial
        end if
        coeff_idx = plan%eval_deriv_idx(3, term_idx)
        if (coeff_idx > 1_i32) matrix(2_i32*ncheck + check_idx, coeff_idx - 1_i32) = &
          matrix(2_i32*ncheck + check_idx, coeff_idx - 1_i32) - monomial
      end do
    end do
  end subroutine build_local_field_matrix

  !> The constant fit is a weighted mean; cache its basis moments instead of rebuilding powers for every RHS.
  subroutine prepare_potential_mean_basis(plan, points, weights, mean_basis)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: points(:, :), weights(:)
    real(dp), intent(out) :: mean_basis(:)
    integer(i32) :: idx, alpha_idx
    real(dp) :: xpow(0:plan%options%order), ypow(0:plan%options%order), zpow(0:plan%options%order)

    mean_basis = 0.0_dp
    do idx = 1_i32, size(weights)
      call build_axis_powers(points(:, idx), plan%options%order, xpow, ypow, zpow)
      do alpha_idx = 1_i32, plan%ncoef
        mean_basis(alpha_idx) = mean_basis(alpha_idx) + weights(idx)* &
                                xpow(plan%alpha(1, alpha_idx))*ypow(plan%alpha(2, alpha_idx))*zpow(plan%alpha(3, alpha_idx))/ &
                                plan%alpha_factorial(alpha_idx)
      end do
    end do
  end subroutine prepare_potential_mean_basis

  subroutine build_minimum_norm_pseudoinverse(matrix, pinv)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), intent(out) :: pinv(:, :)
    integer(i32) :: nrow, ncol, rhs_idx
    real(dp), allocatable :: matrix_t(:, :), q(:, :), r(:, :), z(:), rhs(:), row_scale(:)

    nrow = int(size(matrix, 1), i32)
    ncol = int(size(matrix, 2), i32)
    if (size(pinv, 1) /= ncol .or. size(pinv, 2) /= nrow) error stop 'build_minimum_norm_pseudoinverse dimension mismatch.'

    allocate (matrix_t(ncol, nrow), q(ncol, nrow), r(nrow, nrow), z(nrow), rhs(nrow), row_scale(nrow))
    matrix_t = transpose(matrix)
    ! Multipole rows carry different powers of length. Normalize before QR so
    ! its rank tolerance does not erase high moments in micrometre SI geometry.
    ! Rescale the constraint RHS as well to preserve the original inverse map.
    do rhs_idx = 1_i32, nrow
      row_scale(rhs_idx) = sqrt(sum(matrix_t(:, rhs_idx)**2))
      if (row_scale(rhs_idx) <= tiny(1.0_dp)) row_scale(rhs_idx) = 1.0_dp
      matrix_t(:, rhs_idx) = matrix_t(:, rhs_idx)/row_scale(rhs_idx)
    end do
    q = 0.0d0
    r = 0.0d0
    call factor_tall_matrix_qr(matrix_t, q, r)

    do rhs_idx = 1_i32, nrow
      rhs = 0.0d0
      rhs(rhs_idx) = 1.0_dp/row_scale(rhs_idx)
      call solve_lower_triangular_transpose_system(r, rhs, z)
      pinv(:, rhs_idx) = matmul(q, z)
    end do
  end subroutine build_minimum_norm_pseudoinverse

  subroutine factor_tall_matrix_qr(matrix, q, r)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), intent(out) :: q(:, :), r(:, :)
    integer(i32) :: mrow, ncol, col_idx, basis_idx
    real(dp), allocatable :: v(:)
    real(dp) :: norm_v, corr, base_norm

    mrow = int(size(matrix, 1), i32)
    ncol = int(size(matrix, 2), i32)
    if (size(q, 1) /= mrow .or. size(q, 2) /= ncol) error stop 'factor_tall_matrix_qr q dimension mismatch.'
    if (size(r, 1) /= ncol .or. size(r, 2) /= ncol) error stop 'factor_tall_matrix_qr r dimension mismatch.'

    q = 0.0d0
    r = 0.0d0
    allocate (v(mrow))
    do col_idx = 1_i32, ncol
      v = matrix(:, col_idx)
      base_norm = max(sqrt(sum(v*v)), 1.0d0)
      do basis_idx = 1_i32, col_idx - 1_i32
        r(basis_idx, col_idx) = dot_product(q(:, basis_idx), v)
        v = v - r(basis_idx, col_idx)*q(:, basis_idx)
      end do
      do basis_idx = 1_i32, col_idx - 1_i32
        corr = dot_product(q(:, basis_idx), v)
        r(basis_idx, col_idx) = r(basis_idx, col_idx) + corr
        v = v - corr*q(:, basis_idx)
      end do
      norm_v = sqrt(sum(v*v))
      if (norm_v <= periodic_operator_qr_tol*base_norm) then
        r(col_idx, col_idx) = periodic_operator_qr_tol*base_norm
      else
        r(col_idx, col_idx) = norm_v
        q(:, col_idx) = v/norm_v
      end if
    end do
  end subroutine factor_tall_matrix_qr

  subroutine solve_lower_triangular_transpose_system(matrix, rhs, solution)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), intent(in) :: rhs(:)
    real(dp), intent(out) :: solution(:)
    integer(i32) :: ncol, row_idx, col_idx
    real(dp) :: diag_val

    ncol = int(size(matrix, 1), i32)
    if (size(matrix, 2) /= ncol .or. size(rhs) /= ncol .or. size(solution) /= ncol) then
      error stop 'solve_lower_triangular_transpose_system dimension mismatch.'
    end if

    solution = rhs
    do row_idx = 1_i32, ncol
      do col_idx = 1_i32, row_idx - 1_i32
        solution(row_idx) = solution(row_idx) - matrix(col_idx, row_idx)*solution(col_idx)
      end do
      diag_val = matrix(row_idx, row_idx)
      if (abs(diag_val) <= tiny(1.0d0)) then
        diag_val = sign(periodic_operator_qr_tol, diag_val + periodic_operator_qr_tol)
      end if
      solution(row_idx) = solution(row_idx)/diag_val
    end do
  end subroutine solve_lower_triangular_transpose_system

end module bem_coulomb_fmm_periodic_root_ops
