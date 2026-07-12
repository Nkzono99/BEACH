!> periodic2 root operator の前計算。
module bem_coulomb_fmm_periodic_root_ops
!$ use omp_lib, only: omp_get_num_threads
  use bem_version, only: beach_version
  use bem_kinds, only: dp, i32
  use bem_mpi, only: mpi_context, mpi_initialize, mpi_shutdown, mpi_is_root, mpi_get_rank_size, &
                     mpi_allreduce_sum_real_dp_array, mpi_bcast_i32_array, mpi_bcast_real_dp_array
  use bem_filesystem, only: create_directories, filesystem_success, acquire_file_lock, release_file_lock
  use bem_coulomb_fmm_types, only: fmm_plan_type
  use bem_coulomb_fmm_basis, only: build_axis_powers
  use bem_coulomb_fmm_periodic, only: use_periodic2_m2l_root_oracle, use_periodic2_cached_kneq0, &
                                      use_periodic2_root_operator
  use bem_coulomb_fmm_periodic_ewald, only: add_periodic2_exact_ewald_correction_single_source, &
                                            add_periodic2_exact_ewald_potential_correction_single_source
  use bem_coulomb_fmm_periodic_cache, only: read_periodic_operator_cache, write_periodic_operator_cache, &
                                            periodic_cache_fingerprint, periodic_cache_ok
  use bem_regularized_qr, only: regularized_qr_type, prepare_regularized_qr, solve_regularized_qr
  use bem_coulomb_fmm_tree_utils, only: active_tree_nnode, active_tree_max_depth, active_tree_child_count, &
                                        active_tree_child_idx, active_tree_node_center, active_tree_node_half_size
  implicit none
  private

  integer(i32), parameter :: root_oracle_target_depth = 1_i32
  real(dp), parameter :: root_oracle_proxy_multiplier = 8.0d0
  real(dp), parameter :: root_oracle_check_multiplier = 24.0d0
  real(dp), parameter :: root_oracle_proxy_shell_scale = 1.15d0
  real(dp), parameter :: root_oracle_check_shell_scale = 0.92d0
  real(dp), parameter :: root_oracle_tall_box_ratio = 4.0d0
  real(dp), parameter :: root_oracle_lstsq_ridge = 1.0d-12
  real(dp), parameter :: root_oracle_qr_tol = 1.0d-12

  public :: precompute_periodic_root_operator

contains

  !> periodic2 の root operator を前計算する。
  !! @param[inout] plan FMM 計画。
  subroutine precompute_periodic_root_operator(plan)
    type(fmm_plan_type), intent(inout) :: plan
    type(mpi_context) :: mpi

    if (allocated(plan%periodic_root_target_nodes)) deallocate (plan%periodic_root_target_nodes)
    if (allocated(plan%periodic_root_operator)) deallocate (plan%periodic_root_operator)
    plan%periodic_root_operator_ready = .false.
    plan%periodic_root_target_count = 0_i32

    plan%periodic_cache_hit = .false.
    plan%periodic_cache_fingerprint = ''
    plan%periodic_cache_path = ''
    if (use_periodic2_cached_kneq0(plan)) then
      call mpi_initialize(mpi)
      call precompute_periodic_cached_operator_collective(plan, mpi)
      call mpi_shutdown(mpi)
    else if (use_periodic2_root_operator(plan)) then
      call precompute_periodic_root_oracle_operator(plan)
    end if
  end subroutine precompute_periodic_root_operator

  subroutine precompute_periodic_root_oracle_operator(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: nproxy, ncheck, j, i, target_idx, node_idx, n_target_nodes, anchor_depth, target_count
    real(dp) :: source_center(3), source_half(3), proxy_half(3), target_center(3), target_half(3)
    real(dp), allocatable :: proxy_points(:, :), check_points(:, :)
    real(dp), allocatable :: proxy_to_multipole(:, :), proxy_to_local(:, :)
    real(dp), allocatable :: field_matrix(:, :), field_rhs(:)
    real(dp), allocatable :: coeff(:), proxy_pinv(:, :)
    real(dp), allocatable :: potential_rhs(:)
    real(dp) :: e_res(3), phi_res
    integer(i32), allocatable :: target_nodes(:)
    type(regularized_qr_type) :: field_factorization
    logical :: use_target_tree, cache_loaded
    integer :: cache_lock_fd, lock_status, directory_status

    if (.not. plan%periodic_ewald%ready) return
    if (plan%ncoef <= 1_i32) return
    use_target_tree = plan%target_tree_ready
    n_target_nodes = active_tree_nnode(plan, use_target_tree)
    if (n_target_nodes <= 0_i32) return

    nproxy = max(4_i32*plan%ncoef, int(root_oracle_proxy_multiplier*real(plan%ncoef, dp), i32))
    ncheck = max(8_i32*plan%ncoef, int(root_oracle_check_multiplier*real(plan%ncoef, dp), i32))
    source_center = plan%node_center(:, 1_i32)
    source_half = plan%node_half_size(:, 1_i32)
    proxy_half = source_half
    proxy_half = max(proxy_half, 0.25d0*min(plan%options%periodic_len(1), plan%options%periodic_len(2)))
    anchor_depth = periodic_root_anchor_depth(plan, use_target_tree)
    allocate (target_nodes(max(1_i32, n_target_nodes)))
    target_count = 0_i32
    call collect_periodic_root_targets(plan, use_target_tree, 1_i32, anchor_depth, target_nodes, target_count)
    if (target_count <= 0_i32) then
      deallocate (target_nodes)
      return
    end if

    cache_loaded = .false.
    cache_lock_fd = -1
    if (use_periodic2_cached_kneq0(plan)) then
      call try_load_cached_operator(plan, use_target_tree, target_nodes, target_count, cache_loaded)
      if (cache_loaded) then
        deallocate (target_nodes)
        return
      end if
      call create_directories(trim(plan%options%periodic_cache_dir), directory_status)
      if (directory_status /= filesystem_success) error stop 'Failed to create periodic operator cache directory.'
      call acquire_file_lock(trim(plan%periodic_cache_path)//'.lock', cache_lock_fd, lock_status)
      if (lock_status /= filesystem_success) error stop 'Failed to lock periodic operator cache.'
      call try_load_cached_operator(plan, use_target_tree, target_nodes, target_count, cache_loaded)
      if (cache_loaded) then
        call release_file_lock(cache_lock_fd, lock_status)
        deallocate (target_nodes)
        return
      end if
    end if

    allocate (proxy_points(3, nproxy), check_points(3, ncheck))
    call build_root_surface_points(source_center, proxy_half, nproxy, 0.13d0, root_oracle_proxy_shell_scale, proxy_points)

    allocate (proxy_to_multipole(plan%ncoef, nproxy), proxy_to_local(plan%ncoef, nproxy))
    allocate (field_matrix(3_i32*ncheck, plan%ncoef - 1_i32))
    allocate (field_rhs(3_i32*ncheck), coeff(plan%ncoef - 1_i32))
    if (use_periodic2_cached_kneq0(plan)) allocate (potential_rhs(ncheck))
    allocate (proxy_pinv(nproxy, plan%ncoef))

    call build_proxy_multipole_matrix(plan, source_center, proxy_points, proxy_to_multipole)
    call build_minimum_norm_pseudoinverse(proxy_to_multipole, proxy_pinv)

    plan%periodic_root_target_count = target_count
    allocate (plan%periodic_root_target_nodes(target_count))
    plan%periodic_root_target_nodes = target_nodes(1:target_count)
    deallocate (target_nodes)
    allocate (plan%periodic_root_operator(plan%ncoef, plan%ncoef, target_count))
    plan%periodic_root_operator = 0.0d0

    do target_idx = 1_i32, plan%periodic_root_target_count
      node_idx = plan%periodic_root_target_nodes(target_idx)
      target_center = active_tree_node_center(plan, use_target_tree, node_idx)
      target_half = active_tree_node_half_size(plan, use_target_tree, node_idx)
      call build_root_surface_points(target_center, target_half, ncheck, 0.37d0, root_oracle_check_shell_scale, check_points)
      call build_local_field_matrix(plan, target_center, check_points, field_matrix)
      call prepare_regularized_qr( &
        field_factorization, field_matrix, root_oracle_lstsq_ridge, root_oracle_qr_tol &
        )
      plan%periodic_qr_preparation_count = plan%periodic_qr_preparation_count + 1_i32

      proxy_to_local = 0.0d0
      do j = 1_i32, nproxy
        field_rhs = 0.0d0
        if (use_periodic2_cached_kneq0(plan)) potential_rhs = 0.0_dp
        do i = 1_i32, ncheck
          e_res = 0.0d0
          if (use_periodic2_cached_kneq0(plan)) then
            phi_res = 0.0_dp
            call add_periodic2_exact_ewald_correction_single_source( &
              plan, 1.0d0, proxy_points(:, j), check_points(:, i), e_res &
              )
            call add_periodic2_exact_ewald_potential_correction_single_source( &
              plan, 1.0_dp, proxy_points(:, j), check_points(:, i), phi_res &
              )
            potential_rhs(i) = phi_res
          else
            call add_periodic2_exact_ewald_correction_single_source( &
              plan, 1.0d0, proxy_points(:, j), check_points(:, i), e_res &
              )
          end if
          field_rhs(i) = e_res(1)
          field_rhs(ncheck + i) = e_res(2)
          field_rhs(2_i32*ncheck + i) = e_res(3)
        end do
        call solve_regularized_qr(field_factorization, field_rhs, coeff)
        proxy_to_local(2:plan%ncoef, j) = coeff
        if (use_periodic2_cached_kneq0(plan)) then
          call fit_local_potential_constant( &
            plan, target_center, check_points, potential_rhs, proxy_to_local(:, j) &
            )
        end if
      end do

      plan%periodic_root_operator(:, :, target_idx) = matmul(proxy_to_local, proxy_pinv)
      if (.not. use_periodic2_cached_kneq0(plan)) plan%periodic_root_operator(1_i32, :, target_idx) = 0.0d0
    end do

    plan%periodic_root_operator_ready = .true.
    plan%periodic_operator_local_target_count = plan%periodic_root_target_count
    plan%periodic_operator_build_count = plan%periodic_operator_build_count + 1_i32
    if (use_periodic2_cached_kneq0(plan)) then
      call publish_cached_operator(plan)
      call release_file_lock(cache_lock_fd, lock_status)
      if (lock_status /= filesystem_success) error stop 'Failed to release periodic operator cache lock.'
    end if
  end subroutine precompute_periodic_root_oracle_operator

  subroutine precompute_periodic_cached_operator_collective(plan, mpi)
    type(fmm_plan_type), intent(inout) :: plan
    type(mpi_context), intent(in) :: mpi
    integer(i32) :: nproxy, ncheck, j, i, target_idx, node_idx, n_target_nodes, anchor_depth, target_count
    integer(i32) :: rank, rank_count
    real(dp) :: source_center(3), source_half(3), proxy_half(3), target_center(3), target_half(3)
    real(dp), allocatable :: proxy_points(:, :), check_points(:, :)
    real(dp), allocatable :: proxy_to_multipole(:, :), proxy_to_local(:, :), proxy_pinv(:, :)
    real(dp), allocatable :: field_matrix(:, :), field_rhs(:), coeff(:), potential_rhs(:), flat_operator(:)
    real(dp) :: e_res(3), phi_res
    integer(i32), allocatable :: target_nodes(:)
    type(regularized_qr_type) :: field_factorization
    logical :: use_target_tree, cache_loaded
    integer :: cache_lock_fd, lock_status, directory_status

    if (.not. plan%periodic_ewald%ready) return
    if (plan%ncoef <= 1_i32) return
    call mpi_get_rank_size(rank, rank_count, mpi)
    use_target_tree = plan%target_tree_ready
    n_target_nodes = active_tree_nnode(plan, use_target_tree)
    if (n_target_nodes <= 0_i32) return

    nproxy = max(4_i32*plan%ncoef, int(root_oracle_proxy_multiplier*real(plan%ncoef, dp), i32))
    ncheck = max(8_i32*plan%ncoef, int(root_oracle_check_multiplier*real(plan%ncoef, dp), i32))
    source_center = plan%node_center(:, 1_i32)
    source_half = plan%node_half_size(:, 1_i32)
    proxy_half = source_half
    proxy_half = max(proxy_half, 0.25_dp*min(plan%options%periodic_len(1), plan%options%periodic_len(2)))
    anchor_depth = periodic_root_anchor_depth(plan, use_target_tree)
    allocate (target_nodes(max(1_i32, n_target_nodes)))
    target_count = 0_i32
    call collect_periodic_root_targets(plan, use_target_tree, 1_i32, anchor_depth, target_nodes, target_count)
    if (target_count <= 0_i32) then
      deallocate (target_nodes)
      return
    end if

    call set_cache_identity(plan, use_target_tree, target_nodes, target_count)
    cache_loaded = .false.
    cache_lock_fd = -1
    if (mpi_is_root(mpi)) call try_load_cached_operator(plan, use_target_tree, target_nodes, target_count, cache_loaded)
    call share_cached_operator(plan, mpi, cache_loaded)
    if (cache_loaded) then
      deallocate (target_nodes)
      return
    end if

    if (mpi_is_root(mpi)) then
      call create_directories(trim(plan%options%periodic_cache_dir), directory_status)
      if (directory_status /= filesystem_success) error stop 'Failed to create periodic operator cache directory.'
      call acquire_file_lock(trim(plan%periodic_cache_path)//'.lock', cache_lock_fd, lock_status)
      if (lock_status /= filesystem_success) error stop 'Failed to lock periodic operator cache.'
      call try_load_cached_operator(plan, use_target_tree, target_nodes, target_count, cache_loaded)
      if (cache_loaded) then
        call release_file_lock(cache_lock_fd, lock_status)
        if (lock_status /= filesystem_success) error stop 'Failed to release periodic operator cache lock.'
      end if
    end if
    call share_cached_operator(plan, mpi, cache_loaded)
    if (cache_loaded) then
      deallocate (target_nodes)
      return
    end if

    allocate (proxy_points(3, nproxy), check_points(3, ncheck))
    call build_root_surface_points(source_center, proxy_half, nproxy, 0.13_dp, root_oracle_proxy_shell_scale, proxy_points)
    allocate (proxy_to_multipole(plan%ncoef, nproxy), proxy_to_local(plan%ncoef, nproxy))
    allocate (field_matrix(3_i32*ncheck, plan%ncoef - 1_i32), proxy_pinv(nproxy, plan%ncoef))
    call build_proxy_multipole_matrix(plan, source_center, proxy_points, proxy_to_multipole)
    call build_minimum_norm_pseudoinverse(proxy_to_multipole, proxy_pinv)

    plan%periodic_root_target_count = target_count
    allocate (plan%periodic_root_target_nodes(target_count))
    plan%periodic_root_target_nodes = target_nodes(1:target_count)
    deallocate (target_nodes)
    allocate (plan%periodic_root_operator(plan%ncoef, plan%ncoef, target_count))
    plan%periodic_root_operator = 0.0_dp

    do target_idx = rank + 1_i32, plan%periodic_root_target_count, rank_count
      node_idx = plan%periodic_root_target_nodes(target_idx)
      target_center = active_tree_node_center(plan, use_target_tree, node_idx)
      target_half = active_tree_node_half_size(plan, use_target_tree, node_idx)
      call build_root_surface_points(target_center, target_half, ncheck, 0.37_dp, root_oracle_check_shell_scale, check_points)
      call build_local_field_matrix(plan, target_center, check_points, field_matrix)
      call prepare_regularized_qr( &
        field_factorization, field_matrix, root_oracle_lstsq_ridge, root_oracle_qr_tol &
        )
      plan%periodic_qr_preparation_count = plan%periodic_qr_preparation_count + 1_i32
      plan%periodic_operator_local_target_count = plan%periodic_operator_local_target_count + 1_i32

      proxy_to_local = 0.0_dp
      plan%periodic_operator_thread_count = max(plan%periodic_operator_thread_count, 1_i32)
!$omp parallel default(shared) private(j, i, field_rhs, potential_rhs, coeff, e_res, phi_res)
      allocate (field_rhs(3_i32*ncheck), potential_rhs(ncheck), coeff(plan%ncoef - 1_i32))
!$omp single
!$    plan%periodic_operator_thread_count = max( &
!$                                          plan%periodic_operator_thread_count, int(omp_get_num_threads(), i32) &
!$                                          )
!$omp end single
!$omp do schedule(static)
      do j = 1_i32, nproxy
        field_rhs = 0.0_dp
        potential_rhs = 0.0_dp
        do i = 1_i32, ncheck
          e_res = 0.0_dp
          phi_res = 0.0_dp
          call add_periodic2_exact_ewald_correction_single_source( &
            plan, 1.0_dp, proxy_points(:, j), check_points(:, i), e_res &
            )
          call add_periodic2_exact_ewald_potential_correction_single_source( &
            plan, 1.0_dp, proxy_points(:, j), check_points(:, i), phi_res &
            )
          potential_rhs(i) = phi_res
          field_rhs(i) = e_res(1)
          field_rhs(ncheck + i) = e_res(2)
          field_rhs(2_i32*ncheck + i) = e_res(3)
        end do
        call solve_regularized_qr(field_factorization, field_rhs, coeff)
        proxy_to_local(2:plan%ncoef, j) = coeff
        call fit_local_potential_constant( &
          plan, target_center, check_points, potential_rhs, proxy_to_local(:, j) &
          )
      end do
!$omp end do
      deallocate (field_rhs, potential_rhs, coeff)
!$omp end parallel
      plan%periodic_root_operator(:, :, target_idx) = matmul(proxy_to_local, proxy_pinv)
    end do

    allocate (flat_operator(size(plan%periodic_root_operator)))
    flat_operator = reshape(plan%periodic_root_operator, [size(flat_operator)])
    call mpi_allreduce_sum_real_dp_array(mpi, flat_operator)
    plan%periodic_root_operator = reshape(flat_operator, shape(plan%periodic_root_operator))
    deallocate (flat_operator)
    plan%periodic_root_operator_ready = .true.
    if (mpi_is_root(mpi)) then
      plan%periodic_operator_build_count = plan%periodic_operator_build_count + 1_i32
      call publish_cached_operator(plan)
      call release_file_lock(cache_lock_fd, lock_status)
      if (lock_status /= filesystem_success) error stop 'Failed to release periodic operator cache lock.'
    end if
  end subroutine precompute_periodic_cached_operator_collective

  subroutine share_cached_operator(plan, mpi, loaded)
    type(fmm_plan_type), intent(inout) :: plan
    type(mpi_context), intent(in) :: mpi
    logical, intent(inout) :: loaded
    integer(i32) :: header(2)
    real(dp), allocatable :: flat_operator(:)

    header = 0_i32
    if (mpi_is_root(mpi)) then
      header = [merge(1_i32, 0_i32, loaded), plan%periodic_root_target_count]
    end if
    call mpi_bcast_i32_array(mpi, header, 0_i32)
    loaded = header(1) == 1_i32
    if (.not. loaded) return
    if (header(2) <= 0_i32) error stop 'Cached periodic operator has no target nodes.'

    if (.not. mpi_is_root(mpi)) then
      plan%periodic_root_target_count = header(2)
      allocate (plan%periodic_root_target_nodes(header(2)))
      allocate (plan%periodic_root_operator(plan%ncoef, plan%ncoef, header(2)))
    end if
    call mpi_bcast_i32_array(mpi, plan%periodic_root_target_nodes, 0_i32)
    allocate (flat_operator(size(plan%periodic_root_operator)))
    if (mpi_is_root(mpi)) flat_operator = reshape(plan%periodic_root_operator, [size(flat_operator)])
    call mpi_bcast_real_dp_array(mpi, flat_operator, 0_i32)
    if (.not. mpi_is_root(mpi)) then
      plan%periodic_root_operator = reshape(flat_operator, shape(plan%periodic_root_operator))
      plan%periodic_root_operator_ready = .true.
      plan%periodic_cache_hit = .true.
    end if
    deallocate (flat_operator)
  end subroutine share_cached_operator

  subroutine try_load_cached_operator(plan, use_target_tree, target_nodes, target_count, loaded)
    type(fmm_plan_type), intent(inout) :: plan
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: target_nodes(:), target_count
    logical, intent(out) :: loaded
    integer(i32), allocatable :: cached_nodes(:)
    real(dp), allocatable :: cached_operator(:, :, :)
    integer(i32) :: cache_status

    loaded = .false.
    call set_cache_identity(plan, use_target_tree, target_nodes, target_count)
    call read_periodic_operator_cache( &
      trim(plan%periodic_cache_path), trim(plan%periodic_cache_fingerprint), cached_nodes, cached_operator, cache_status, &
      expected_ncoef=plan%ncoef, expected_ntarget=target_count &
      )
    if (cache_status /= periodic_cache_ok) return
    if (size(cached_nodes) /= target_count .or. any(cached_nodes /= target_nodes(1:target_count)) .or. &
        size(cached_operator, 1) /= plan%ncoef .or. size(cached_operator, 2) /= plan%ncoef .or. &
        size(cached_operator, 3) /= target_count) then
      deallocate (cached_nodes, cached_operator)
      return
    end if
    plan%periodic_root_target_count = target_count
    call move_alloc(cached_nodes, plan%periodic_root_target_nodes)
    call move_alloc(cached_operator, plan%periodic_root_operator)
    plan%periodic_root_operator_ready = .true.
    plan%periodic_cache_hit = .true.
    loaded = .true.
  end subroutine try_load_cached_operator

  subroutine set_cache_identity(plan, use_target_tree, target_nodes, target_count)
    type(fmm_plan_type), intent(inout) :: plan
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: target_nodes(:), target_count
    integer(i32), allocatable :: integer_values(:)
    real(dp), allocatable :: real_values(:)
    integer(i32) :: target_idx, node_idx, real_pos
    character(len=256) :: tag

    allocate (integer_values(8 + target_count), real_values(11 + 6*target_count))
    integer_values(1:8) = [ &
      plan%options%order, plan%ncoef, plan%options%periodic_axes, plan%options%periodic_image_layers, &
      plan%options%periodic_ewald_layers, merge(1_i32, 0_i32, plan%panel_source), target_count &
      ]
    integer_values(9:) = target_nodes(1:target_count)
    real_values(1:11) = [ &
      plan%options%periodic_len, plan%options%softening, plan%options%periodic_ewald_alpha, &
      plan%options%periodic_generation_tolerance, plan%node_center(:, 1_i32), plan%node_half_size(:, 1_i32) &
      ]
    real_pos = 12_i32
    do target_idx = 1_i32, target_count
      node_idx = target_nodes(target_idx)
      real_values(real_pos:real_pos + 2) = active_tree_node_center(plan, use_target_tree, node_idx)
      real_values(real_pos + 3:real_pos + 5) = active_tree_node_half_size(plan, use_target_tree, node_idx)
      real_pos = real_pos + 6_i32
    end do
    if (plan%panel_source) then
      tag = 'periodic-kneq0-generator-v7|real64|kernel=triangle_p0|zero=state_subtracted|normalization=core|'
    else
      tag = 'periodic-kneq0-generator-v7|real64|kernel=point_softened|zero=state_subtracted|normalization=core|'
    end if
    tag = trim(tag)//trim(beach_version)
    plan%periodic_cache_fingerprint = periodic_cache_fingerprint(tag, integer_values, real_values)
    plan%periodic_cache_path = trim(plan%options%periodic_cache_dir)//'/operator-'// &
                               trim(plan%periodic_cache_fingerprint)//'.bin'
  end subroutine set_cache_identity

  subroutine publish_cached_operator(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: cache_status
    integer :: directory_status

    call create_directories(trim(plan%options%periodic_cache_dir), directory_status)
    if (directory_status /= filesystem_success) error stop 'Failed to create periodic operator cache directory.'
    call write_periodic_operator_cache( &
      trim(plan%periodic_cache_path), trim(plan%periodic_cache_fingerprint), &
      plan%periodic_root_target_nodes, plan%periodic_root_operator, cache_status &
      )
    if (cache_status /= periodic_cache_ok) error stop 'Failed to publish periodic operator cache.'
  end subroutine publish_cached_operator

  pure integer(i32) function periodic_root_anchor_depth(plan, use_target_tree)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree
    real(dp) :: target_half(3), periodic_span, target_span_ratio

    periodic_root_anchor_depth = min(active_tree_max_depth(plan, use_target_tree), root_oracle_target_depth)
    target_half = active_tree_node_half_size(plan, use_target_tree, 1_i32)
    periodic_span = max(minval(plan%options%periodic_len), tiny(1.0d0))
    target_span_ratio = maxval(2.0d0*target_half)/periodic_span
    if (target_span_ratio > root_oracle_tall_box_ratio) then
      periodic_root_anchor_depth = min(active_tree_max_depth(plan, use_target_tree), periodic_root_anchor_depth + 1_i32)
    end if
  end function periodic_root_anchor_depth

  recursive subroutine collect_periodic_root_targets(plan, use_target_tree, node_idx, anchor_depth, target_nodes, target_count)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: node_idx, anchor_depth
    integer(i32), intent(inout) :: target_nodes(:)
    integer(i32), intent(inout) :: target_count
    integer(i32) :: node_depth, child_k, child_count

    if (node_idx <= 0_i32) return
    if (use_target_tree) then
      node_depth = plan%target_node_depth(node_idx)
    else
      node_depth = plan%node_depth(node_idx)
    end if
    child_count = active_tree_child_count(plan, use_target_tree, node_idx)
    if (child_count <= 0_i32 .or. node_depth >= anchor_depth) then
      target_count = target_count + 1_i32
      target_nodes(target_count) = node_idx
      return
    end if

    do child_k = 1_i32, child_count
      call collect_periodic_root_targets( &
        plan, use_target_tree, active_tree_child_idx(plan, use_target_tree, child_k, node_idx), &
        anchor_depth, target_nodes, target_count &
        )
    end do
  end subroutine collect_periodic_root_targets

  subroutine build_root_surface_points(center, half_size, npoint, offset, scale, points)
    real(dp), intent(in) :: center(3), half_size(3), offset, scale
    integer(i32), intent(in) :: npoint
    real(dp), intent(out) :: points(3, npoint)
    integer(i32) :: idx, face
    real(dp) :: f1, f2, u, v, h(3)
    real(dp), parameter :: g1 = 0.7548776662466927d0
    real(dp), parameter :: g2 = 0.5698402909980532d0

    h = scale*half_size
    do idx = 1_i32, npoint
      f1 = modulo(offset + real(idx, dp)*g1, 1.0d0)
      f2 = modulo(offset + real(idx, dp)*g2, 1.0d0)
      u = 2.0d0*(0.05d0 + 0.9d0*f1) - 1.0d0
      v = 2.0d0*(0.05d0 + 0.9d0*f2) - 1.0d0
      face = mod(idx - 1_i32, 6_i32) + 1_i32
      select case (face)
      case (1_i32)
        points(:, idx) = center + [h(1), u*h(2), v*h(3)]
      case (2_i32)
        points(:, idx) = center + [-h(1), u*h(2), v*h(3)]
      case (3_i32)
        points(:, idx) = center + [u*h(1), h(2), v*h(3)]
      case (4_i32)
        points(:, idx) = center + [u*h(1), -h(2), v*h(3)]
      case (5_i32)
        points(:, idx) = center + [u*h(1), v*h(2), h(3)]
      case default
        points(:, idx) = center + [u*h(1), v*h(2), -h(3)]
      end select
    end do
  end subroutine build_root_surface_points

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

  subroutine fit_local_potential_constant(plan, target_center, check_points, potential, local_coeff)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: target_center(3), check_points(:, :)
    real(dp), intent(in) :: potential(:)
    real(dp), intent(inout) :: local_coeff(:)
    integer(i32) :: check_idx, alpha_idx, ncheck
    real(dp) :: d(3), monomial, predicted
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    ncheck = int(size(check_points, 2), i32)
    local_coeff(1) = 0.0_dp
    do check_idx = 1_i32, ncheck
      d = check_points(:, check_idx) - target_center
      call build_axis_powers(d, plan%options%order, xpow, ypow, zpow)
      predicted = 0.0_dp
      do alpha_idx = 2_i32, plan%ncoef
        monomial = xpow(plan%alpha(1, alpha_idx))*ypow(plan%alpha(2, alpha_idx))* &
                   zpow(plan%alpha(3, alpha_idx))/plan%alpha_factorial(alpha_idx)
        predicted = predicted + local_coeff(alpha_idx)*monomial
      end do
      local_coeff(1) = local_coeff(1) + potential(check_idx) - predicted
    end do
    local_coeff(1) = local_coeff(1)/real(ncheck, dp)
  end subroutine fit_local_potential_constant

  subroutine build_minimum_norm_pseudoinverse(matrix, pinv)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), intent(out) :: pinv(:, :)
    integer(i32) :: nrow, ncol, rhs_idx
    real(dp), allocatable :: matrix_t(:, :), q(:, :), r(:, :), z(:), rhs(:)

    nrow = int(size(matrix, 1), i32)
    ncol = int(size(matrix, 2), i32)
    if (size(pinv, 1) /= ncol .or. size(pinv, 2) /= nrow) error stop 'build_minimum_norm_pseudoinverse dimension mismatch.'

    allocate (matrix_t(ncol, nrow), q(ncol, nrow), r(nrow, nrow), z(nrow), rhs(nrow))
    matrix_t = transpose(matrix)
    q = 0.0d0
    r = 0.0d0
    call factor_tall_matrix_qr(matrix_t, q, r)

    do rhs_idx = 1_i32, nrow
      rhs = 0.0d0
      rhs(rhs_idx) = 1.0d0
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
      if (norm_v <= root_oracle_qr_tol*base_norm) then
        r(col_idx, col_idx) = root_oracle_qr_tol*base_norm
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
      if (abs(diag_val) <= tiny(1.0d0)) diag_val = sign(root_oracle_qr_tol, diag_val + root_oracle_qr_tol)
      solution(row_idx) = solution(row_idx)/diag_val
    end do
  end subroutine solve_lower_triangular_transpose_system

  subroutine solve_square_system(matrix, rhs, solution)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), intent(in) :: rhs(:)
    real(dp), intent(out) :: solution(:)
    integer(i32) :: ncol, pivot_row, row_idx, col_idx, swap_idx
    real(dp), allocatable :: work(:, :), rhs_work(:), row_tmp(:)
    real(dp) :: pivot_abs, factor

    ncol = int(size(matrix, 1), i32)
    if (size(matrix, 2) /= ncol .or. size(rhs) /= ncol .or. size(solution) /= ncol) then
      error stop 'solve_square_system dimension mismatch.'
    end if

    allocate (work(ncol, ncol), rhs_work(ncol), row_tmp(ncol))
    work = matrix
    rhs_work = rhs

    do col_idx = 1_i32, ncol
      pivot_row = col_idx
      pivot_abs = abs(work(col_idx, col_idx))
      do row_idx = col_idx + 1_i32, ncol
        if (abs(work(row_idx, col_idx)) > pivot_abs) then
          pivot_abs = abs(work(row_idx, col_idx))
          pivot_row = row_idx
        end if
      end do
      if (pivot_abs <= 1.0d-20) error stop 'periodic root oracle linear system is singular.'

      if (pivot_row /= col_idx) then
        row_tmp = work(col_idx, :)
        work(col_idx, :) = work(pivot_row, :)
        work(pivot_row, :) = row_tmp
        factor = rhs_work(col_idx)
        rhs_work(col_idx) = rhs_work(pivot_row)
        rhs_work(pivot_row) = factor
      end if

      factor = work(col_idx, col_idx)
      work(col_idx, col_idx:ncol) = work(col_idx, col_idx:ncol)/factor
      rhs_work(col_idx) = rhs_work(col_idx)/factor
      do row_idx = col_idx + 1_i32, ncol
        factor = work(row_idx, col_idx)
        if (abs(factor) <= tiny(1.0d0)) cycle
        work(row_idx, col_idx:ncol) = work(row_idx, col_idx:ncol) - factor*work(col_idx, col_idx:ncol)
        rhs_work(row_idx) = rhs_work(row_idx) - factor*rhs_work(col_idx)
      end do
    end do

    solution = rhs_work
    do row_idx = ncol, 1_i32, -1_i32
      do swap_idx = row_idx + 1_i32, ncol
        solution(row_idx) = solution(row_idx) - work(row_idx, swap_idx)*solution(swap_idx)
      end do
    end do
  end subroutine solve_square_system

end module bem_coulomb_fmm_periodic_root_ops
