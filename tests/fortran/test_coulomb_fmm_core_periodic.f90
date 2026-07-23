!> periodic2 FMM コアの有限画像精度テスト。
program test_coulomb_fmm_core_periodic
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, build_plan, update_state, eval_point, &
                                  destroy_plan, destroy_state
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true
  implicit none

  call test_init(2)

  call test_begin('periodic2_field_accuracy')
  call test_periodic2_field_accuracy()
  call test_end()

  call test_begin('periodic2_canonical_seam_accuracy')
  call test_periodic2_canonical_seam_accuracy()
  call test_end()

  call test_summary()

contains

  subroutine test_periodic2_field_accuracy()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 6), e_fmm(3), e_ref(3)
    real(dp) :: norm_ref, rel_err, max_rel_err
    integer(i32) :: i, valid_count

    allocate (plan, state)
    call make_periodic_sources(src_pos, q)
    options%theta = 0.55_dp
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%use_periodic2 = .true.
    options%periodic_axes = [1_i32, 2_i32]
    options%periodic_len = [1.0_dp, 1.0_dp]
    options%periodic_image_layers = 1_i32
    options%target_box_min = [0.0_dp, 0.0_dp, -1.0_dp]
    options%target_box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    queries(:, 1) = [0.15_dp, 0.15_dp, -0.60_dp]
    queries(:, 2) = [0.85_dp, 0.20_dp, -0.20_dp]
    queries(:, 3) = [0.20_dp, 0.80_dp, 0.10_dp]
    queries(:, 4) = [0.75_dp, 0.75_dp, 0.50_dp]
    queries(:, 5) = [0.55_dp, 0.35_dp, -0.75_dp]
    queries(:, 6) = [0.35_dp, 0.60_dp, 0.85_dp]

    max_rel_err = 0.0_dp
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      call eval_point(plan, state, queries(:, i), e_fmm)
      call direct_field_periodic2( &
        plan%src_pos, q, queries(:, i), options%target_box_min, options%target_box_max, &
        options%periodic_axes, options%periodic_image_layers, e_ref &
        )
      norm_ref = sqrt(sum(e_ref*e_ref))
      if (norm_ref <= 1.0e-16_dp) cycle
      rel_err = sqrt(sum((e_fmm - e_ref)*(e_fmm - e_ref)))/norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 6_i32, 'periodic2 core test lost valid samples')
    write (*, '(A,I0,A,ES12.5)') &
      'test_periodic2_field_accuracy: valid_count=', valid_count, ', max_rel_err=', max_rel_err
    call assert_true(max_rel_err <= 8.0e-2_dp, 'periodic2 core FMM relative error exceeds 8e-2')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_periodic2_field_accuracy

  subroutine test_periodic2_canonical_seam_accuracy()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 4), e_fmm(3), e_ref(3)
    real(dp) :: norm_ref, rel_err, max_rel_err, bb_width_x, bb_width_y
    integer(i32) :: i, valid_count

    allocate (plan, state, src_pos(3, 6), q(6))
    src_pos(:, 1) = [0.02_dp, 0.03_dp, -0.30_dp]
    src_pos(:, 2) = [0.04_dp, 0.97_dp, 0.20_dp]
    src_pos(:, 3) = [0.96_dp, 0.02_dp, -0.10_dp]
    src_pos(:, 4) = [0.98_dp, 0.98_dp, 0.40_dp]
    src_pos(:, 5) = [0.50_dp, 0.50_dp, 0.00_dp]
    src_pos(:, 6) = [0.03_dp, 0.95_dp, -0.50_dp]
    q = [-1.0e-12_dp, 1.1e-12_dp, -0.9e-12_dp, 1.2e-12_dp, -1.3e-12_dp, 1.0e-12_dp]

    options%theta = 0.55_dp
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%use_periodic2 = .true.
    options%periodic_axes = [1_i32, 2_i32]
    options%periodic_len = [1.0_dp, 1.0_dp]
    options%periodic_image_layers = 1_i32
    options%periodic_far_correction = 'none'
    options%target_box_min = [0.0_dp, 0.0_dp, -1.0_dp]
    options%target_box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    bb_width_x = maxval(plan%src_pos(1, :)) - minval(plan%src_pos(1, :))
    bb_width_y = maxval(plan%src_pos(2, :)) - minval(plan%src_pos(2, :))
    call assert_true(bb_width_x < 0.60_dp, 'canonical seam should compact x bounding box')
    call assert_true(bb_width_y < 0.60_dp, 'canonical seam should compact y bounding box')

    queries(:, 1) = [0.15_dp, 0.15_dp, -0.50_dp]
    queries(:, 2) = [0.85_dp, 0.85_dp, 0.30_dp]
    queries(:, 3) = [0.50_dp, 0.50_dp, -0.80_dp]
    queries(:, 4) = [0.30_dp, 0.70_dp, 0.60_dp]

    max_rel_err = 0.0_dp
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      call eval_point(plan, state, queries(:, i), e_fmm)
      call direct_field_periodic2( &
        plan%src_pos, q, queries(:, i), options%target_box_min, options%target_box_max, &
        options%periodic_axes, options%periodic_image_layers, e_ref &
        )
      norm_ref = sqrt(sum(e_ref*e_ref))
      if (norm_ref <= 1.0e-16_dp) cycle
      rel_err = sqrt(sum((e_fmm - e_ref)*(e_fmm - e_ref)))/norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count >= 3_i32, 'canonical seam test lost valid samples')
    write (*, '(A,I0,A,ES12.5)') &
      'test_periodic2_canonical_seam_accuracy: valid_count=', valid_count, ', max_rel_err=', max_rel_err
    call assert_true(max_rel_err <= 8.0e-2_dp, 'canonical seam FMM relative error exceeds 8e-2')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_periodic2_canonical_seam_accuracy

  subroutine make_periodic_sources(src_pos, q)
    real(dp), allocatable, intent(out) :: src_pos(:, :)
    real(dp), allocatable, intent(out) :: q(:)

    allocate (src_pos(3, 8), q(8))
    src_pos(:, 1) = [0.20_dp, 0.20_dp, -0.40_dp]
    src_pos(:, 2) = [0.80_dp, 0.20_dp, -0.20_dp]
    src_pos(:, 3) = [0.25_dp, 0.75_dp, 0.10_dp]
    src_pos(:, 4) = [0.75_dp, 0.70_dp, 0.35_dp]
    src_pos(:, 5) = [0.40_dp, 0.35_dp, -0.70_dp]
    src_pos(:, 6) = [0.60_dp, 0.55_dp, 0.55_dp]
    src_pos(:, 7) = [0.30_dp, 0.60_dp, 0.85_dp]
    src_pos(:, 8) = [0.65_dp, 0.30_dp, -0.85_dp]
    q = [ &
        -1.0e-12_dp, 1.1e-12_dp, -1.2e-12_dp, 1.3e-12_dp, &
        -1.4e-12_dp, 1.5e-12_dp, -1.6e-12_dp, 1.3e-12_dp &
        ]
  end subroutine make_periodic_sources

  subroutine direct_field_periodic2(src_pos, q, r, box_min, box_max, periodic_axes, nimg, e)
    real(dp), intent(in) :: src_pos(:, :), q(:), r(3), box_min(3), box_max(3)
    integer(i32), intent(in) :: periodic_axes(2), nimg
    real(dp), intent(out) :: e(3)
    integer(i32) :: idx, img_i, img_j, axis1, axis2
    real(dp) :: shifted(3), d(3), r2, inv_r3, len1, len2

    axis1 = periodic_axes(1)
    axis2 = periodic_axes(2)
    len1 = box_max(axis1) - box_min(axis1)
    len2 = box_max(axis2) - box_min(axis2)
    e = 0.0_dp
    do idx = 1_i32, int(size(q), i32)
      do img_i = -nimg, nimg
        do img_j = -nimg, nimg
          shifted = src_pos(:, idx)
          shifted(axis1) = shifted(axis1) + real(img_i, dp)*len1
          shifted(axis2) = shifted(axis2) + real(img_j, dp)*len2
          d = r - shifted
          r2 = sum(d*d)
          if (r2 <= tiny(1.0_dp)) cycle
          inv_r3 = 1.0_dp/(sqrt(r2)*r2)
          e = e + q(idx)*inv_r3*d
        end do
      end do
    end do
  end subroutine direct_field_periodic2

end program test_coulomb_fmm_core_periodic
