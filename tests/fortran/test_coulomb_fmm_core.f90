!> 疎結合 Coulomb FMM コアの単体挙動を検証する。
program test_coulomb_fmm_core
  use bem_constants, only: k_coulomb
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_templates, only: make_sphere
  use bem_field, only: electric_field_at
  use bem_field_solver, only: field_solver_type
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, build_plan, update_state, eval_point, &
                                   destroy_plan, destroy_state
  use bem_coulomb_fmm_periodic_ewald, only: add_periodic2_exact_ewald_correction_all_sources
  use test_support, only: assert_true, assert_close_dp, assert_allclose_1d, assert_equal_i32
  implicit none

  call test_p2m_m2m_root_moments()
  call test_free_field_accuracy()
  call test_softened_free_field_accuracy()
  call test_periodic2_field_accuracy()
  call test_periodic2_m2l_root_trunc_correction_effect()
  call test_periodic2_m2l_root_oracle_correction_effect()
  call test_periodic2_m2l_root_trunc_charged_wall_closure()
  call test_target_box_dual_tree()
  call test_state_update_reuse()
  call test_state_eval_profile_counts()
  call test_field_solver_core_adapter()
  call test_field_solver_core_softened_adapter()
  call test_field_solver_core_periodic2_m2l_root_trunc_adapter()

contains

  subroutine test_p2m_m2m_root_moments()
    type(fmm_plan_type) :: plan
    type(fmm_state_type) :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:), expected(:)
    integer(i32) :: alpha_idx

    call make_free_sources(src_pos, q)
    options%theta = 0.55d0
    options%leaf_max = 2_i32
    options%order = 4_i32
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    allocate (expected(plan%ncoef))
    call accumulate_direct_moments(plan, src_pos, q, plan%node_center(:, 1_i32), expected)
    do alpha_idx = 1_i32, plan%ncoef
      call assert_close_dp( &
        state%multipole(alpha_idx, 1_i32), expected(alpha_idx), &
        max(1.0d-12, 1.0d-11 * abs(expected(alpha_idx))), 'root multipole moment mismatch' &
      )
    end do

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_p2m_m2m_root_moments

  subroutine test_free_field_accuracy()
    type(fmm_plan_type) :: plan
    type(fmm_state_type) :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 8), e_fmm(3), e_ref(3)
    real(dp) :: norm_ref, rel_err, max_rel_err
    integer(i32) :: i, valid_count

    call make_free_sources(src_pos, q)
    options%theta = 0.55d0
    options%leaf_max = 4_i32
    options%order = 4_i32
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    queries(:, 1) = [1.8d0, 1.2d0, 1.5d0]
    queries(:, 2) = [-1.7d0, 1.3d0, 1.1d0]
    queries(:, 3) = [1.6d0, -1.4d0, 1.4d0]
    queries(:, 4) = [-1.5d0, -1.6d0, 1.2d0]
    queries(:, 5) = [1.9d0, 0.8d0, -1.7d0]
    queries(:, 6) = [-1.8d0, 1.0d0, -1.5d0]
    queries(:, 7) = [1.4d0, -1.9d0, -1.3d0]
    queries(:, 8) = [-1.6d0, -1.2d0, -1.8d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      call eval_point(plan, state, queries(:, i), e_fmm)
      call direct_field_free(src_pos, q, queries(:, i), e_ref)
      norm_ref = sqrt(sum(e_ref * e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_ref) * (e_fmm - e_ref))) / norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 8_i32, 'free core test lost valid samples')
    call assert_true(max_rel_err <= 3.0d-3, 'free core FMM relative error exceeds 3e-3')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_free_field_accuracy

  subroutine test_softened_free_field_accuracy()
    type(fmm_plan_type) :: plan
    type(fmm_state_type) :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 4), e_fmm(3), e_ref(3)
    real(dp) :: norm_ref, rel_err, max_rel_err
    integer(i32) :: i, valid_count

    call make_free_sources(src_pos, q)
    options%theta = 0.55d0
    options%leaf_max = 4_i32
    options%order = 4_i32
    options%softening = 5.0d-2
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    queries(:, 1) = [1.8d0, 1.2d0, 1.5d0]
    queries(:, 2) = [-1.7d0, 1.3d0, 1.1d0]
    queries(:, 3) = [1.9d0, 0.8d0, -1.7d0]
    queries(:, 4) = [-1.6d0, -1.2d0, -1.8d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      call eval_point(plan, state, queries(:, i), e_fmm)
      call direct_field_free(src_pos, q, queries(:, i), e_ref, options%softening)
      norm_ref = sqrt(sum(e_ref * e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_ref) * (e_fmm - e_ref))) / norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 4_i32, 'softened free core test lost valid samples')
    call assert_true(max_rel_err <= 3.0d-3, 'softened free core FMM relative error exceeds 3e-3')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_softened_free_field_accuracy

  subroutine test_periodic2_field_accuracy()
    type(fmm_plan_type) :: plan
    type(fmm_state_type) :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 6), e_fmm(3), e_ref(3)
    real(dp) :: norm_ref, rel_err, max_rel_err
    integer(i32) :: i, valid_count

    call make_periodic_sources(src_pos, q)
    options%theta = 0.55d0
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%use_periodic2 = .true.
    options%periodic_axes = [1_i32, 2_i32]
    options%periodic_len = [1.0d0, 1.0d0]
    options%periodic_image_layers = 1_i32
    options%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    queries(:, 1) = [0.15d0, 0.15d0, -0.60d0]
    queries(:, 2) = [0.85d0, 0.20d0, -0.20d0]
    queries(:, 3) = [0.20d0, 0.80d0, 0.10d0]
    queries(:, 4) = [0.75d0, 0.75d0, 0.50d0]
    queries(:, 5) = [0.55d0, 0.35d0, -0.75d0]
    queries(:, 6) = [0.35d0, 0.60d0, 0.85d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      call eval_point(plan, state, queries(:, i), e_fmm)
      call direct_field_periodic2(src_pos, q, queries(:, i), options%target_box_min, options%target_box_max, &
                                  options%periodic_axes, options%periodic_image_layers, e_ref)
      norm_ref = sqrt(sum(e_ref * e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_ref) * (e_fmm - e_ref))) / norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 6_i32, 'periodic2 core test lost valid samples')
    call assert_true(max_rel_err <= 5.0d-3, 'periodic2 core FMM relative error exceeds 5e-3')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_periodic2_field_accuracy

  subroutine test_periodic2_m2l_root_trunc_correction_effect()
    type(fmm_plan_type) :: plan_base, plan_root
    type(fmm_state_type) :: state_base, state_root
    type(fmm_options_type) :: options_base, options_root
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 6), e_ref(3), e_base(3), e_root(3), d_er(3)
    real(dp) :: norm_ref, rel_base, rel_root
    real(dp) :: mean_rel_base, mean_rel_root, max_delta_base_root
    integer(i32) :: i, valid_count

    call make_periodic_sources(src_pos, q)
    options_base%theta = 0.55d0
    options_base%leaf_max = 2_i32
    options_base%order = 4_i32
    options_base%use_periodic2 = .true.
    options_base%periodic_axes = [1_i32, 2_i32]
    options_base%periodic_len = [1.0d0, 1.0d0]
    options_base%periodic_image_layers = 1_i32
    options_base%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options_base%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    options_root = options_base
    options_root%periodic_far_correction = 'm2l_root_trunc'
    options_root%periodic_ewald_layers = 4_i32
    call build_plan(plan_base, src_pos, options_base)
    call build_plan(plan_root, src_pos, options_root)
    call update_state(plan_base, state_base, q)
    call update_state(plan_root, state_root, q)

    queries(:, 1) = [0.15d0, 0.15d0, -0.60d0]
    queries(:, 2) = [0.85d0, 0.20d0, -0.20d0]
    queries(:, 3) = [0.20d0, 0.80d0, 0.10d0]
    queries(:, 4) = [0.75d0, 0.75d0, 0.50d0]
    queries(:, 5) = [0.55d0, 0.35d0, -0.75d0]
    queries(:, 6) = [0.35d0, 0.60d0, 0.85d0]

    mean_rel_base = 0.0d0
    mean_rel_root = 0.0d0
    max_delta_base_root = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      call direct_field_periodic2(src_pos, q, queries(:, i), options_base%target_box_min, options_base%target_box_max, &
                                  options_base%periodic_axes, 5_i32, e_ref)
      call eval_point(plan_base, state_base, queries(:, i), e_base)
      call eval_point(plan_root, state_root, queries(:, i), e_root)
      d_er = e_root - e_base
      max_delta_base_root = max(max_delta_base_root, sqrt(sum(d_er * d_er)))

      norm_ref = sqrt(sum(e_ref * e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_base = sqrt(sum((e_base - e_ref) * (e_base - e_ref))) / norm_ref
      rel_root = sqrt(sum((e_root - e_ref) * (e_root - e_ref))) / norm_ref
      mean_rel_base = mean_rel_base + rel_base
      mean_rel_root = mean_rel_root + rel_root
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 6_i32, 'periodic2 m2l_root_trunc core test lost valid samples')
    mean_rel_base = mean_rel_base / real(valid_count, dp)
    mean_rel_root = mean_rel_root / real(valid_count, dp)
    call assert_true(max_delta_base_root > 1.0d-18, 'm2l_root_trunc core correction should affect periodic2 field')
    call assert_true(mean_rel_root <= 1.2d0 * mean_rel_base, 'm2l_root_trunc core correction degrades periodic2 accuracy too much')

    call destroy_state(state_base)
    call destroy_state(state_root)
    call destroy_plan(plan_base)
    call destroy_plan(plan_root)
  end subroutine test_periodic2_m2l_root_trunc_correction_effect

  subroutine test_periodic2_m2l_root_trunc_charged_wall_closure()
    type(fmm_plan_type) :: plan
    type(fmm_state_type) :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 5), e_ref(3), e_fmm(3)
    real(dp) :: norm_ref, rel_err, max_rel_err
    integer(i32) :: i, valid_count

    call make_periodic_sources_nonneutral(src_pos, q)
    options%theta = 0.55d0
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%use_periodic2 = .true.
    options%periodic_axes = [1_i32, 2_i32]
    options%periodic_len = [1.0d0, 1.0d0]
    options%periodic_image_layers = 1_i32
    options%periodic_far_correction = 'm2l_root_trunc'
    options%periodic_ewald_layers = 4_i32
    options%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    queries(:, 1) = [0.15d0, 0.15d0, -0.60d0]
    queries(:, 2) = [0.85d0, 0.20d0, -0.20d0]
    queries(:, 3) = [0.20d0, 0.80d0, 0.10d0]
    queries(:, 4) = [0.75d0, 0.75d0, 0.50d0]
    queries(:, 5) = [0.55d0, 0.35d0, -0.75d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      call direct_field_periodic2(src_pos, q, queries(:, i), options%target_box_min, options%target_box_max, &
                                  options%periodic_axes, 5_i32, e_ref)
      call eval_point(plan, state, queries(:, i), e_fmm)
      norm_ref = sqrt(sum(e_ref * e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_ref) * (e_fmm - e_ref))) / norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 5_i32, 'periodic2 m2l_root_trunc charged-wall test lost valid samples')
    call assert_true(max_rel_err <= 8.0d-2, 'm2l_root_trunc charged-wall closure accuracy exceeds 8e-2')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_periodic2_m2l_root_trunc_charged_wall_closure

  subroutine test_periodic2_m2l_root_oracle_correction_effect()
    type(fmm_plan_type) :: plan_base, plan_oracle
    type(fmm_state_type) :: state_base, state_oracle
    type(fmm_options_type) :: options_base, options_oracle
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 6), e_ref(3), e_base(3), e_oracle(3), d_er(3)
    real(dp) :: norm_ref, rel_base, rel_oracle
    real(dp) :: mean_rel_base, mean_rel_oracle, max_delta_base_oracle
    integer(i32) :: i, valid_count

    call make_periodic_sources(src_pos, q)
    options_base%theta = 0.55d0
    options_base%leaf_max = 2_i32
    options_base%order = 4_i32
    options_base%use_periodic2 = .true.
    options_base%periodic_axes = [1_i32, 2_i32]
    options_base%periodic_len = [1.0d0, 1.0d0]
    options_base%periodic_image_layers = 1_i32
    options_base%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options_base%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    options_oracle = options_base
    options_oracle%periodic_far_correction = 'm2l_root_oracle'
    options_oracle%periodic_ewald_layers = 4_i32
    call build_plan(plan_base, src_pos, options_base)
    call build_plan(plan_oracle, src_pos, options_oracle)
    call update_state(plan_base, state_base, q)
    call update_state(plan_oracle, state_oracle, q)
    state_oracle%profile_enabled = .true.

    queries(:, 1) = [0.15d0, 0.15d0, -0.60d0]
    queries(:, 2) = [0.85d0, 0.20d0, -0.20d0]
    queries(:, 3) = [0.20d0, 0.80d0, 0.10d0]
    queries(:, 4) = [0.75d0, 0.75d0, 0.50d0]
    queries(:, 5) = [0.55d0, 0.35d0, -0.75d0]
    queries(:, 6) = [0.35d0, 0.60d0, 1.20d0]

    call assert_true(plan_oracle%periodic_ewald%ready, 'm2l_root_oracle should precompute periodic Ewald data')
    call assert_true(plan_oracle%periodic_root_operator_ready, 'm2l_root_oracle should build a root operator')

    mean_rel_base = 0.0d0
    mean_rel_oracle = 0.0d0
    max_delta_base_oracle = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      call direct_field_periodic2(src_pos, q, queries(:, i), options_base%target_box_min, options_base%target_box_max, &
                                  options_base%periodic_axes, options_base%periodic_image_layers, e_ref)
      call add_periodic2_exact_ewald_correction_all_sources(plan_oracle, state_oracle, queries(:, i), e_ref)
      call eval_point(plan_base, state_base, queries(:, i), e_base)
      call eval_point(plan_oracle, state_oracle, queries(:, i), e_oracle)
      d_er = e_oracle - e_base
      max_delta_base_oracle = max(max_delta_base_oracle, sqrt(sum(d_er * d_er)))

      norm_ref = sqrt(sum(e_ref * e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_base = sqrt(sum((e_base - e_ref) * (e_base - e_ref))) / norm_ref
      rel_oracle = sqrt(sum((e_oracle - e_ref) * (e_oracle - e_ref))) / norm_ref
      mean_rel_base = mean_rel_base + rel_base
      mean_rel_oracle = mean_rel_oracle + rel_oracle
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 6_i32, 'periodic2 m2l_root_oracle core test lost valid samples')
    mean_rel_base = mean_rel_base / real(valid_count, dp)
    mean_rel_oracle = mean_rel_oracle / real(valid_count, dp)
    call assert_true(max_delta_base_oracle > 1.0d-18, 'm2l_root_oracle should affect periodic2 field')
    call assert_true(mean_rel_oracle <= 1.2d0 * mean_rel_base, 'm2l_root_oracle degrades periodic2 accuracy too much')
    call assert_true(state_oracle%eval_ewald_count >= 1_i32, 'm2l_root_oracle fallback should record Ewald oracle usage')
    call assert_true(state_oracle%eval_fallback_count >= 1_i32, 'm2l_root_oracle test should exercise fallback evaluation')

    call destroy_state(state_base)
    call destroy_state(state_oracle)
    call destroy_plan(plan_base)
    call destroy_plan(plan_oracle)
  end subroutine test_periodic2_m2l_root_oracle_correction_effect

  subroutine test_target_box_dual_tree()
    type(fmm_plan_type) :: plan
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)

    call make_periodic_sources(src_pos, q)
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    call build_plan(plan, src_pos, options)

    call assert_true(plan%target_tree_ready, 'core free/use_box should enable dual target tree')
    call assert_allclose_1d(plan%target_node_center(:, 1_i32), [0.5d0, 0.5d0, 0.0d0], 1.0d-12, 'target root center mismatch')
    call assert_allclose_1d(plan%target_node_half_size(:, 1_i32), [0.5d0, 0.5d0, 1.0d0], 1.0d-12, 'target root half-size mismatch')

    call destroy_plan(plan)
  end subroutine test_target_box_dual_tree

  subroutine test_state_update_reuse()
    type(fmm_plan_type) :: plan
    type(fmm_state_type) :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q1(:), q2(:)
    real(dp) :: e1(3), e2(3)
    integer(i32) :: pair_count, build_count

    call make_periodic_sources(src_pos, q1)
    q2 = 2.0d0 * q1
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%use_periodic2 = .true.
    options%periodic_axes = [1_i32, 2_i32]
    options%periodic_len = [1.0d0, 1.0d0]
    options%periodic_image_layers = 1_i32
    options%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    call build_plan(plan, src_pos, options)
    pair_count = plan%m2l_pair_count
    build_count = plan%m2l_build_count

    call update_state(plan, state, q1)
    call eval_point(plan, state, [0.35d0, 0.45d0, 0.55d0], e1)
    call update_state(plan, state, q2)
    call eval_point(plan, state, [0.35d0, 0.45d0, 0.55d0], e2)

    call assert_equal_i32(plan%m2l_pair_count, pair_count, 'state update should preserve M2L pair count')
    call assert_equal_i32(plan%m2l_build_count, build_count, 'state update should not rebuild the M2L cache')
    call assert_equal_i32(state%update_count, 2_i32, 'state update count mismatch')
    call assert_true(sqrt(sum((e2 - 2.0d0 * e1) * (e2 - 2.0d0 * e1))) <= 1.0d-10, 'field should scale linearly with charge')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_state_update_reuse

  subroutine test_state_eval_profile_counts()
    type(fmm_plan_type) :: plan
    type(fmm_state_type) :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: e(3)

    call make_periodic_sources(src_pos, q)
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%use_periodic2 = .true.
    options%periodic_axes = [1_i32, 2_i32]
    options%periodic_len = [1.0d0, 1.0d0]
    options%periodic_image_layers = 1_i32
    options%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)
    state%profile_enabled = .true.

    call eval_point(plan, state, [0.35d0, 0.45d0, 0.55d0], e)
    call eval_point(plan, state, [0.40d0, 0.50d0, -0.25d0], e)

    call assert_equal_i32(state%eval_count, 2_i32, 'profile eval_count mismatch')
    call assert_equal_i32(state%eval_fallback_count, 0_i32, 'profile should stay on target tree for this fixture')
    call assert_true(state%eval_near_source_count > 0, 'profile should count near sources')
    call assert_equal_i32( &
      int(state%eval_direct_kernel_count, i32), int(state%eval_near_source_count, i32), &
      'periodic2 direct-kernel count should match retained near source-image entries' &
    )

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_state_eval_profile_counts

  subroutine test_field_solver_core_adapter()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: r(3), e_direct(3), e_fmm(3), rel_err, norm_ref, max_rel_err
    integer(i32) :: i, valid_count

    call make_sphere(mesh_fmm, radius=0.35d0, n_lon=12_i32, n_lat=6_i32, center=[0.0d0, 0.0d0, 0.0d0])
    mesh_fmm%q_elem = 1.0d-12

    sim = sim_config()
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'free'
    sim%softening = 0.0d0
    call solver%init(mesh_fmm, sim)
    call assert_true(solver%fmm_use_core, 'softening=0 free FMM should use the core path')
    call solver%refresh(mesh_fmm)

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, 4_i32
      r = [1.2d0 + 0.1d0 * real(i - 1_i32, dp), -0.8d0 + 0.2d0 * real(i - 1_i32, dp), 0.9d0 - 0.1d0 * real(i - 1_i32, dp)]
      call electric_field_at(mesh_fmm, r, 0.0d0, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)
      norm_ref = sqrt(sum(e_direct * e_direct))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_direct) * (e_fmm - e_direct))) / norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 4_i32, 'core adapter test lost valid samples')
    call assert_true(max_rel_err <= 5.0d-3, 'core adapter relative error exceeds 5e-3')
  end subroutine test_field_solver_core_adapter

  subroutine test_field_solver_core_softened_adapter()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: r(3), e_direct(3), e_fmm(3), rel_err, norm_ref, max_rel_err
    integer(i32) :: i, valid_count

    call make_sphere(mesh_fmm, radius=0.35d0, n_lon=12_i32, n_lat=6_i32, center=[0.0d0, 0.0d0, 0.0d0])
    mesh_fmm%q_elem = 1.0d-12

    sim = sim_config()
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'free'
    sim%softening = 1.0d-4
    call solver%init(mesh_fmm, sim)
    call assert_true(solver%fmm_use_core, 'softening>0 free FMM should use the core path')
    call solver%refresh(mesh_fmm)

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, 4_i32
      r = [1.2d0 + 0.1d0 * real(i - 1_i32, dp), -0.8d0 + 0.2d0 * real(i - 1_i32, dp), 0.9d0 - 0.1d0 * real(i - 1_i32, dp)]
      call electric_field_at(mesh_fmm, r, sim%softening, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)
      norm_ref = sqrt(sum(e_direct * e_direct))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_direct) * (e_fmm - e_direct))) / norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 4_i32, 'softened core adapter test lost valid samples')
    call assert_true(max_rel_err <= 5.0d-3, 'softened core adapter relative error exceeds 5e-3')
  end subroutine test_field_solver_core_softened_adapter

  subroutine test_field_solver_core_periodic2_m2l_root_trunc_adapter()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver_default = field_solver_type()
    type(field_solver_type) :: solver_root = field_solver_type()
    type(sim_config) :: sim
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 4), e_raw(3), e_ref(3), e_default(3), e_root(3)
    real(dp) :: norm_ref, rel_default, rel_root, mean_rel_default, mean_rel_root, max_delta_default_root
    integer(i32) :: i, valid_count, ref_layers

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=8_i32, n_lat=4_i32, center=[0.5d0, 0.5d0, 0.0d0])
    do i = 1_i32, mesh_fmm%nelem
      if (mod(i, 2_i32) == 0_i32) then
        mesh_fmm%q_elem(i) = -1.0d-12
      else
        mesh_fmm%q_elem(i) = 1.0d-12
      end if
    end do

    sim = sim_config()
    sim%softening = 0.0d0
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_image_layers = 1_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver_default%init(mesh_fmm, sim)
    call solver_default%refresh(mesh_fmm)

    sim%field_periodic_far_correction = 'm2l_root_trunc'
    sim%field_periodic_ewald_layers = 4_i32
    call solver_root%init(mesh_fmm, sim)
    call solver_root%refresh(mesh_fmm)

    call assert_true(solver_default%fmm_use_core, 'softening=0 periodic2 default FMM should use the core path')
    call assert_true( &
      solver_root%fmm_use_core, &
      'softening=0 periodic2 m2l_root_trunc FMM should use the core path' &
    )
    call assert_true( &
      trim(solver_default%periodic_far_correction) == 'm2l_root_trunc', &
      'periodic2 default should normalize to m2l_root_trunc' &
    )
    call assert_true( &
      trim(solver_default%fmm_core_options%periodic_far_correction) == 'm2l_root_trunc', &
      'core adapter should pass normalized m2l_root_trunc into FMM options' &
    )

    call mesh_centers_as_sources(mesh_fmm, src_pos, q)
    queries(:, 1) = [0.15d0, 0.15d0, -0.60d0]
    queries(:, 2) = [0.75d0, 0.75d0, 0.50d0]

    ref_layers = solver_default%periodic_image_layers + solver_default%periodic_ewald_layers
    mean_rel_default = 0.0d0
    mean_rel_root = 0.0d0
    max_delta_default_root = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, 2_i32
      call direct_field_periodic2(src_pos, q, queries(:, i), sim%box_min, sim%box_max, [1_i32, 2_i32], ref_layers, e_raw)
      e_ref = k_coulomb * e_raw
      call solver_default%eval_e(mesh_fmm, queries(:, i), e_default)
      call solver_root%eval_e(mesh_fmm, queries(:, i), e_root)
      max_delta_default_root = max(max_delta_default_root, sqrt(sum((e_default - e_root) * (e_default - e_root))))

      norm_ref = sqrt(sum(e_ref * e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_default = sqrt(sum((e_default - e_ref) * (e_default - e_ref))) / norm_ref
      rel_root = sqrt(sum((e_root - e_ref) * (e_root - e_ref))) / norm_ref
      mean_rel_default = mean_rel_default + rel_default
      mean_rel_root = mean_rel_root + rel_root
      valid_count = valid_count + 1_i32
    end do

    call assert_true( &
      valid_count == 2_i32, &
      'core periodic2 default m2l_root_trunc adapter test lost valid samples' &
    )
    mean_rel_default = mean_rel_default / real(valid_count, dp)
    mean_rel_root = mean_rel_root / real(valid_count, dp)
    call assert_true( &
      max_delta_default_root <= 1.0d-18 * max(1.0d0, sqrt(sum(e_ref * e_ref))), &
      'default periodic2 and explicit m2l_root_trunc should agree at the adapter level' &
    )
    call assert_true( &
      mean_rel_default <= 8.0d-2 .and. mean_rel_root <= 8.0d-2, &
      'core adapter periodic2 m2l_root_trunc accuracy exceeds 8e-2' &
    )
  end subroutine test_field_solver_core_periodic2_m2l_root_trunc_adapter

  subroutine make_free_sources(src_pos, q)
    real(dp), allocatable, intent(out) :: src_pos(:, :)
    real(dp), allocatable, intent(out) :: q(:)
    real(dp), parameter :: axis_vals(3) = [-0.45d0, 0.0d0, 0.45d0]
    integer(i32) :: ix, iy, iz, idx

    allocate (src_pos(3, 27), q(27))
    idx = 0_i32
    do ix = 1_i32, 3_i32
      do iy = 1_i32, 3_i32
        do iz = 1_i32, 3_i32
          idx = idx + 1_i32
          src_pos(:, idx) = [axis_vals(ix), axis_vals(iy), axis_vals(iz)]
          q(idx) = real((-1_i32)**idx, dp) * 1.0d-12 * (1.0d0 + 0.05d0 * real(mod(idx, 5_i32), dp))
        end do
      end do
    end do
  end subroutine make_free_sources

  subroutine make_periodic_sources(src_pos, q)
    real(dp), allocatable, intent(out) :: src_pos(:, :)
    real(dp), allocatable, intent(out) :: q(:)

    allocate (src_pos(3, 8), q(8))
    src_pos(:, 1) = [0.20d0, 0.20d0, -0.40d0]
    src_pos(:, 2) = [0.80d0, 0.20d0, -0.20d0]
    src_pos(:, 3) = [0.25d0, 0.75d0, 0.10d0]
    src_pos(:, 4) = [0.75d0, 0.70d0, 0.35d0]
    src_pos(:, 5) = [0.40d0, 0.35d0, -0.70d0]
    src_pos(:, 6) = [0.60d0, 0.55d0, 0.55d0]
    src_pos(:, 7) = [0.30d0, 0.60d0, 0.85d0]
    src_pos(:, 8) = [0.65d0, 0.30d0, -0.85d0]
    q = [ &
      -1.0d-12, 1.1d-12, -1.2d-12, 1.3d-12, &
      -1.4d-12, 1.5d-12, -1.6d-12, 1.3d-12 &
    ]
  end subroutine make_periodic_sources

  subroutine make_periodic_sources_nonneutral(src_pos, q)
    real(dp), allocatable, intent(out) :: src_pos(:, :)
    real(dp), allocatable, intent(out) :: q(:)

    allocate (src_pos(3, 4), q(4))
    src_pos(:, 1) = [0.20d0, 0.20d0, -0.30d0]
    src_pos(:, 2) = [0.78d0, 0.30d0, 0.15d0]
    src_pos(:, 3) = [0.32d0, 0.75d0, 0.45d0]
    src_pos(:, 4) = [0.68d0, 0.68d0, -0.55d0]
    q(1) = 1.0d-12
    q(2) = 0.4d-12
    q(3) = -0.2d-12
    q(4) = 0.7d-12
  end subroutine make_periodic_sources_nonneutral

  subroutine accumulate_direct_moments(plan, src_pos, q, center, moments)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: src_pos(:, :), q(:), center(3)
    real(dp), intent(out) :: moments(:)
    integer(i32) :: idx, alpha_idx
    real(dp) :: d(3)
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    moments = 0.0d0
    do idx = 1_i32, int(size(q), i32)
      d = src_pos(:, idx) - center
      call build_powers(d, plan%options%order, xpow, ypow, zpow)
      do alpha_idx = 1_i32, plan%ncoef
        moments(alpha_idx) = moments(alpha_idx) + q(idx) * xpow(plan%alpha(1, alpha_idx)) * ypow(plan%alpha(2, alpha_idx)) &
                             * zpow(plan%alpha(3, alpha_idx)) / plan%alpha_factorial(alpha_idx)
      end do
    end do
  end subroutine accumulate_direct_moments

  subroutine direct_field_free(src_pos, q, r, e, softening)
    real(dp), intent(in) :: src_pos(:, :), q(:), r(3)
    real(dp), intent(out) :: e(3)
    real(dp), intent(in), optional :: softening
    integer(i32) :: idx
    real(dp) :: d(3), r2, inv_r3, soft2

    e = 0.0d0
    soft2 = 0.0d0
    if (present(softening)) soft2 = softening * softening
    do idx = 1_i32, int(size(q), i32)
      d = r - src_pos(:, idx)
      r2 = sum(d * d) + soft2
      if (r2 <= tiny(1.0d0)) cycle
      inv_r3 = 1.0d0 / (sqrt(r2) * r2)
      e = e + q(idx) * inv_r3 * d
    end do
  end subroutine direct_field_free

  subroutine mesh_centers_as_sources(mesh, src_pos, q)
    type(mesh_type), intent(in) :: mesh
    real(dp), allocatable, intent(out) :: src_pos(:, :)
    real(dp), allocatable, intent(out) :: q(:)
    integer(i32) :: idx

    allocate (src_pos(3, mesh%nelem), q(mesh%nelem))
    do idx = 1_i32, mesh%nelem
      src_pos(:, idx) = [mesh%center_x(idx), mesh%center_y(idx), mesh%center_z(idx)]
      q(idx) = mesh%q_elem(idx)
    end do
  end subroutine mesh_centers_as_sources

  subroutine direct_field_periodic2(src_pos, q, r, box_min, box_max, periodic_axes, nimg, e, softening)
    real(dp), intent(in) :: src_pos(:, :), q(:), r(3), box_min(3), box_max(3)
    integer(i32), intent(in) :: periodic_axes(2), nimg
    real(dp), intent(out) :: e(3)
    real(dp), intent(in), optional :: softening
    integer(i32) :: idx, img_i, img_j, axis1, axis2
    real(dp) :: shifted(3), d(3), r2, inv_r3, l1, l2, soft2

    axis1 = periodic_axes(1)
    axis2 = periodic_axes(2)
    l1 = box_max(axis1) - box_min(axis1)
    l2 = box_max(axis2) - box_min(axis2)
    e = 0.0d0
    soft2 = 0.0d0
    if (present(softening)) soft2 = softening * softening
    do idx = 1_i32, int(size(q), i32)
      do img_i = -nimg, nimg
        do img_j = -nimg, nimg
          shifted = src_pos(:, idx)
          shifted(axis1) = shifted(axis1) + real(img_i, dp) * l1
          shifted(axis2) = shifted(axis2) + real(img_j, dp) * l2
          d = r - shifted
          r2 = sum(d * d) + soft2
          if (r2 <= tiny(1.0d0)) cycle
          inv_r3 = 1.0d0 / (sqrt(r2) * r2)
          e = e + q(idx) * inv_r3 * d
        end do
      end do
    end do
  end subroutine direct_field_periodic2

  subroutine build_powers(d, order, xpow, ypow, zpow)
    real(dp), intent(in) :: d(3)
    integer(i32), intent(in) :: order
    real(dp), intent(out) :: xpow(0:order), ypow(0:order), zpow(0:order)
    integer(i32) :: k

    xpow(0) = 1.0d0
    ypow(0) = 1.0d0
    zpow(0) = 1.0d0
    do k = 1_i32, order
      xpow(k) = xpow(k - 1_i32) * d(1)
      ypow(k) = ypow(k - 1_i32) * d(2)
      zpow(k) = zpow(k - 1_i32) * d(3)
    end do
  end subroutine build_powers

end program test_coulomb_fmm_core
