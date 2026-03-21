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
  use bem_coulomb_fmm_periodic_ewald, only: add_periodic2_exact_ewald_correction_all_sources, &
                                            add_periodic2_exact_ewald_correction_single_source
  use test_support, only: assert_true, assert_close_dp, assert_allclose_1d, assert_equal_i32
  implicit none

  call test_p2m_m2m_root_moments()
  call test_free_field_accuracy()
  call test_softened_free_field_accuracy()
  call test_periodic2_field_accuracy()
  call test_periodic2_target_oracle_operator_residual_accuracy()
  call test_periodic2_m2l_root_oracle_correction_effect()
  call test_periodic2_nonneutral_charged_wall_outside_box()
  call test_target_box_dual_tree()
  call test_state_update_reuse()
  call test_field_solver_core_adapter()
  call test_field_solver_core_softened_adapter()
  call test_field_solver_core_periodic2_m2l_root_oracle_adapter()

contains

  subroutine test_p2m_m2m_root_moments()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:), expected(:)
    integer(i32) :: alpha_idx

    allocate (plan, state)
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
        max(1.0d-12, 1.0d-11*abs(expected(alpha_idx))), 'root multipole moment mismatch' &
        )
    end do

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_p2m_m2m_root_moments

  subroutine test_free_field_accuracy()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 8), e_fmm(3), e_ref(3)
    real(dp) :: norm_ref, rel_err, max_rel_err
    integer(i32) :: i, valid_count

    allocate (plan, state)
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
      norm_ref = sqrt(sum(e_ref*e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_ref)*(e_fmm - e_ref)))/norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 8_i32, 'free core test lost valid samples')
    call assert_true(max_rel_err <= 3.0d-3, 'free core FMM relative error exceeds 3e-3')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_free_field_accuracy

  subroutine test_softened_free_field_accuracy()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 4), e_fmm(3), e_ref(3)
    real(dp) :: norm_ref, rel_err, max_rel_err
    integer(i32) :: i, valid_count

    allocate (plan, state)
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
      norm_ref = sqrt(sum(e_ref*e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_ref)*(e_fmm - e_ref)))/norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 4_i32, 'softened free core test lost valid samples')
    call assert_true(max_rel_err <= 3.0d-3, 'softened free core FMM relative error exceeds 3e-3')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_softened_free_field_accuracy

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
      if (trim(plan%options%periodic_far_correction) == 'm2l_root_oracle') then
        call add_periodic2_exact_ewald_correction_all_sources(plan, state, queries(:, i), e_ref)
      end if
      norm_ref = sqrt(sum(e_ref*e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_ref)*(e_fmm - e_ref)))/norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 6_i32, 'periodic2 core test lost valid samples')
    call assert_true(max_rel_err <= 5.0d-2, 'periodic2 core FMM relative error exceeds 5e-2')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_periodic2_field_accuracy

  subroutine test_periodic2_target_oracle_operator_residual_accuracy()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:), probes(:, :), target_local(:)
    real(dp) :: e_operator(3), e_ref(3)
    real(dp) :: norm_ref, rel_err, max_rel_err, mean_rel_err, global_rel_err
    real(dp) :: sq_err_sum, sq_ref_sum
    integer(i32) :: i, nprobe, valid_count, ntarget_sample, sample_idx, target_idx, node_idx

    allocate (plan, state)
    call make_periodic_sources(src_pos, q)
    options%theta = 0.55d0
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%use_periodic2 = .true.
    options%periodic_axes = [1_i32, 2_i32]
    options%periodic_len = [1.0d0, 1.0d0]
    options%periodic_image_layers = 1_i32
    options%periodic_far_correction = 'm2l_root_oracle'
    options%periodic_ewald_layers = 4_i32
    options%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    call assert_true(plan%periodic_root_operator_ready, 'm2l_root_oracle should build target operators')
    call assert_true(plan%periodic_root_target_count > 1_i32, 'target oracle should cover multiple target nodes')
    call assert_true( &
      size(plan%periodic_root_operator, 3) == plan%periodic_root_target_count, &
      'target operator count should match stored target nodes' &
      )
    nprobe = 4_i32*plan%ncoef
    ntarget_sample = min(plan%periodic_root_target_count, 8_i32)
    allocate (probes(3, nprobe), target_local(plan%ncoef))

    max_rel_err = 0.0d0
    mean_rel_err = 0.0d0
    sq_err_sum = 0.0d0
    sq_ref_sum = 0.0d0
    valid_count = 0_i32
    do sample_idx = 1_i32, ntarget_sample
      if (ntarget_sample == 1_i32) then
        target_idx = 1_i32
      else
        target_idx = 1_i32 + int( &
                     real(plan%periodic_root_target_count - 1_i32, dp)*real(sample_idx - 1_i32, dp)/ &
                     real(ntarget_sample - 1_i32, dp), i32 &
                     )
      end if
      node_idx = plan%periodic_root_target_nodes(target_idx)
      call build_test_root_sample_points( &
        plan%target_node_center(:, node_idx), plan%target_node_half_size(:, node_idx), nprobe, 0.61d0, probes &
        )
      target_local = matmul(plan%periodic_root_operator(:, :, target_idx), state%multipole(:, 1_i32))
      do i = 1_i32, nprobe
        call eval_local_field_from_coeff(plan, plan%target_node_center(:, node_idx), target_local, probes(:, i), e_operator)
        e_ref = 0.0d0
        call add_periodic2_exact_ewald_correction_all_sources(plan, state, probes(:, i), e_ref)
        norm_ref = sqrt(sum(e_ref*e_ref))
        if (norm_ref <= 1.0d-18) cycle
        rel_err = sqrt(sum((e_operator - e_ref)*(e_operator - e_ref)))/norm_ref
        max_rel_err = max(max_rel_err, rel_err)
        mean_rel_err = mean_rel_err + rel_err
        sq_err_sum = sq_err_sum + sum((e_operator - e_ref)*(e_operator - e_ref))
        sq_ref_sum = sq_ref_sum + sum(e_ref*e_ref)
        valid_count = valid_count + 1_i32
      end do
    end do

    call assert_true( &
      valid_count >= max(1_i32, ntarget_sample*nprobe/2_i32), &
      'periodic2 target oracle operator test lost too many valid samples' &
      )
    mean_rel_err = mean_rel_err/real(valid_count, dp)
    global_rel_err = sqrt(sq_err_sum/max(sq_ref_sum, tiny(1.0d0)))
    write (*, '(A,I0,A,ES12.5,A,ES12.5,A,ES12.5)') &
      'test_periodic2_target_oracle_operator_residual_accuracy: valid_count=', valid_count, &
      ', global_rel_err=', global_rel_err, ', mean_rel_err=', mean_rel_err, ', max_rel_err=', max_rel_err
    call assert_true(global_rel_err <= 5.0d-2, 'periodic2 target oracle operator global relative error exceeds 5e-2')
    call assert_true(max_rel_err <= 1.0d-1, 'periodic2 target oracle operator pointwise relative error exceeds 1e-1')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_periodic2_target_oracle_operator_residual_accuracy

  subroutine test_periodic2_m2l_root_oracle_correction_effect()
    type(fmm_plan_type), allocatable :: plan_oracle
    type(fmm_state_type), allocatable :: state_oracle
    type(fmm_options_type) :: options_oracle
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 6), e_ref(3), e_oracle(3)
    real(dp) :: norm_ref, rel_oracle, max_rel_oracle
    integer(i32) :: i, valid_count

    allocate (plan_oracle, state_oracle)
    call make_periodic_sources(src_pos, q)
    options_oracle%theta = 0.55d0
    options_oracle%leaf_max = 2_i32
    options_oracle%order = 4_i32
    options_oracle%use_periodic2 = .true.
    options_oracle%periodic_axes = [1_i32, 2_i32]
    options_oracle%periodic_len = [1.0d0, 1.0d0]
    options_oracle%periodic_image_layers = 1_i32
    options_oracle%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options_oracle%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    options_oracle%periodic_far_correction = 'm2l_root_oracle'
    options_oracle%periodic_ewald_layers = 4_i32
    call build_plan(plan_oracle, src_pos, options_oracle)
    call update_state(plan_oracle, state_oracle, q)

    queries(:, 1) = [0.15d0, 0.15d0, -0.60d0]
    queries(:, 2) = [0.85d0, 0.20d0, -0.20d0]
    queries(:, 3) = [0.20d0, 0.80d0, 0.10d0]
    queries(:, 4) = [0.75d0, 0.75d0, 0.50d0]
    queries(:, 5) = [0.55d0, 0.35d0, -0.75d0]
    queries(:, 6) = [0.35d0, 0.60d0, 1.20d0]

    call assert_true(plan_oracle%periodic_ewald%ready, 'm2l_root_oracle should precompute periodic Ewald data')
    call assert_true(plan_oracle%periodic_root_operator_ready, 'm2l_root_oracle should build target operators')

    max_rel_oracle = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      call direct_field_periodic2(src_pos, q, queries(:, i), options_oracle%target_box_min, options_oracle%target_box_max, &
                                  options_oracle%periodic_axes, options_oracle%periodic_image_layers, e_ref)
      call add_periodic2_exact_ewald_correction_all_sources(plan_oracle, state_oracle, queries(:, i), e_ref)
      call eval_point(plan_oracle, state_oracle, queries(:, i), e_oracle)

      norm_ref = sqrt(sum(e_ref*e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_oracle = sqrt(sum((e_oracle - e_ref)*(e_oracle - e_ref)))/norm_ref
      max_rel_oracle = max(max_rel_oracle, rel_oracle)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 6_i32, 'periodic2 m2l_root_oracle core test lost valid samples')
    call assert_true(max_rel_oracle <= 8.0d-2, 'm2l_root_oracle periodic2 accuracy exceeds 8e-2')

    call destroy_state(state_oracle)
    call destroy_plan(plan_oracle)
  end subroutine test_periodic2_m2l_root_oracle_correction_effect

  subroutine test_periodic2_nonneutral_charged_wall_outside_box()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 2), e_fmm(3), e_ref(3)
    real(dp) :: norm_ref, rel_err, max_rel_err, total_charge
    integer(i32) :: i, idx, valid_count

    allocate (plan, state)
    call make_periodic_sources_nonneutral(src_pos, q)
    total_charge = sum(q)
    options%theta = 0.55d0
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%use_periodic2 = .true.
    options%periodic_axes = [1_i32, 2_i32]
    options%periodic_len = [1.0d0, 1.0d0]
    options%periodic_image_layers = 1_i32
    options%periodic_far_correction = 'm2l_root_oracle'
    options%periodic_ewald_layers = 4_i32
    options%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    queries(:, 1) = [0.45d0, 0.40d0, 1.30d0]
    queries(:, 2) = [0.45d0, 0.40d0, -1.25d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      call eval_point(plan, state, queries(:, i), e_fmm)
      call direct_field_periodic2(src_pos, q, queries(:, i), options%target_box_min, options%target_box_max, &
                                  options%periodic_axes, options%periodic_image_layers, e_ref)
      do idx = 1_i32, int(size(q), i32)
        call add_periodic2_exact_ewald_correction_single_source(plan, q(idx), src_pos(:, idx), queries(:, i), e_ref)
      end do
      call add_expected_charged_wall_field( &
        total_charge, options%target_box_min(3), options%target_box_max(3), plan%periodic_ewald%cell_area, 3_i32, &
        queries(:, i), e_ref &
        )
      norm_ref = sqrt(sum(e_ref*e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_ref)*(e_fmm - e_ref)))/norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 2_i32, 'periodic2 nonneutral charged-wall test lost valid samples')
    call assert_true(max_rel_err <= 1.0d-12, 'periodic2 nonneutral charged-wall correction mismatch')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_periodic2_nonneutral_charged_wall_outside_box

  subroutine test_target_box_dual_tree()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)

    allocate (plan)
    call make_periodic_sources(src_pos, q)
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    call build_plan(plan, src_pos, options)

    call assert_true( &
      plan%target_tree_ready, &
      'core free/use_box should enable dual target tree' &
      )
    call assert_allclose_1d( &
      plan%target_node_center(:, 1_i32), [0.5d0, 0.5d0, 0.0d0], 1.0d-12, &
      'target root center mismatch' &
      )
    call assert_allclose_1d( &
      plan%target_node_half_size(:, 1_i32), [0.5d0, 0.5d0, 1.0d0], 1.0d-12, &
      'target root half-size mismatch' &
      )

    call destroy_plan(plan)
  end subroutine test_target_box_dual_tree

  subroutine test_state_update_reuse()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q1(:), q2(:)
    real(dp) :: e1(3), e2(3)
    integer(i32) :: pair_count, build_count

    allocate (plan, state)
    call make_periodic_sources(src_pos, q1)
    q2 = 2.0d0*q1
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
    call assert_true(sqrt(sum((e2 - 2.0d0*e1)*(e2 - 2.0d0*e1))) <= 1.0d-10, 'field should scale linearly with charge')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_state_update_reuse

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
      r = [1.2d0 + 0.1d0*real(i - 1_i32, dp), -0.8d0 + 0.2d0*real(i - 1_i32, dp), 0.9d0 - 0.1d0*real(i - 1_i32, dp)]
      call electric_field_at(mesh_fmm, r, 0.0d0, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)
      norm_ref = sqrt(sum(e_direct*e_direct))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))/norm_ref
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
      r = [1.2d0 + 0.1d0*real(i - 1_i32, dp), -0.8d0 + 0.2d0*real(i - 1_i32, dp), 0.9d0 - 0.1d0*real(i - 1_i32, dp)]
      call electric_field_at(mesh_fmm, r, sim%softening, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)
      norm_ref = sqrt(sum(e_direct*e_direct))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))/norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 4_i32, 'softened core adapter test lost valid samples')
    call assert_true(max_rel_err <= 5.0d-3, 'softened core adapter relative error exceeds 5e-3')
  end subroutine test_field_solver_core_softened_adapter

  subroutine test_field_solver_core_periodic2_m2l_root_oracle_adapter()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver_default = field_solver_type()
    type(field_solver_type) :: solver_root = field_solver_type()
    type(sim_config) :: sim
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 4), e_raw(3), e_ref(3), e_default(3), e_root(3)
    real(dp) :: norm_ref, rel_default, rel_root, mean_rel_default, mean_rel_root, max_delta_default_root
    integer(i32) :: i, valid_count

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

    sim%field_periodic_far_correction = 'm2l_root_oracle'
    sim%field_periodic_ewald_layers = 4_i32
    call solver_root%init(mesh_fmm, sim)
    call solver_root%refresh(mesh_fmm)

    call assert_true(solver_default%fmm_use_core, 'softening=0 periodic2 default FMM should use the core path')
    call assert_true( &
      solver_root%fmm_use_core, &
      'softening=0 periodic2 m2l_root_oracle FMM should use the core path' &
      )
    call assert_true( &
      trim(solver_default%periodic_far_correction) == 'm2l_root_oracle', &
      'periodic2 default should normalize to m2l_root_oracle' &
      )
    call assert_true( &
      trim(solver_default%fmm_core_options%periodic_far_correction) == 'm2l_root_oracle', &
      'core adapter should pass normalized m2l_root_oracle into FMM options' &
      )

    call mesh_centers_as_sources(mesh_fmm, src_pos, q)
    queries(:, 1) = [0.15d0, 0.15d0, -0.60d0]
    queries(:, 2) = [0.75d0, 0.75d0, 0.50d0]

    mean_rel_default = 0.0d0
    mean_rel_root = 0.0d0
    max_delta_default_root = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, 2_i32
      call direct_field_periodic2( &
        src_pos, q, queries(:, i), sim%box_min, sim%box_max, [1_i32, 2_i32], solver_default%periodic_image_layers, e_raw &
        )
      call add_periodic2_exact_ewald_correction_all_sources( &
        solver_root%fmm_core_plan, solver_root%fmm_core_state, queries(:, i), e_raw &
        )
      e_ref = k_coulomb*e_raw
      call solver_default%eval_e(mesh_fmm, queries(:, i), e_default)
      call solver_root%eval_e(mesh_fmm, queries(:, i), e_root)
      max_delta_default_root = max(max_delta_default_root, sqrt(sum((e_default - e_root)*(e_default - e_root))))

      norm_ref = sqrt(sum(e_ref*e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_default = sqrt(sum((e_default - e_ref)*(e_default - e_ref)))/norm_ref
      rel_root = sqrt(sum((e_root - e_ref)*(e_root - e_ref)))/norm_ref
      mean_rel_default = mean_rel_default + rel_default
      mean_rel_root = mean_rel_root + rel_root
      valid_count = valid_count + 1_i32
    end do

    call assert_true( &
      valid_count == 2_i32, &
      'core periodic2 default m2l_root_oracle adapter test lost valid samples' &
      )
    mean_rel_default = mean_rel_default/real(valid_count, dp)
    mean_rel_root = mean_rel_root/real(valid_count, dp)
    call assert_true( &
      max_delta_default_root <= 1.0d-18*max(1.0d0, sqrt(sum(e_ref*e_ref))), &
      'default periodic2 and explicit m2l_root_oracle should agree at the adapter level' &
      )
    call assert_true( &
      mean_rel_default <= 8.0d-2 .and. mean_rel_root <= 8.0d-2, &
      'core adapter periodic2 m2l_root_oracle accuracy exceeds 8e-2' &
      )
  end subroutine test_field_solver_core_periodic2_m2l_root_oracle_adapter

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
          q(idx) = real((-1_i32)**idx, dp)*1.0d-12*(1.0d0 + 0.05d0*real(mod(idx, 5_i32), dp))
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
        moments(alpha_idx) = moments(alpha_idx) + q(idx)*xpow(plan%alpha(1, alpha_idx))*ypow(plan%alpha(2, alpha_idx)) &
                             *zpow(plan%alpha(3, alpha_idx))/plan%alpha_factorial(alpha_idx)
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
    if (present(softening)) soft2 = softening*softening
    do idx = 1_i32, int(size(q), i32)
      d = r - src_pos(:, idx)
      r2 = sum(d*d) + soft2
      if (r2 <= tiny(1.0d0)) cycle
      inv_r3 = 1.0d0/(sqrt(r2)*r2)
      e = e + q(idx)*inv_r3*d
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
    if (present(softening)) soft2 = softening*softening
    do idx = 1_i32, int(size(q), i32)
      do img_i = -nimg, nimg
        do img_j = -nimg, nimg
          shifted = src_pos(:, idx)
          shifted(axis1) = shifted(axis1) + real(img_i, dp)*l1
          shifted(axis2) = shifted(axis2) + real(img_j, dp)*l2
          d = r - shifted
          r2 = sum(d*d) + soft2
          if (r2 <= tiny(1.0d0)) cycle
          inv_r3 = 1.0d0/(sqrt(r2)*r2)
          e = e + q(idx)*inv_r3*d
        end do
      end do
    end do
  end subroutine direct_field_periodic2

  subroutine add_expected_charged_wall_field(total_charge, z_low, z_high, cell_area, axis_free, target, e)
    real(dp), intent(in) :: total_charge, z_low, z_high, cell_area, target(3)
    integer(i32), intent(in) :: axis_free
    real(dp), intent(inout) :: e(3)
    real(dp), parameter :: two_pi_dp = 2.0d0*acos(-1.0d0)
    real(dp) :: pref, tol

    if (abs(total_charge) <= tiny(1.0d0) .or. cell_area <= 0.0d0) return
    tol = 1.0d-12*max(1.0d0, max(abs(z_low), abs(z_high)))
    pref = two_pi_dp*total_charge/cell_area
    if (target(axis_free) < z_low - tol) then
      e(axis_free) = e(axis_free) + pref
    else if (target(axis_free) > z_high + tol) then
      e(axis_free) = e(axis_free) - pref
    end if
  end subroutine add_expected_charged_wall_field

  subroutine build_test_root_sample_points(center, half_size, npoint, offset, points)
    real(dp), intent(in) :: center(3), half_size(3), offset
    integer(i32), intent(in) :: npoint
    real(dp), intent(out) :: points(3, npoint)
    integer(i32) :: idx
    real(dp) :: f1, f2, f3
    real(dp), parameter :: g1 = 0.7548776662466927d0
    real(dp), parameter :: g2 = 0.5698402909980532d0
    real(dp), parameter :: g3 = 0.4301597090019468d0

    do idx = 1_i32, npoint
      f1 = modulo(offset + real(idx, dp)*g1, 1.0d0)
      f2 = modulo(offset + real(idx, dp)*g2, 1.0d0)
      f3 = modulo(offset + real(idx, dp)*g3, 1.0d0)
      points(1, idx) = center(1) + (2.0d0*(0.05d0 + 0.9d0*f1) - 1.0d0)*half_size(1)
      points(2, idx) = center(2) + (2.0d0*(0.05d0 + 0.9d0*f2) - 1.0d0)*half_size(2)
      points(3, idx) = center(3) + (2.0d0*(0.05d0 + 0.9d0*f3) - 1.0d0)*half_size(3)
    end do
  end subroutine build_test_root_sample_points

  subroutine eval_local_field_from_coeff(plan, center, local_coeff, target, e)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: center(3), local_coeff(:), target(3)
    real(dp), intent(out) :: e(3)
    integer(i32) :: term_idx
    real(dp) :: dr(3), monomial
    real(dp) :: xpow(0:plan%options%order), ypow(0:plan%options%order), zpow(0:plan%options%order)

    if (size(local_coeff) /= plan%ncoef) error stop 'eval_local_field_from_coeff size mismatch.'

    e = 0.0d0
    if (plan%eval_term_count <= 0_i32) return

    dr = target - center
    call build_powers(dr, plan%options%order, xpow, ypow, zpow)
    do term_idx = 1_i32, plan%eval_term_count
      monomial = xpow(plan%eval_exp(1, term_idx))*ypow(plan%eval_exp(2, term_idx))*zpow(plan%eval_exp(3, term_idx))* &
                 plan%eval_inv_factorial(term_idx)
      e(1) = e(1) - local_coeff(plan%eval_deriv_idx(1, term_idx))*monomial
      e(2) = e(2) - local_coeff(plan%eval_deriv_idx(2, term_idx))*monomial
      e(3) = e(3) - local_coeff(plan%eval_deriv_idx(3, term_idx))*monomial
    end do
  end subroutine eval_local_field_from_coeff

  subroutine build_powers(d, order, xpow, ypow, zpow)
    real(dp), intent(in) :: d(3)
    integer(i32), intent(in) :: order
    real(dp), intent(out) :: xpow(0:order), ypow(0:order), zpow(0:order)
    integer(i32) :: k

    xpow(0) = 1.0d0
    ypow(0) = 1.0d0
    zpow(0) = 1.0d0
    do k = 1_i32, order
      xpow(k) = xpow(k - 1_i32)*d(1)
      ypow(k) = ypow(k - 1_i32)*d(2)
      zpow(k) = zpow(k - 1_i32)*d(3)
    end do
  end subroutine build_powers

end program test_coulomb_fmm_core
