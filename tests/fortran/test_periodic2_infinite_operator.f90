program test_periodic2_infinite_operator
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, build_plan, update_state, &
                                  eval_point, eval_potential_point, destroy_plan, destroy_state
  use bem_coulomb_fmm_periodic_ewald, only: add_periodic2_exact_ewald_kneq0_correction_single_source, &
                                            add_periodic2_kneq0_potential_correction_single_source
  use test_support, only: test_begin, test_end, test_summary, assert_true
  implicit none

  type(fmm_options_type) :: options
  type(fmm_plan_type) :: plan0, plan1
  type(fmm_state_type) :: state0, state1
  real(dp) :: src_pos(3, 4), q(4), targets(3, 3)
  real(dp) :: field0(3), field1(3), field_ref(3), phi0, phi1, phi_ref
  real(dp) :: spectral_field(3), spectral_phi, teacher_field_error, teacher_potential_error
  real(dp) :: field_error, potential_error, max_field_error, max_potential_error
  real(dp) :: layer_field_delta, layer_potential_delta, translation_field(3), translation_phi
  integer(i32) :: target_idx
  character(len=512) :: cache_path0, cache_path1

  src_pos(:, 1) = [0.15_dp, 0.20_dp, -0.10_dp]
  src_pos(:, 2) = [0.75_dp, 0.25_dp, 0.08_dp]
  src_pos(:, 3) = [0.30_dp, 0.80_dp, 0.12_dp]
  src_pos(:, 4) = [0.82_dp, 0.72_dp, -0.06_dp]
  q = [1.0_dp, -0.8_dp, 0.6_dp, -0.8_dp]
  targets(:, 1) = [0.12_dp, 0.65_dp, -0.22_dp]
  targets(:, 2) = [0.55_dp, 0.48_dp, 0.18_dp]
  targets(:, 3) = [0.91_dp, 0.08_dp, 0.31_dp]
  call configure_options(options)

  call test_begin('cached_kneq0_matches_reference_and_image_layers')
  options%periodic_image_layers = 1_i32
  call build_plan(plan0, src_pos, options)
  call update_state(plan0, state0, q)
  cache_path0 = plan0%periodic_cache_path
  options%periodic_image_layers = 2_i32
  call build_plan(plan1, src_pos, options)
  call update_state(plan1, state1, q)
  cache_path1 = plan1%periodic_cache_path
  max_field_error = 0.0_dp
  max_potential_error = 0.0_dp
  layer_field_delta = 0.0_dp
  layer_potential_delta = 0.0_dp
  do target_idx = 1_i32, size(targets, 2)
    call eval_point(plan0, state0, targets(:, target_idx), field0)
    call eval_potential_point(plan0, state0, targets(:, target_idx), phi0)
    call eval_point(plan1, state1, targets(:, target_idx), field1)
    call eval_potential_point(plan1, state1, targets(:, target_idx), phi1)
    call exact_kneq0_reference(plan1, q, targets(:, target_idx), field_ref, phi_ref)
    field_error = sqrt(sum((field1 - field_ref)**2))/max(1.0e-10_dp, sqrt(sum(field_ref*field_ref)))
    potential_error = abs(phi1 - phi_ref)/max(1.0e-10_dp, abs(phi_ref))
    max_field_error = max(max_field_error, field_error)
    max_potential_error = max(max_potential_error, potential_error)
    layer_field_delta = max( &
                        layer_field_delta, sqrt(sum((field1 - field0)**2))/max(1.0e-10_dp, sqrt(sum(field_ref*field_ref))) &
                        )
    layer_potential_delta = max(layer_potential_delta, abs(phi1 - phi0)/max(1.0e-10_dp, abs(phi_ref)))
  end do
  write (*, '(a,4(es12.5,1x))') &
    'cached kneq0 errors(field,potential,layer-field,layer-potential)=', &
    max_field_error, max_potential_error, layer_field_delta, layer_potential_delta
  call assert_true(max_field_error < 8.0e-2_dp, 'cached kneq0 field error exceeds 8e-2')
  call assert_true(max_potential_error < 8.0e-2_dp, 'cached kneq0 potential error exceeds 8e-2')
  call assert_true(layer_field_delta < 8.0e-2_dp, 'cached kneq0 field changes with image layers')
  call assert_true(layer_potential_delta < 8.0e-2_dp, 'cached kneq0 potential gauge changes with image layers')
  q(4) = -0.6_dp
  call update_state(plan1, state1, q)
  do target_idx = 1_i32, size(targets, 2)
    call eval_point(plan1, state1, targets(:, target_idx), field1)
    call eval_potential_point(plan1, state1, targets(:, target_idx), phi1)
    call exact_kneq0_reference(plan1, q, targets(:, target_idx), field_ref, phi_ref)
    field_error = sqrt(sum((field1 - field_ref)**2))/max(1.0e-10_dp, sqrt(sum(field_ref*field_ref)))
    potential_error = abs(phi1 - phi_ref)/max(1.0e-10_dp, abs(phi_ref))
    max_field_error = max(max_field_error, field_error)
    max_potential_error = max(max_potential_error, potential_error)
  end do
  call assert_true(max_field_error < 8.0e-2_dp, 'non-neutral cached kneq0 field error exceeds 8e-2')
  call assert_true(max_potential_error < 8.0e-2_dp, 'non-neutral cached kneq0 potential error exceeds 8e-2')
  call exact_kneq0_reference(plan1, q, targets(:, 2), field_ref, phi_ref)
  call spectral_point_reference(src_pos, q, targets(:, 2), 16_i32, spectral_field, spectral_phi)
  teacher_field_error = sqrt(sum((field_ref - spectral_field)**2))/max(1.0e-10_dp, sqrt(sum(spectral_field**2)))
  teacher_potential_error = abs(phi_ref - spectral_phi)/max(1.0e-10_dp, abs(spectral_phi))
  write (*, '(a,2(es12.5,1x))') 'Ewald kneq0 teacher errors(field,potential)=', &
    teacher_field_error, teacher_potential_error
  call assert_true(teacher_field_error < 8.0e-2_dp, 'Ewald kneq0 field teacher differs from spectral reference')
  call assert_true(teacher_potential_error < 8.0e-2_dp, 'Ewald kneq0 potential teacher differs from spectral reference')
  call test_end()

  call test_begin('periodic_translation_invariance')
  call eval_point(plan1, state1, targets(:, 2) + [1.0_dp, -1.0_dp, 0.0_dp], translation_field)
  call eval_potential_point(plan1, state1, targets(:, 2) + [1.0_dp, -1.0_dp, 0.0_dp], translation_phi)
  call eval_point(plan1, state1, targets(:, 2), field1)
  call eval_potential_point(plan1, state1, targets(:, 2), phi1)
  call assert_true(sqrt(sum((translation_field - field1)**2)) < 1.0e-12_dp*max(1.0_dp, sqrt(sum(field1*field1))), &
                   'cached kneq0 field is not periodic')
  call assert_true(abs(translation_phi - phi1) < 1.0e-12_dp*max(1.0_dp, abs(phi1)), &
                   'cached kneq0 potential is not periodic')
  call test_end()

  call destroy_state(state0)
  call destroy_state(state1)
  call destroy_plan(plan0)
  call destroy_plan(plan1)
  call delete_cache(cache_path0)
  call delete_cache(cache_path1)
  call test_summary()

contains

  subroutine configure_options(config)
    type(fmm_options_type), intent(out) :: config
    config%theta = 0.4_dp
    config%leaf_max = 2_i32
    config%order = 4_i32
    config%use_periodic2 = .true.
    config%periodic_axes = [1_i32, 2_i32]
    config%periodic_len = [1.0_dp, 1.0_dp]
    config%periodic_far_correction = 'cached_kneq0'
    config%periodic_ewald_layers = 3_i32
    config%periodic_cache_dir = '.'
    config%periodic_generation_tolerance = 1.0e-8_dp
    config%target_box_min = [0.0_dp, 0.0_dp, -0.5_dp]
    config%target_box_max = [1.0_dp, 1.0_dp, 0.5_dp]
  end subroutine configure_options

  subroutine exact_kneq0_reference(plan, charge, target, field, potential)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: charge(:), target(3)
    real(dp), intent(out) :: field(3), potential
    real(dp) :: source(3), delta(3), r2, inv_r3
    integer(i32) :: source_idx, image1, image2

    field = 0.0_dp
    potential = 0.0_dp
    do source_idx = 1_i32, plan%nsrc
      do image1 = 1_i32, size(plan%shift_axis1)
        do image2 = 1_i32, size(plan%shift_axis2)
          source = plan%src_pos(:, source_idx)
          source(plan%options%periodic_axes(1)) = source(plan%options%periodic_axes(1)) + plan%shift_axis1(image1)
          source(plan%options%periodic_axes(2)) = source(plan%options%periodic_axes(2)) + plan%shift_axis2(image2)
          delta = target - source
          r2 = sum(delta*delta)
          inv_r3 = 1.0_dp/(sqrt(r2)*r2)
          field = field + charge(source_idx)*inv_r3*delta
          potential = potential + charge(source_idx)/sqrt(r2)
        end do
      end do
      call add_periodic2_exact_ewald_kneq0_correction_single_source( &
        plan, charge(source_idx), plan%src_pos(:, source_idx), target, field &
        )
      call add_periodic2_kneq0_potential_correction_single_source( &
        plan, charge(source_idx), plan%src_pos(:, source_idx), target, potential &
        )
    end do
  end subroutine exact_kneq0_reference

  subroutine spectral_point_reference(positions, charge, target, mode_layers, field, potential)
    real(dp), intent(in) :: positions(:, :), charge(:), target(3)
    integer(i32), intent(in) :: mode_layers
    real(dp), intent(out) :: field(3), potential
    integer(i32) :: source_idx, h1, h2
    real(dp) :: k1, k2, kmag, phase, dz, decay, pref
    real(dp), parameter :: two_pi = 2.0_dp*acos(-1.0_dp)

    field = 0.0_dp
    potential = 0.0_dp
    do source_idx = 1_i32, size(charge)
      dz = target(3) - positions(3, source_idx)
      do h1 = -mode_layers, mode_layers
        k1 = two_pi*real(h1, dp)
        do h2 = -mode_layers, mode_layers
          if (h1 == 0_i32 .and. h2 == 0_i32) cycle
          k2 = two_pi*real(h2, dp)
          kmag = sqrt(k1*k1 + k2*k2)
          phase = k1*(target(1) - positions(1, source_idx)) + k2*(target(2) - positions(2, source_idx))
          decay = exp(-kmag*abs(dz))
          pref = two_pi*charge(source_idx)*decay
          potential = potential + pref*cos(phase)/kmag
          field(1) = field(1) + pref*k1*sin(phase)/kmag
          field(2) = field(2) + pref*k2*sin(phase)/kmag
          if (dz > 0.0_dp) then
            field(3) = field(3) + pref*cos(phase)
          else if (dz < 0.0_dp) then
            field(3) = field(3) - pref*cos(phase)
          end if
        end do
      end do
    end do
  end subroutine spectral_point_reference

  subroutine delete_cache(path)
    character(len=*), intent(in) :: path
    call delete_file(path)
    call delete_file(trim(path)//'.lock')
  end subroutine delete_cache

  subroutine delete_file(path)
    character(len=*), intent(in) :: path
    logical :: exists
    integer :: unit
    inquire (file=trim(path), exist=exists)
    if (.not. exists) return
    open (newunit=unit, file=trim(path), status='old')
    close (unit, status='delete')
  end subroutine delete_file
end program test_periodic2_infinite_operator
