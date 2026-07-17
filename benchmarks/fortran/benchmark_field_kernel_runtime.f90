program benchmark_field_kernel_runtime
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type
  use bem_app_config, only: app_config, load_app_config, build_mesh_from_config
  use bem_field_solver, only: field_solver_type
  use bem_performance_profile, only: perf_wall_time_seconds
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer(i32), parameter :: default_target_count = 2048_i32
  type(app_config) :: cfg
  type(mesh_type) :: mesh
  type(field_solver_type) :: solver
  character(len=512) :: config_path, argument
  integer(i32) :: target_count, i
  integer :: argument_status
  real(dp), allocatable :: volume_targets(:, :), near_targets(:, :)
  real(dp) :: t0, mesh_seconds, init_seconds, refresh_seconds
  real(dp) :: e_volume_seconds, e_near_seconds, phi_volume_seconds, phi_near_seconds
  real(dp) :: e_volume_checksum, e_near_checksum, phi_volume_checksum, phi_near_checksum
  real(dp) :: warmup_field(3), warmup_potential

  call get_command_argument(1, config_path, status=argument_status)
  if (argument_status /= 0 .or. len_trim(config_path) == 0) then
    error stop 'usage: benchmark_field_kernel_runtime <beach.toml> [target_count]'
  end if
  target_count = default_target_count
  call get_command_argument(2, argument, status=argument_status)
  if (argument_status == 0 .and. len_trim(argument) > 0) then
    read (argument, *) target_count
  end if
  if (target_count <= 0_i32) error stop 'target_count must be positive'

  call load_app_config(trim(config_path), cfg)

  t0 = perf_wall_time_seconds()
  call build_mesh_from_config(cfg, mesh)
  mesh_seconds = perf_wall_time_seconds() - t0

  t0 = perf_wall_time_seconds()
  call solver%init(mesh, cfg%sim, cfg%field, cfg%periodic2, cfg%panel)
  init_seconds = perf_wall_time_seconds() - t0

  do i = 1_i32, mesh%nelem
    mesh%q_elem(i) = merge(1.0_dp, -1.0_dp, modulo(i, 2_i32) == 0_i32)* &
                     1.0e-17_dp*(1.0_dp + 0.25_dp*sin(real(i, dp)))
  end do
  t0 = perf_wall_time_seconds()
  call solver%refresh(mesh)
  refresh_seconds = perf_wall_time_seconds() - t0

  allocate (volume_targets(3, target_count), near_targets(3, target_count))
  call build_targets(cfg, mesh, volume_targets, near_targets)

  call solver%eval_e(mesh, volume_targets(:, 1), warmup_field)
  call solver%eval_potential(mesh, cfg%sim, volume_targets(:, 1), warmup_potential)
  call time_field_evaluations(solver, mesh, volume_targets, e_volume_seconds, e_volume_checksum)
  call time_field_evaluations(solver, mesh, near_targets, e_near_seconds, e_near_checksum)
  call time_potential_evaluations(solver, mesh, cfg, volume_targets, phi_volume_seconds, phi_volume_checksum)
  call time_potential_evaluations(solver, mesh, cfg, near_targets, phi_near_seconds, phi_near_checksum)

  if (.not. all(ieee_is_finite([mesh_seconds, init_seconds, refresh_seconds, &
                                e_volume_seconds, e_near_seconds, phi_volume_seconds, phi_near_seconds, &
                                e_volume_checksum, e_near_checksum, phi_volume_checksum, phi_near_checksum]))) then
    error stop 'benchmark produced a non-finite timing or checksum'
  end if

  write (*, '(a)') 'source_model,nelem,target_count,mesh_s,init_s,refresh_s,'// &
    'eval_e_volume_s,eval_e_near_s,eval_phi_volume_s,eval_phi_near_s,'// &
    'e_volume_checksum,e_near_checksum,phi_volume_checksum,phi_near_checksum'
  write (*, '(a,a,a,i0,a,i0,11(a,es24.16))') 'field_kernel_benchmark_csv=', trim(cfg%panel%source_model), ',', &
    mesh%nelem, ',', target_count, ',', mesh_seconds, ',', init_seconds, ',', refresh_seconds, ',', &
    e_volume_seconds, ',', e_near_seconds, ',', phi_volume_seconds, ',', phi_near_seconds, ',', &
    e_volume_checksum, ',', e_near_checksum, ',', phi_volume_checksum, ',', phi_near_checksum

contains

  subroutine build_targets(app, local_mesh, volume, near)
    type(app_config), intent(in) :: app
    type(mesh_type), intent(in) :: local_mesh
    real(dp), intent(out) :: volume(:, :), near(:, :)
    real(dp), parameter :: alpha(3) = [0.6180339887498949_dp, 0.4142135623730950_dp, 0.7320508075688772_dp]
    real(dp) :: box_span(3), fraction(3)
    integer(i32) :: target_idx, element_idx

    box_span = app%sim%box_max - app%sim%box_min
    do target_idx = 1_i32, int(size(volume, 2), i32)
      fraction = modulo(0.5_dp + alpha*real(target_idx, dp), 1.0_dp)
      fraction(3) = 0.01_dp + 0.98_dp*fraction(3)
      volume(:, target_idx) = app%sim%box_min + box_span*fraction

      element_idx = 1_i32 + modulo((target_idx - 1_i32)*7919_i32, local_mesh%nelem)
      near(:, target_idx) = local_mesh%centers(:, element_idx) + &
                            0.05_dp*local_mesh%h_elem(element_idx)*local_mesh%normals(:, element_idx)
    end do
  end subroutine build_targets

  subroutine time_field_evaluations(local_solver, local_mesh, targets, seconds, checksum)
    type(field_solver_type), intent(inout) :: local_solver
    type(mesh_type), intent(in) :: local_mesh
    real(dp), intent(in) :: targets(:, :)
    real(dp), intent(out) :: seconds, checksum
    real(dp) :: field(3), start_time
    integer(i32) :: target_idx

    checksum = 0.0_dp
    start_time = perf_wall_time_seconds()
    do target_idx = 1_i32, int(size(targets, 2), i32)
      call local_solver%eval_e(local_mesh, targets(:, target_idx), field)
      checksum = checksum + sum(abs(field))
    end do
    seconds = perf_wall_time_seconds() - start_time
  end subroutine time_field_evaluations

  subroutine time_potential_evaluations(local_solver, local_mesh, app, targets, seconds, checksum)
    type(field_solver_type), intent(inout) :: local_solver
    type(mesh_type), intent(in) :: local_mesh
    type(app_config), intent(in) :: app
    real(dp), intent(in) :: targets(:, :)
    real(dp), intent(out) :: seconds, checksum
    real(dp) :: potential, start_time
    integer(i32) :: target_idx

    checksum = 0.0_dp
    start_time = perf_wall_time_seconds()
    do target_idx = 1_i32, int(size(targets, 2), i32)
      call local_solver%eval_potential(local_mesh, app%sim, targets(:, target_idx), potential)
      checksum = checksum + abs(potential)
    end do
    seconds = perf_wall_time_seconds() - start_time
  end subroutine time_potential_evaluations
end program benchmark_field_kernel_runtime
