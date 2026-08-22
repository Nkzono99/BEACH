!> 光線追跡による光電子放出位置と速度のsampling。
module bem_photoelectron_injection
  use, intrinsic :: iso_fortran_env, only: error_unit
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_types, only: mesh_type, sim_config, hit_info, bc_periodic
  use bem_boundary, only: apply_box_boundary
  use bem_collision, only: collision_query_grid_stalled, collision_query_image_limit, &
                           collision_query_index_range, collision_query_invalid_segment, collision_query_ok, find_first_hit
  use bem_string_utils, only: lower_ascii
  use bem_injection_flux, only: sample_flux_weighted_normal_component
  use bem_injection_geometry, only: compute_face_area_from_bounds, resolve_face_axes, resolve_face_geometry
  use bem_injection_random, only: sample_standard_normal
  implicit none
  private

  real(dp), parameter :: default_velocity_sigma_cutoff = 6.0_dp

  public :: sample_photo_raycast_particles

contains

  !> 光線を注入面からレイキャストし、最初の命中要素から光電子を放出する。
  !! @param[in] mesh 交差判定に使うメッシュ。
  !! @param[in] sim ボックス境界条件とバッチ時間を含むシミュレーション設定。
  !! @param[in] inject_face 照射面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[in] pos_low 照射開口の下限座標 `(x,y,z)` [m]。
  !! @param[in] pos_high 照射開口の上限座標 `(x,y,z)` [m]。
  !! @param[in] ray_direction レイ進行方向ベクトル（正規化前でも可）。
  !! @param[in] m_particle 粒子1個あたりの質量 [kg]。
  !! @param[in] temperature_k 放出温度 [K]。
  !! @param[in] normal_drift_speed 放出法線方向のシフト速度 [m/s]。
  !! @param[in] emit_current_density_a_m2 レイ垂直面基準の放出電流面密度 [A/m^2]。
  !! @param[in] q_particle 粒子1個あたりの電荷 [C]。
  !! @param[in] rays_per_batch このrankで発射するレイ本数。
  !! @param[in] global_rays_per_batch 全rank合計のレイ本数（省略時は `rays_per_batch`）。
  !! @param[in] vmin_normal 放出法線速度の下限 [m/s]（省略時は 0）。
  !! @param[out] x 放出位置配列 `x(3,rays_per_batch)` [m]。
  !! @param[out] v 放出速度配列 `v(3,rays_per_batch)` [m/s]。
  !! @param[out] w 各マクロ粒子重み `w(rays_per_batch)`。
  !! @param[out] n_emit 実際に放出された粒子数（`<= rays_per_batch`）。
  !! @param[out] emit_elem_idx 放出元要素ID `emit_elem_idx(rays_per_batch)`（省略可）。
  !! @param[out] collision_failure_status 不完全な衝突照会の status（省略時は OpenMP 終了後に停止）。
  !! @param[out] collision_failure_ray 不完全な照会を返した最小 ray index。
  !! @param[out] collision_failure_bounce 不完全な照会を返した bounce index。
  subroutine sample_photo_raycast_particles( &
    mesh, sim, inject_face, pos_low, pos_high, ray_direction, m_particle, temperature_k, normal_drift_speed, &
    emit_current_density_a_m2, q_particle, rays_per_batch, x, v, w, n_emit, emit_elem_idx, global_rays_per_batch, &
    vmin_normal, collision_failure_status, collision_failure_ray, collision_failure_bounce &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: pos_low(3), pos_high(3)
    real(dp), intent(in) :: ray_direction(3)
    real(dp), intent(in) :: m_particle, temperature_k, normal_drift_speed
    real(dp), intent(in) :: emit_current_density_a_m2, q_particle
    integer(i32), intent(in) :: rays_per_batch
    real(dp), intent(out) :: x(:, :)
    real(dp), intent(out) :: v(:, :)
    real(dp), intent(out) :: w(:)
    integer(i32), intent(out) :: n_emit
    integer(i32), intent(out), optional :: emit_elem_idx(:)
    integer(i32), intent(in), optional :: global_rays_per_batch
    real(dp), intent(in), optional :: vmin_normal
    integer(i32), intent(out), optional :: collision_failure_status, collision_failure_ray, collision_failure_bounce

    real(dp), parameter :: eps = 1.0d-12
    integer(i32) :: i, total_rays, bounce_count, collision_status
    integer(i32) :: query_failure_status, query_failure_ray, query_failure_bounce
    integer :: axis_n, axis_t1, axis_t2
    real(dp) :: boundary_value, inward_normal(3), launch_dir(3), launch_dir_norm, inward_dot
    real(dp) :: launch_area, projected_area, w_hit, sigma
    real(dp) :: ray_pos(3), ray_dir(3), seg_end(3), boundary_probe(3), boundary_dir(3)
    real(dp) :: surf_normal(3), tangent1(3), tangent2(3)
    real(dp), allocatable :: u(:, :), hit_pos(:, :), hit_normal(:, :)
    integer(i32), allocatable :: hit_elem(:)
    logical, allocatable :: ray_emitted(:)
    logical :: reached_boundary, alive, escaped_boundary
    logical :: use_periodic2_mode
    type(hit_info) :: hit

    if (size(x, 1) /= 3 .or. size(v, 1) /= 3) error stop "photo_raycast particle arrays must have first dimension 3"
    if (size(x, 2) < rays_per_batch .or. size(v, 2) < rays_per_batch) then
      error stop "photo_raycast x/v arrays are smaller than rays_per_batch"
    end if
    if (size(w) < rays_per_batch) error stop "photo_raycast w array is smaller than rays_per_batch"
    if (present(emit_elem_idx)) then
      if (size(emit_elem_idx) < rays_per_batch) error stop "photo_raycast emit_elem_idx is smaller than rays_per_batch"
      emit_elem_idx = -1_i32
    end if
    if (rays_per_batch <= 0_i32) error stop "rays_per_batch must be > 0"
    total_rays = rays_per_batch
    if (present(global_rays_per_batch)) total_rays = global_rays_per_batch
    if (total_rays <= 0_i32) error stop "global_rays_per_batch must be > 0"
    if (.not. sim%use_box) error stop "photo_raycast requires sim.use_box = true"
    if (sim%batch_duration <= 0.0_dp) error stop "photo_raycast requires sim.batch_duration > 0"
    if (m_particle <= 0.0_dp) error stop "m_particle must be > 0"
    if (temperature_k < 0.0_dp) error stop "temperature_k must be >= 0"
    if (emit_current_density_a_m2 <= 0.0_dp) error stop "emit_current_density_a_m2 must be > 0"
    if (abs(q_particle) <= 0.0_dp) error stop "q_particle must be non-zero"

    call resolve_face_geometry(sim%box_min, sim%box_max, inject_face, axis_n, boundary_value, inward_normal)
    call resolve_face_axes(inject_face, axis_t1, axis_t2)

    launch_dir = ray_direction
    launch_dir_norm = sqrt(sum(launch_dir*launch_dir))
    if (launch_dir_norm <= 0.0_dp) error stop "ray_direction norm must be > 0"
    launch_dir = launch_dir/launch_dir_norm
    inward_dot = dot_product(launch_dir, inward_normal)
    if (inward_dot <= 0.0_dp) error stop "ray_direction must point inward from inject_face"

    launch_area = compute_face_area_from_bounds(inject_face, pos_low, pos_high)
    if (launch_area <= 0.0_dp) error stop "photo_raycast opening area must be positive"
    projected_area = launch_area*abs(inward_dot)
    w_hit = emit_current_density_a_m2*projected_area*sim%batch_duration/(abs(q_particle)*real(total_rays, dp))
    if (w_hit <= 0.0_dp) error stop "photo_raycast produced invalid w_hit"
    sigma = sqrt(k_boltzmann*temperature_k/m_particle)

    n_emit = 0_i32
    x = 0.0_dp
    v = 0.0_dp
    w = 0.0_dp
    query_failure_status = collision_query_ok
    query_failure_ray = huge(0_i32)
    query_failure_bounce = 0_i32
    if (present(collision_failure_status)) collision_failure_status = collision_query_ok
    if (present(collision_failure_ray)) collision_failure_ray = query_failure_ray
    if (present(collision_failure_bounce)) collision_failure_bounce = query_failure_bounce

    use_periodic2_mode = trim(lower_ascii(sim%field_bc_mode)) == 'periodic2'

    allocate (u(2, rays_per_batch))
    allocate (hit_pos(3, rays_per_batch), hit_normal(3, rays_per_batch), hit_elem(rays_per_batch), ray_emitted(rays_per_batch))
    call random_number(u)
    hit_pos = 0.0_dp
    hit_normal = 0.0_dp
    hit_elem = -1_i32
    ray_emitted = .false.

    !$omp parallel do default(none) schedule(static) &
    !$omp shared(rays_per_batch, axis_n, axis_t1, axis_t2, boundary_value, pos_low, pos_high, u, launch_dir, sim, mesh, &
    !$omp        use_periodic2_mode, ray_emitted, hit_pos, hit_normal, hit_elem, query_failure_status, &
    !$omp        query_failure_ray, query_failure_bounce) &
    !$omp private(i, ray_pos, ray_dir, seg_end, reached_boundary, alive, bounce_count, hit, surf_normal, boundary_probe, &
    !$omp         boundary_dir, escaped_boundary, collision_status)
    do i = 1_i32, rays_per_batch
      ray_pos = 0.0_dp
      ray_pos(axis_n) = boundary_value
      ray_pos(axis_t1) = pos_low(axis_t1) + (pos_high(axis_t1) - pos_low(axis_t1))*u(1, i)
      ray_pos(axis_t2) = pos_low(axis_t2) + (pos_high(axis_t2) - pos_low(axis_t2))*u(2, i)
      ray_dir = launch_dir
      ray_pos = ray_pos + ray_dir*eps

      alive = .true.
      bounce_count = 0_i32
      do while (alive .and. bounce_count <= sim%raycast_max_bounce)
        call step_ray_to_boundary(sim%box_min, sim%box_max, ray_pos, ray_dir, seg_end, reached_boundary)
        if (.not. reached_boundary) exit

        call find_first_hit( &
          mesh, ray_pos, seg_end, hit, sim=sim, box_min=sim%box_min, box_max=sim%box_max, &
          require_elem_inside=.true., status=collision_status &
          )
        if (collision_status /= collision_query_ok) then
          !$omp critical (beach_photo_collision_query_failure)
          if (query_failure_status == collision_query_ok .or. i < query_failure_ray .or. &
              (i == query_failure_ray .and. bounce_count < query_failure_bounce)) then
            query_failure_status = collision_status
            query_failure_ray = i
            query_failure_bounce = bounce_count
          end if
          !$omp end critical (beach_photo_collision_query_failure)
          exit
        end if
        if (hit%has_hit) then
          surf_normal = mesh%normals(:, hit%elem_idx)
          if (dot_product(surf_normal, ray_dir) > 0.0_dp) surf_normal = -surf_normal
          if (use_periodic2_mode) then
            hit_pos(:, i) = hit%pos_wrapped + surf_normal*eps
            call canonicalize_periodic_emission_position(sim, hit_pos(:, i))
          else
            hit_pos(:, i) = hit%pos + surf_normal*eps
          end if
          hit_normal(:, i) = surf_normal
          hit_elem(i) = hit%elem_idx
          ray_emitted(i) = .true.
          exit
        end if

        boundary_probe = seg_end + ray_dir*eps
        boundary_dir = ray_dir
        escaped_boundary = .false.
        call apply_box_boundary(sim, boundary_probe, boundary_dir, alive, escaped_boundary)
        if (.not. alive) exit
        ray_dir = boundary_dir/sqrt(sum(boundary_dir*boundary_dir))
        ray_pos = boundary_probe + ray_dir*eps
        bounce_count = bounce_count + 1_i32
      end do
    end do
    !$omp end parallel do

    if (query_failure_status /= collision_query_ok) then
      call finalize_photo_collision_query( &
        query_failure_status, query_failure_ray, query_failure_bounce, &
        collision_failure_status, collision_failure_ray, collision_failure_bounce &
        )
      return
    end if

    do i = 1_i32, rays_per_batch
      if (.not. ray_emitted(i)) cycle
      if (n_emit >= int(size(w), i32)) error stop "photo_raycast emitted particle buffer overflow"
      n_emit = n_emit + 1_i32
      surf_normal = hit_normal(:, i)
      call build_tangent_basis(surf_normal, tangent1, tangent2)
      if (present(vmin_normal)) then
        call sample_photo_emission_velocity( &
          sigma, normal_drift_speed, surf_normal, tangent1, tangent2, v(:, n_emit), vmin_normal=vmin_normal &
          )
      else
        call sample_photo_emission_velocity(sigma, normal_drift_speed, surf_normal, tangent1, tangent2, v(:, n_emit))
      end if
      x(:, n_emit) = hit_pos(:, i)
      w(n_emit) = w_hit
      if (present(emit_elem_idx)) emit_elem_idx(n_emit) = hit_elem(i)
    end do
  end subroutine sample_photo_raycast_particles

  !> 法線offset後の周期軸をprimary cellのstrict interiorへ戻す。
  pure subroutine canonicalize_periodic_emission_position(sim, position)
    type(sim_config), intent(in) :: sim
    real(dp), intent(inout) :: position(3)
    integer(i32) :: axis
    real(dp) :: span, scale, inset

    do axis = 1_i32, 3_i32
      if (sim%bc_low(axis) /= bc_periodic .or. sim%bc_high(axis) /= bc_periodic) cycle
      span = sim%box_max(axis) - sim%box_min(axis)
      scale = max(abs(sim%box_min(axis)), abs(sim%box_max(axis)), span, tiny(1.0_dp))
      inset = max(64.0_dp*epsilon(1.0_dp)*scale, spacing(scale))
      inset = min(0.25_dp*span, inset)
      position(axis) = sim%box_min(axis) + modulo(position(axis) - sim%box_min(axis), span)
      position(axis) = min(max(position(axis), sim%box_min(axis) + inset), sim%box_max(axis) - inset)
    end do
  end subroutine canonicalize_periodic_emission_position

  !> photo ray の不完全な衝突照会を返し、status 未要求なら OpenMP 外で停止する。
  subroutine finalize_photo_collision_query( &
    query_status, query_ray, query_bounce, collision_failure_status, collision_failure_ray, collision_failure_bounce &
    )
    integer(i32), intent(in) :: query_status, query_ray, query_bounce
    integer(i32), intent(out), optional :: collision_failure_status, collision_failure_ray, collision_failure_bounce
    character(len=16) :: status_name
    character(len=256) :: failure_message

    if (present(collision_failure_status)) collision_failure_status = query_status
    if (present(collision_failure_ray)) collision_failure_ray = query_ray
    if (present(collision_failure_bounce)) collision_failure_bounce = query_bounce
    if (present(collision_failure_status)) return

    select case (query_status)
    case (collision_query_image_limit)
      status_name = 'image_limit'
    case (collision_query_index_range)
      status_name = 'index_range'
    case (collision_query_invalid_segment)
      status_name = 'invalid_segment'
    case (collision_query_grid_stalled)
      status_name = 'grid_stalled'
    case default
      status_name = 'unknown'
    end select
    write (failure_message, '(a,i0,a,i0,a,a,a,i0)') &
      'photo_raycast collision query incomplete: ray=', query_ray, ' bounce=', query_bounce, &
      ' status=', trim(status_name), ' code=', query_status
    write (error_unit, '(a)') trim(failure_message)
    flush (error_unit)
    error stop 1
  end subroutine finalize_photo_collision_query

  !> レイを現在位置から最初のボックス境界まで進める。
  !! @param[in] box_min ボックス下限座標 `(x,y,z)` [m]。
  !! @param[in] box_max ボックス上限座標 `(x,y,z)` [m]。
  !! @param[in] x0 レイの現在位置 [m]。
  !! @param[in] ray_dir レイ進行方向（単位ベクトル）。
  !! @param[out] x1 境界到達位置 [m]。
  !! @param[out] reached_boundary 境界到達位置が求まった場合 `.true.`。
  subroutine step_ray_to_boundary(box_min, box_max, x0, ray_dir, x1, reached_boundary)
    real(dp), intent(in) :: box_min(3), box_max(3)
    real(dp), intent(in) :: x0(3), ray_dir(3)
    real(dp), intent(out) :: x1(3)
    logical, intent(out) :: reached_boundary

    real(dp), parameter :: eps = 1.0d-14
    integer :: axis
    real(dp) :: t_axis, t_hit

    t_hit = huge(1.0_dp)
    do axis = 1, 3
      if (ray_dir(axis) > eps) then
        t_axis = (box_max(axis) - x0(axis))/ray_dir(axis)
      else if (ray_dir(axis) < -eps) then
        t_axis = (box_min(axis) - x0(axis))/ray_dir(axis)
      else
        cycle
      end if
      if (t_axis > eps .and. t_axis < t_hit) t_hit = t_axis
    end do

    if (t_hit >= huge(1.0_dp)*0.5_dp) then
      reached_boundary = .false.
      x1 = x0
      return
    end if

    reached_boundary = .true.
    x1 = x0 + ray_dir*t_hit
    x1 = min(box_max, max(box_min, x1))
  end subroutine step_ray_to_boundary

  !> 面法線ベクトルから接線2軸を構築する。
  !! @param[in] normal 法線ベクトル。
  !! @param[out] tangent1 第1接線ベクトル。
  !! @param[out] tangent2 第2接線ベクトル。
  subroutine build_tangent_basis(normal, tangent1, tangent2)
    real(dp), intent(in) :: normal(3)
    real(dp), intent(out) :: tangent1(3), tangent2(3)

    real(dp) :: n(3), ref(3), norm_n, norm_t1

    norm_n = sqrt(sum(normal*normal))
    if (norm_n <= 0.0_dp) error stop "surface normal norm must be > 0"
    n = normal/norm_n

    if (abs(n(1)) < 0.9_dp) then
      ref = [1.0_dp, 0.0_dp, 0.0_dp]
    else
      ref = [0.0_dp, 1.0_dp, 0.0_dp]
    end if

    tangent1 = cross3(n, ref)
    norm_t1 = sqrt(sum(tangent1*tangent1))
    if (norm_t1 <= 0.0_dp) error stop "failed to build tangent basis"
    tangent1 = tangent1/norm_t1
    tangent2 = cross3(n, tangent1)
  end subroutine build_tangent_basis

  !> 光電子放出速度を局所法線座標でサンプルする。
  !! @param[in] sigma 熱速度標準偏差 [m/s]。
  !! @param[in] normal_drift_speed 放出法線方向のシフト速度 [m/s]。
  !! @param[in] normal 放出法線ベクトル（単位化済み）。
  !! @param[in] tangent1 第1接線ベクトル（単位化済み）。
  !! @param[in] tangent2 第2接線ベクトル（単位化済み）。
  !! @param[in] vmin_normal 放出法線速度の下限 [m/s]（省略時は 0）。
  !! @param[out] vel サンプルした速度ベクトル [m/s]。
  subroutine sample_photo_emission_velocity(sigma, normal_drift_speed, normal, tangent1, tangent2, vel, vmin_normal)
    real(dp), intent(in) :: sigma, normal_drift_speed
    real(dp), intent(in) :: normal(3), tangent1(3), tangent2(3)
    real(dp), intent(out) :: vel(3)
    real(dp), intent(in), optional :: vmin_normal

    real(dp) :: vn(1), z(2, 1), vt1, vt2

    if (present(vmin_normal)) then
      call sample_flux_weighted_normal_component( &
        normal_drift_speed, sigma, vn, vmin_normal=vmin_normal, sigma_cutoff=default_velocity_sigma_cutoff &
        )
    else
      call sample_flux_weighted_normal_component(normal_drift_speed, sigma, vn, sigma_cutoff=default_velocity_sigma_cutoff)
    end if
    vt1 = 0.0_dp
    vt2 = 0.0_dp
    if (sigma > 0.0_dp) then
      call sample_standard_normal(z, sigma_cutoff=default_velocity_sigma_cutoff)
      vt1 = sigma*z(1, 1)
      vt2 = sigma*z(2, 1)
    end if
    vel = normal*vn(1) + tangent1*vt1 + tangent2*vt2
  end subroutine sample_photo_emission_velocity

  !> 3次元外積を返す。
  pure function cross3(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross3

end module bem_photoelectron_injection
