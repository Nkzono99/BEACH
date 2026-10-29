!> 光線追跡による光電子放出位置と速度のsampling。
module bem_photoelectron_injection
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: k_boltzmann
  use bem_types, only: mesh_type, photo_raycast_config, hit_info, bc_periodic
  use bem_boundary, only: apply_ray_box_boundary
  use bem_collision, only: find_first_hit
  use bem_status_codes, only: collision_query_ok
  use bem_string_utils, only: lower_ascii
  use bem_injection_flux, only: sample_flux_weighted_normal_component
  use bem_injection_geometry, only: compute_face_area_from_bounds, resolve_face_axes, resolve_face_geometry
  use bem_injection_random, only: sample_standard_normal
  use bem_counter_rng, only: counter_uniform, counter_uniforms
  implicit none
  private

  real(dp), parameter :: default_velocity_sigma_cutoff = 6.0_dp

  public :: sample_photo_emission_velocity, build_tangent_basis
  public :: canonicalize_periodic_emission_position
  public :: cast_photo_ray, sample_photo_ray_slot, photo_ray_macro_weight

contains

  !> 面上の点 `ray_start` から `launch_dir` へ光線を飛ばし、box 境界の反射を経て最初に当たる panel を返す。
  !! @param[out] emitted 要素に当たったとき `.true.`。
  !! @param[out] collision_status 衝突照会の status（`collision_query_ok` 以外なら未完了）。
  !! @param[out] bounce_count 未完了時の bounce index。
  subroutine cast_photo_ray( &
    mesh, sim, ray_start, launch_dir, use_periodic2_mode, hit_pos, hit_normal, hit_elem, emitted, collision_status, &
    bounce_count &
    )
    type(mesh_type), intent(in) :: mesh
    class(photo_raycast_config), intent(in) :: sim
    real(dp), intent(in) :: ray_start(3), launch_dir(3)
    logical, intent(in) :: use_periodic2_mode
    real(dp), intent(out) :: hit_pos(3), hit_normal(3)
    integer(i32), intent(out) :: hit_elem, collision_status, bounce_count
    logical, intent(out) :: emitted
    real(dp), parameter :: eps = 1.0d-12
    real(dp) :: ray_pos(3), ray_dir(3), seg_end(3), boundary_probe(3), boundary_dir(3), surf_normal(3)
    logical :: reached_boundary, alive
    type(hit_info) :: hit

    hit_pos = 0.0_dp
    hit_normal = 0.0_dp
    hit_elem = -1_i32
    emitted = .false.
    collision_status = collision_query_ok
    ray_dir = launch_dir
    ray_pos = ray_start + ray_dir*eps
    alive = .true.
    bounce_count = 0_i32
    do while (alive .and. bounce_count <= sim%raycast_max_bounce)
      call step_ray_to_boundary(sim%box_min, sim%box_max, ray_pos, ray_dir, seg_end, reached_boundary)
      if (.not. reached_boundary) exit

      call find_first_hit( &
        mesh, ray_pos, seg_end, hit, sim=sim, box_min=sim%box_min, box_max=sim%box_max, &
        require_elem_inside=.true., status=collision_status &
        )
      if (collision_status /= collision_query_ok) return
      if (hit%has_hit) then
        surf_normal = mesh%normals(:, hit%elem_idx)
        if (dot_product(surf_normal, ray_dir) > 0.0_dp) surf_normal = -surf_normal
        if (use_periodic2_mode) then
          hit_pos = hit%pos_wrapped + surf_normal*eps
          call canonicalize_periodic_emission_position(sim, hit_pos)
        else
          hit_pos = hit%pos + surf_normal*eps
        end if
        hit_normal = surf_normal
        hit_elem = hit%elem_idx
        emitted = .true.
        return
      end if

      boundary_probe = seg_end + ray_dir*eps
      boundary_dir = ray_dir
      call apply_ray_box_boundary(sim, boundary_probe, boundary_dir, alive)
      if (.not. alive) exit
      ray_dir = boundary_dir/sqrt(sum(boundary_dir*boundary_dir))
      ray_pos = boundary_probe + ray_dir*eps
      bounce_count = bounce_count + 1_i32
    end do
  end subroutine cast_photo_ray

  !> 1 本の光線で放出されるマクロ粒子の重み（全 rank 合計 `total_rays` 本で放出電流を分担する）。
  pure real(dp) function photo_ray_macro_weight( &
    inject_face, pos_low, pos_high, ray_direction, emit_current_density_a_m2, q_particle, batch_duration, total_rays &
    ) result(w_hit)
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: pos_low(3), pos_high(3), ray_direction(3), emit_current_density_a_m2, q_particle
    real(dp), intent(in) :: batch_duration
    integer(i32), intent(in) :: total_rays
    real(dp) :: launch_dir(3), inward_dot, launch_area, projected_area
    integer :: axis_n, axis_t1, axis_t2

    call resolve_face_axes(inject_face, axis_t1, axis_t2)
    axis_n = 6 - axis_t1 - axis_t2
    launch_dir = ray_direction/sqrt(sum(ray_direction*ray_direction))
    inward_dot = launch_dir(axis_n)
    launch_area = compute_face_area_from_bounds(inject_face, pos_low, pos_high)
    projected_area = launch_area*abs(inward_dot)
    w_hit = emit_current_density_a_m2*projected_area*batch_duration/(abs(q_particle)*real(max(1_i32, total_rays), dp))
  end function photo_ray_macro_weight

  !> 1 本の光線を counter-based 乱数で発射し、当たれば光電子を 1 つ返す。
  !! 面内の発射位置は `ray_counter` の draw 1, 2、放出速度は `velocity_counter` で決まるので、どの rank・thread が
  !! 同じ slot を生成しても同じ粒子になる。
  subroutine sample_photo_ray_slot( &
    mesh, sim, inject_face, pos_low, pos_high, ray_direction, m_particle, temperature_k, normal_drift_speed, w_hit, &
    rng_seed, ray_counter, velocity_counter, x, v, w, emit_elem, emitted, collision_status, bounce_count, vmin_normal &
    )
    type(mesh_type), intent(in) :: mesh
    class(photo_raycast_config), intent(in) :: sim
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: pos_low(3), pos_high(3), ray_direction(3)
    real(dp), intent(in) :: m_particle, temperature_k, normal_drift_speed, w_hit
    integer(i32), intent(in) :: rng_seed
    integer(i64), intent(in) :: ray_counter(:), velocity_counter(:)
    real(dp), intent(out) :: x(3), v(3), w
    integer(i32), intent(out) :: emit_elem, collision_status, bounce_count
    logical, intent(out) :: emitted
    real(dp), intent(in), optional :: vmin_normal
    integer :: axis_n, axis_t1, axis_t2
    real(dp) :: boundary_value, inward_normal(3), launch_dir(3), ray_pos(3), hit_pos(3), hit_normal(3)
    real(dp) :: tangent1(3), tangent2(3), sigma, ray_uniform(2)
    logical :: use_periodic2_mode

    x = 0.0_dp
    v = 0.0_dp
    w = 0.0_dp
    emit_elem = -1_i32
    call resolve_face_geometry(sim%box_min, sim%box_max, inject_face, axis_n, boundary_value, inward_normal)
    call resolve_face_axes(inject_face, axis_t1, axis_t2)
    launch_dir = ray_direction/sqrt(sum(ray_direction*ray_direction))
    use_periodic2_mode = trim(lower_ascii(sim%field_bc_mode)) == 'periodic2'
    ray_pos = 0.0_dp
    ray_pos(axis_n) = boundary_value
    call counter_uniforms(rng_seed, ray_counter, 1_i64, ray_uniform)
    ray_pos(axis_t1) = pos_low(axis_t1) + (pos_high(axis_t1) - pos_low(axis_t1))*ray_uniform(1)
    ray_pos(axis_t2) = pos_low(axis_t2) + (pos_high(axis_t2) - pos_low(axis_t2))*ray_uniform(2)
    call cast_photo_ray( &
      mesh, sim, ray_pos, launch_dir, use_periodic2_mode, hit_pos, hit_normal, emit_elem, emitted, collision_status, &
      bounce_count &
      )
    if (collision_status /= collision_query_ok .or. .not. emitted) return
    sigma = sqrt(k_boltzmann*temperature_k/m_particle)
    call build_tangent_basis(hit_normal, tangent1, tangent2)
    if (present(vmin_normal)) then
      call sample_photo_emission_velocity( &
        sigma, normal_drift_speed, hit_normal, tangent1, tangent2, v, vmin_normal=vmin_normal, rng_seed=rng_seed, &
        rng_counter=velocity_counter &
        )
    else
      call sample_photo_emission_velocity( &
        sigma, normal_drift_speed, hit_normal, tangent1, tangent2, v, rng_seed=rng_seed, rng_counter=velocity_counter &
        )
    end if
    x = hit_pos
    w = w_hit
  end subroutine sample_photo_ray_slot

  !> 法線offset後の周期軸をprimary cellのstrict interiorへ戻す。
  pure subroutine canonicalize_periodic_emission_position(sim, position)
    class(photo_raycast_config), intent(in) :: sim
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
  !! @param[in] rng_seed counter-based 乱数の seed。
  !! @param[in] rng_counter この放出の counter（法線成分は draw 1、接線成分は draw 100 以降）。
  !! @param[out] vel サンプルした速度ベクトル [m/s]。
  subroutine sample_photo_emission_velocity( &
    sigma, normal_drift_speed, normal, tangent1, tangent2, vel, vmin_normal, rng_seed, rng_counter &
    )
    real(dp), intent(in) :: sigma, normal_drift_speed
    real(dp), intent(in) :: normal(3), tangent1(3), tangent2(3)
    real(dp), intent(out) :: vel(3)
    real(dp), intent(in), optional :: vmin_normal

    integer(i32), intent(in) :: rng_seed
    integer(i64), intent(in) :: rng_counter(:)

    real(dp) :: vn(1), z(2, 1), vt1, vt2

    call sample_flux_weighted_normal_component( &
      normal_drift_speed, sigma, vn, vmin_normal=vmin_normal, sigma_cutoff=default_velocity_sigma_cutoff, &
      uniform=[counter_uniform(rng_seed, rng_counter, 1_i64)] &
      )
    vt1 = 0.0_dp
    vt2 = 0.0_dp
    if (sigma > 0.0_dp) then
      call sample_standard_normal(z, sigma_cutoff=default_velocity_sigma_cutoff, rng_seed=rng_seed, rng_counter=rng_counter)
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
