!> 平面/穴あき平面/円板/リング/箱/円柱/球テンプレートから三角形メッシュを生成するユーティリティ。
module bem_templates
    use bem_kinds, only: dp, i32
    use bem_types, only: mesh_type
    use bem_mesh, only: init_mesh
    implicit none
contains

    !> XY平面を `nx*ny` 分割し、各セルを2三角形へ分割したメッシュを生成する。
  !! @param[out] mesh 生成した平面三角形メッシュ。
  !! @param[in] size_x X方向の平面サイズ [m]（省略時 1.0）。
  !! @param[in] size_y Y方向の平面サイズ [m]（省略時 1.0）。
  !! @param[in] nx X方向分割数（省略時 1）。
  !! @param[in] ny Y方向分割数（省略時 1）。
  !! @param[in] center 平面中心座標 `(x,y,z)` [m]（省略時原点）。
    subroutine make_plane(mesh, size_x, size_y, nx, ny, center)
        type(mesh_type), intent(out) :: mesh
        real(dp), intent(in), optional :: size_x, size_y
        integer(i32), intent(in), optional :: nx, ny
        real(dp), intent(in), optional :: center(3)
        real(dp) :: sx, sy, c(3), x0, y0, dx, dy
        integer(i32) :: nx0, ny0, ix, iy, itri, nelem
        real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)

        sx = 1.0d0; sy = 1.0d0; nx0 = 1; ny0 = 1; c = 0.0d0
        if (present(size_x)) sx = size_x
        if (present(size_y)) sy = size_y
        if (present(nx)) nx0 = nx
        if (present(ny)) ny0 = ny
        if (present(center)) c = center
        if (nx0 <= 0 .or. ny0 <= 0) error stop "nx and ny must be positive"

        nelem = 2*nx0*ny0
        allocate (v0(3, nelem), v1(3, nelem), v2(3, nelem))
        dx = sx/real(nx0, dp); dy = sy/real(ny0, dp)
        x0 = c(1) - 0.5d0*sx; y0 = c(2) - 0.5d0*sy

        itri = 0
        do ix = 0, nx0 - 1
            do iy = 0, ny0 - 1
                itri = itri + 1
                v0(:, itri) = [x0 + dx*ix, y0 + dy*iy, c(3)]
                v1(:, itri) = [x0 + dx*(ix + 1), y0 + dy*iy, c(3)]
                v2(:, itri) = [x0 + dx*(ix + 1), y0 + dy*(iy + 1), c(3)]
                itri = itri + 1
                v0(:, itri) = [x0 + dx*ix, y0 + dy*iy, c(3)]
                v1(:, itri) = [x0 + dx*(ix + 1), y0 + dy*(iy + 1), c(3)]
                v2(:, itri) = [x0 + dx*ix, y0 + dy*(iy + 1), c(3)]
            end do
        end do
        call init_mesh(mesh, v0, v1, v2)
    end subroutine make_plane

    !> XY平面上の円板を極座標分割し、外周へ向かって三角形化したメッシュを生成する。
  !! @param[out] mesh 生成した円板メッシュ。
  !! @param[in] radius 円板半径 [m]（省略時 0.5）。
  !! @param[in] n_theta 周方向分割数（省略時 24）。
  !! @param[in] n_r 半径方向分割数（省略時 8）。
  !! @param[in] center 円板中心座標 `(x,y,z)` [m]（省略時原点）。
    subroutine make_disk(mesh, radius, n_theta, n_r, center)
        type(mesh_type), intent(out) :: mesh
        real(dp), intent(in), optional :: radius, center(3)
        integer(i32), intent(in), optional :: n_theta, n_r
        real(dp) :: r, c(3)
        integer(i32) :: nt, nr0

        r = 0.5d0; nt = 24; nr0 = 8; c = 0.0d0
        if (present(radius)) r = radius
        if (present(n_theta)) nt = n_theta
        if (present(n_r)) nr0 = n_r
        if (present(center)) c = center
        call make_annulus(mesh, radius=r, inner_radius=0.0d0, n_theta=nt, n_r=nr0, center=c)
    end subroutine make_disk

    !> XY平面上の同心リングを極座標分割し、三角形メッシュを生成する。
  !! @param[out] mesh 生成したリング（環状）メッシュ。
  !! @param[in] radius 外半径 [m]（省略時 0.5）。
  !! @param[in] inner_radius 内半径 [m]（省略時 0.25）。
  !! @param[in] n_theta 周方向分割数（省略時 24）。
  !! @param[in] n_r 半径方向分割数（省略時 4）。
  !! @param[in] center リング中心座標 `(x,y,z)` [m]（省略時原点）。
    subroutine make_annulus(mesh, radius, inner_radius, n_theta, n_r, center)
        type(mesh_type), intent(out) :: mesh
        real(dp), intent(in), optional :: radius, inner_radius, center(3)
        integer(i32), intent(in), optional :: n_theta, n_r
        real(dp) :: r_out, r_in, c(3), pi, t0, t1, r0, r1, dr
        integer(i32) :: nt, nr0, ir, it, itri, nelem
        real(dp) :: p00(3), p01(3), p10(3), p11(3), center_p(3)
        real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)

        r_out = 0.5d0; r_in = 0.25d0; nt = 24; nr0 = 4; c = 0.0d0
        if (present(radius)) r_out = radius
        if (present(inner_radius)) r_in = inner_radius
        if (present(n_theta)) nt = n_theta
        if (present(n_r)) nr0 = n_r
        if (present(center)) c = center
        if (nt < 3 .or. nr0 <= 0) error stop "n_theta >= 3 and n_r > 0"
        if (r_out <= 0.0d0) error stop "radius must be > 0"
        if (r_in < 0.0d0) error stop "inner_radius must be >= 0"
        if (r_in >= r_out) error stop "inner_radius must be smaller than radius"

        if (r_in <= 0.0d0) then
            nelem = nt*(2*nr0 - 1)
        else
            nelem = 2*nt*nr0
        end if
        allocate (v0(3, nelem), v1(3, nelem), v2(3, nelem))

        pi = acos(-1.0d0)
        dr = (r_out - r_in)/real(nr0, dp)
        center_p = [c(1), c(2), c(3)]
        itri = 0
        do ir = 0, nr0 - 1
            r0 = r_in + dr*real(ir, dp)
            r1 = r_in + dr*real(ir + 1, dp)
            do it = 0, nt - 1
                t0 = 2.0d0*pi*real(it, dp)/real(nt, dp)
                t1 = 2.0d0*pi*real(mod(it + 1, nt), dp)/real(nt, dp)
                p00 = [c(1) + r0*cos(t0), c(2) + r0*sin(t0), c(3)]
                p01 = [c(1) + r0*cos(t1), c(2) + r0*sin(t1), c(3)]
                p10 = [c(1) + r1*cos(t0), c(2) + r1*sin(t0), c(3)]
                p11 = [c(1) + r1*cos(t1), c(2) + r1*sin(t1), c(3)]
                if (ir == 0 .and. r_in <= 0.0d0) then
                    call push_tri(v0, v1, v2, itri, center_p, p10, p11)
                else
                    call push_tri(v0, v1, v2, itri, p00, p10, p11)
                    call push_tri(v0, v1, v2, itri, p00, p11, p01)
                end if
            end do
        end do

        call init_mesh(mesh, v0, v1, v2)
    end subroutine make_annulus

    !> XY平面の長方形プレートから円形穴を除いたメッシュを生成する。
  !! 穴境界は `n_theta` 分割の多角形近似で表し、外周は長方形境界に一致させる。
  !! @param[out] mesh 生成した穴あきプレートメッシュ。
  !! @param[in] size_x X方向サイズ [m]（省略時 1.0）。
  !! @param[in] size_y Y方向サイズ [m]（省略時 1.0）。
  !! @param[in] radius 穴半径 [m]（省略時 0.2）。
  !! @param[in] n_theta 穴の周方向分割数（省略時 48）。
  !! @param[in] n_r 穴縁から外周までの半径方向分割数（省略時 4）。
  !! @param[in] center プレート中心座標 `(x,y,z)` [m]（省略時原点）。
    subroutine make_plate_hole(mesh, size_x, size_y, radius, n_theta, n_r, center)
        type(mesh_type), intent(out) :: mesh
        real(dp), intent(in), optional :: size_x, size_y, radius, center(3)
        integer(i32), intent(in), optional :: n_theta, n_r
        real(dp) :: sx, sy, r_hole, c(3), pi, angle, f0, f1
        real(dp) :: x_min, x_max, y_min, y_max, min_margin
        real(dp) :: p00(3), p01(3), p10(3), p11(3), corners(3, 3)
        integer(i32) :: nt, nr0, it, ir, itri, nelem, n_corner, k
        integer(i32), allocatable :: edge_ids(:)
        real(dp), allocatable :: hole_pts(:, :), outer_pts(:, :)
        real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)

        sx = 1.0d0; sy = 1.0d0; r_hole = 0.2d0; nt = 48; nr0 = 4; c = 0.0d0
        if (present(size_x)) sx = size_x
        if (present(size_y)) sy = size_y
        if (present(radius)) r_hole = radius
        if (present(n_theta)) nt = n_theta
        if (present(n_r)) nr0 = n_r
        if (present(center)) c = center
        if (sx <= 0.0d0 .or. sy <= 0.0d0) error stop "size_x and size_y must be positive"
        if (nt < 3 .or. nr0 <= 0) error stop "n_theta >= 3 and n_r > 0"

        x_min = c(1) - 0.5d0*sx
        x_max = c(1) + 0.5d0*sx
        y_min = c(2) - 0.5d0*sy
        y_max = c(2) + 0.5d0*sy
        min_margin = min(min(c(1) - x_min, x_max - c(1)), min(c(2) - y_min, y_max - c(2)))
        if (r_hole <= 0.0d0) error stop "radius must be > 0 for plate_hole"
        if (r_hole >= min_margin) error stop "plate_hole radius must be smaller than half-width/half-height"

        allocate (hole_pts(3, nt + 1), outer_pts(3, nt + 1), edge_ids(nt + 1))
        pi = acos(-1.0d0)
        do it = 0, nt
            angle = 2.0d0*pi*real(it, dp)/real(nt, dp)
            hole_pts(:, it + 1) = [c(1) + r_hole*cos(angle), c(2) + r_hole*sin(angle), c(3)]
            call ray_to_rectangle( &
              c(1), c(2), cos(angle), sin(angle), x_min, x_max, y_min, y_max, c(3), outer_pts(:, it + 1), edge_ids(it + 1) &
            )
        end do

        nelem = 2*nt*nr0
        do it = 1, nt
            nelem = nelem + transition_corner_count(edge_ids(it), edge_ids(it + 1))
        end do
        allocate (v0(3, nelem), v1(3, nelem), v2(3, nelem))

        itri = 0
        do it = 1, nt
            do ir = 0, nr0 - 1
                f0 = real(ir, dp)/real(nr0, dp)
                f1 = real(ir + 1, dp)/real(nr0, dp)
                p00 = hole_pts(:, it) + f0*(outer_pts(:, it) - hole_pts(:, it))
                p10 = hole_pts(:, it) + f1*(outer_pts(:, it) - hole_pts(:, it))
                p01 = hole_pts(:, it + 1) + f0*(outer_pts(:, it + 1) - hole_pts(:, it + 1))
                p11 = hole_pts(:, it + 1) + f1*(outer_pts(:, it + 1) - hole_pts(:, it + 1))
                call push_tri(v0, v1, v2, itri, p00, p10, p11)
                call push_tri(v0, v1, v2, itri, p00, p11, p01)
            end do

            call transition_corners( &
              edge_ids(it), edge_ids(it + 1), x_min, x_max, y_min, y_max, c(3), corners, n_corner &
            )
            if (n_corner > 0) then
                if (n_corner == 1) then
                    call push_tri(v0, v1, v2, itri, outer_pts(:, it), corners(:, 1), outer_pts(:, it + 1))
                else
                    call push_tri(v0, v1, v2, itri, outer_pts(:, it), corners(:, 1), corners(:, 2))
                    do k = 2, n_corner - 1
                        call push_tri(v0, v1, v2, itri, outer_pts(:, it), corners(:, k), corners(:, k + 1))
                    end do
                    call push_tri(v0, v1, v2, itri, outer_pts(:, it), corners(:, n_corner), outer_pts(:, it + 1))
                end if
            end if
        end do

        if (itri /= nelem) error stop "plate_hole triangle bookkeeping mismatch"
        call init_mesh(mesh, v0, v1, v2)
    end subroutine make_plate_hole

    !> 直方体6面を分割数に応じて三角形化し、外向き法線向きでメッシュを生成する。
  !! @param[out] mesh 生成した直方体表面メッシュ。
  !! @param[in] size 直方体サイズ `(sx,sy,sz)` [m]（省略時 `[1,1,1]`）。
  !! @param[in] center 直方体中心座標 `(x,y,z)` [m]（省略時原点）。
  !! @param[in] nx X方向分割数（省略時 1）。
  !! @param[in] ny Y方向分割数（省略時 1）。
  !! @param[in] nz Z方向分割数（省略時 1）。
    subroutine make_box(mesh, size, center, nx, ny, nz)
        type(mesh_type), intent(out) :: mesh
        real(dp), intent(in), optional :: size(3), center(3)
        integer(i32), intent(in), optional :: nx, ny, nz
        real(dp) :: s(3), c(3), hx, hy, hz, dx, dy, dz, x, y, z
        integer(i32) :: nx0, ny0, nz0, nelem, itri, ix, iy, iz
        real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)

        s = [1.0d0, 1.0d0, 1.0d0]; c = 0.0d0; nx0 = 1; ny0 = 1; nz0 = 1
        if (present(size)) s = size
        if (present(center)) c = center
        if (present(nx)) nx0 = nx
        if (present(ny)) ny0 = ny
        if (present(nz)) nz0 = nz
        if (nx0 <= 0 .or. ny0 <= 0 .or. nz0 <= 0) error stop "nx, ny, nz must be positive"

        hx = 0.5d0*s(1); hy = 0.5d0*s(2); hz = 0.5d0*s(3)
        dx = s(1)/real(nx0, dp); dy = s(2)/real(ny0, dp); dz = s(3)/real(nz0, dp)
        nelem = 4*(nx0*ny0 + ny0*nz0 + nx0*nz0)
        allocate (v0(3, nelem), v1(3, nelem), v2(3, nelem)); itri = 0

        do iz = 1, 2
            z = merge(c(3) - hz, c(3) + hz, iz == 1)
            do ix = 0, nx0 - 1
                do iy = 0, ny0 - 1
                    if (iz == 1) then
                        call push_tri(v0, v1, v2, itri, &
                                      [c(1) - hx + dx*ix, c(2) - hy + dy*iy, z], &
                                      [c(1) - hx + dx*(ix + 1), c(2) - hy + dy*(iy + 1), z], &
                                      [c(1) - hx + dx*(ix + 1), c(2) - hy + dy*iy, z])
                        call push_tri(v0, v1, v2, itri, &
                                      [c(1) - hx + dx*ix, c(2) - hy + dy*iy, z], &
                                      [c(1) - hx + dx*ix, c(2) - hy + dy*(iy + 1), z], &
                                      [c(1) - hx + dx*(ix + 1), c(2) - hy + dy*(iy + 1), z])
                    else
                        call push_tri(v0, v1, v2, itri, &
                                      [c(1) - hx + dx*ix, c(2) - hy + dy*iy, z], &
                                      [c(1) - hx + dx*(ix + 1), c(2) - hy + dy*iy, z], &
                                      [c(1) - hx + dx*(ix + 1), c(2) - hy + dy*(iy + 1), z])
                        call push_tri(v0, v1, v2, itri, &
                                      [c(1) - hx + dx*ix, c(2) - hy + dy*iy, z], &
                                      [c(1) - hx + dx*(ix + 1), c(2) - hy + dy*(iy + 1), z], &
                                      [c(1) - hx + dx*ix, c(2) - hy + dy*(iy + 1), z])
                    end if
                end do
            end do
        end do

        do ix = 1, 2
            x = merge(c(1) - hx, c(1) + hx, ix == 1)
            do iy = 0, ny0 - 1
                do iz = 0, nz0 - 1
                    if (ix == 1) then
                        call push_tri(v0, v1, v2, itri, &
                                      [x, c(2) - hy + dy*iy, c(3) - hz + dz*iz], &
                                      [x, c(2) - hy + dy*(iy + 1), c(3) - hz + dz*(iz + 1)], &
                                      [x, c(2) - hy + dy*(iy + 1), c(3) - hz + dz*iz])
                        call push_tri(v0, v1, v2, itri, &
                                      [x, c(2) - hy + dy*iy, c(3) - hz + dz*iz], &
                                      [x, c(2) - hy + dy*iy, c(3) - hz + dz*(iz + 1)], &
                                      [x, c(2) - hy + dy*(iy + 1), c(3) - hz + dz*(iz + 1)])
                    else
                        call push_tri(v0, v1, v2, itri, &
                                      [x, c(2) - hy + dy*iy, c(3) - hz + dz*iz], &
                                      [x, c(2) - hy + dy*(iy + 1), c(3) - hz + dz*iz], &
                                      [x, c(2) - hy + dy*(iy + 1), c(3) - hz + dz*(iz + 1)])
                        call push_tri(v0, v1, v2, itri, &
                                      [x, c(2) - hy + dy*iy, c(3) - hz + dz*iz], &
                                      [x, c(2) - hy + dy*(iy + 1), c(3) - hz + dz*(iz + 1)], &
                                      [x, c(2) - hy + dy*iy, c(3) - hz + dz*(iz + 1)])
                    end if
                end do
            end do
        end do

        do iy = 1, 2
            y = merge(c(2) - hy, c(2) + hy, iy == 1)
            do ix = 0, nx0 - 1
                do iz = 0, nz0 - 1
                    if (iy == 1) then
                        call push_tri(v0, v1, v2, itri, &
                                      [c(1) - hx + dx*ix, y, c(3) - hz + dz*iz], &
                                      [c(1) - hx + dx*(ix + 1), y, c(3) - hz + dz*iz], &
                                      [c(1) - hx + dx*(ix + 1), y, c(3) - hz + dz*(iz + 1)])
                        call push_tri(v0, v1, v2, itri, &
                                      [c(1) - hx + dx*ix, y, c(3) - hz + dz*iz], &
                                      [c(1) - hx + dx*(ix + 1), y, c(3) - hz + dz*(iz + 1)], &
                                      [c(1) - hx + dx*ix, y, c(3) - hz + dz*(iz + 1)])
                    else
                        call push_tri(v0, v1, v2, itri, &
                                      [c(1) - hx + dx*ix, y, c(3) - hz + dz*iz], &
                                      [c(1) - hx + dx*(ix + 1), y, c(3) - hz + dz*(iz + 1)], &
                                      [c(1) - hx + dx*(ix + 1), y, c(3) - hz + dz*iz])
                        call push_tri(v0, v1, v2, itri, &
                                      [c(1) - hx + dx*ix, y, c(3) - hz + dz*iz], &
                                      [c(1) - hx + dx*ix, y, c(3) - hz + dz*(iz + 1)], &
                                      [c(1) - hx + dx*(ix + 1), y, c(3) - hz + dz*(iz + 1)])
                    end if
                end do
            end do
        end do

        call init_mesh(mesh, v0, v1, v2)
    end subroutine make_box

    !> 円柱側面を分割生成し、必要に応じて上下面キャップを追加したメッシュを生成する。
  !! @param[out] mesh 生成した円柱メッシュ。
  !! @param[in] radius 円柱半径 [m]（省略時 0.5）。
  !! @param[in] height 円柱高さ [m]（省略時 1.0）。
  !! @param[in] n_theta 周方向分割数（省略時 24）。
  !! @param[in] n_z 軸方向分割数（省略時 1）。
  !! @param[in] cap 上下面キャップをまとめて指定する後方互換フラグ（省略時 `.true.`）。
  !! @param[in] center 円柱中心座標 `(x,y,z)` [m]（省略時原点）。
  !! @param[in] cap_top 上面キャップを生成するか（省略時 `.true.`）。
  !! @param[in] cap_bottom 下面キャップを生成するか（省略時 `.true.`）。
    subroutine make_cylinder(mesh, radius, height, n_theta, n_z, cap, center, cap_top, cap_bottom)
        type(mesh_type), intent(out) :: mesh
        real(dp), intent(in), optional :: radius, height, center(3)
        integer(i32), intent(in), optional :: n_theta, n_z
        logical, intent(in), optional :: cap, cap_top, cap_bottom
        real(dp) :: r, h, c(3), t0, t1, z0, z1, pi
        integer(i32) :: nt, nz0, iz, it, nelem, itri
        logical :: cap_top0, cap_bottom0
        real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)

        pi = acos(-1.0d0)
        r = 0.5d0; h = 1.0d0; nt = 24; nz0 = 1; cap_top0 = .true.; cap_bottom0 = .true.; c = 0.0d0
        if (present(radius)) r = radius
        if (present(height)) h = height
        if (present(n_theta)) nt = n_theta
        if (present(n_z)) nz0 = n_z
        if (present(cap)) then
            cap_top0 = cap
            cap_bottom0 = cap
        end if
        if (present(cap_top)) cap_top0 = cap_top
        if (present(cap_bottom)) cap_bottom0 = cap_bottom
        if (present(center)) c = center
        if (nt < 3 .or. nz0 <= 0) error stop "n_theta >= 3 and n_z > 0"

        nelem = 2*nt*nz0
        if (cap_bottom0) nelem = nelem + nt
        if (cap_top0) nelem = nelem + nt
        allocate (v0(3, nelem), v1(3, nelem), v2(3, nelem)); itri = 0

        do iz = 0, nz0 - 1
            z0 = c(3) - 0.5d0*h + h*real(iz, dp)/real(nz0, dp)
            z1 = c(3) - 0.5d0*h + h*real(iz + 1, dp)/real(nz0, dp)
            do it = 0, nt - 1
                t0 = 2.0d0*pi*real(it, dp)/real(nt, dp)
                t1 = 2.0d0*pi*real(mod(it + 1, nt), dp)/real(nt, dp)
                call push_tri(v0, v1, v2, itri, &
                              [c(1) + r*cos(t0), c(2) + r*sin(t0), z0], &
                              [c(1) + r*cos(t1), c(2) + r*sin(t1), z0], &
                              [c(1) + r*cos(t1), c(2) + r*sin(t1), z1])
                call push_tri(v0, v1, v2, itri, &
                              [c(1) + r*cos(t0), c(2) + r*sin(t0), z0], &
                              [c(1) + r*cos(t1), c(2) + r*sin(t1), z1], &
                              [c(1) + r*cos(t0), c(2) + r*sin(t0), z1])
            end do
        end do

        if (cap_bottom0) then
            do it = 0, nt - 1
                t0 = 2.0d0*pi*real(it, dp)/real(nt, dp)
                t1 = 2.0d0*pi*real(mod(it + 1, nt), dp)/real(nt, dp)
                call push_tri(v0, v1, v2, itri, &
                              [c(1), c(2), c(3) - 0.5d0*h], &
                              [c(1) + r*cos(t1), c(2) + r*sin(t1), c(3) - 0.5d0*h], &
                              [c(1) + r*cos(t0), c(2) + r*sin(t0), c(3) - 0.5d0*h])
            end do
        end if
        if (cap_top0) then
            do it = 0, nt - 1
                t0 = 2.0d0*pi*real(it, dp)/real(nt, dp)
                t1 = 2.0d0*pi*real(mod(it + 1, nt), dp)/real(nt, dp)
                call push_tri(v0, v1, v2, itri, &
                              [c(1), c(2), c(3) + 0.5d0*h], &
                              [c(1) + r*cos(t0), c(2) + r*sin(t0), c(3) + 0.5d0*h], &
                              [c(1) + r*cos(t1), c(2) + r*sin(t1), c(3) + 0.5d0*h])
            end do
        end if

        call init_mesh(mesh, v0, v1, v2)
    end subroutine make_cylinder

    !> 経度・緯度分割に基づき球面三角形メッシュを生成する。
  !! @param[out] mesh 生成した球面三角形メッシュ。
  !! @param[in] radius 球半径 [m]（省略時 0.5）。
  !! @param[in] n_lon 経度方向分割数（省略時 24）。
  !! @param[in] n_lat 緯度方向分割数（省略時 12）。
  !! @param[in] center 球中心座標 `(x,y,z)` [m]（省略時原点）。
    subroutine make_sphere(mesh, radius, n_lon, n_lat, center)
        type(mesh_type), intent(out) :: mesh
        real(dp), intent(in), optional :: radius, center(3)
        integer(i32), intent(in), optional :: n_lon, n_lat
        real(dp) :: r, c(3), t0, t1, p0, p1, pi
        integer(i32) :: nl, nphi, ilon, ilat, nelem, itri
        real(dp) :: a(3), b(3), c0(3), d(3)
        real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)

        pi = acos(-1.0d0)
        r = 0.5d0; nl = 24; nphi = 12; c = 0.0d0
        if (present(radius)) r = radius
        if (present(n_lon)) nl = n_lon
        if (present(n_lat)) nphi = n_lat
        if (present(center)) c = center
        if (nl < 3 .or. nphi < 2) error stop "n_lon >= 3 and n_lat >= 2"

        nelem = 2*nl*(nphi - 1)
        allocate (v0(3, nelem), v1(3, nelem), v2(3, nelem)); itri = 0

        do ilat = 0, nphi - 1
            p0 = pi*real(ilat, dp)/real(nphi, dp)
            p1 = pi*real(ilat + 1, dp)/real(nphi, dp)
            do ilon = 0, nl - 1
                t0 = 2.0d0*pi*real(ilon, dp)/real(nl, dp)
                t1 = 2.0d0*pi*real(mod(ilon + 1, nl), dp)/real(nl, dp)
                call sph(r, c, t0, p0, a)
                call sph(r, c, t1, p0, b)
                call sph(r, c, t0, p1, c0)
                call sph(r, c, t1, p1, d)
                if (ilat == 0) then
                    call push_tri(v0, v1, v2, itri, a, c0, d)
                else if (ilat == nphi - 1) then
                    call push_tri(v0, v1, v2, itri, a, d, b)
                else
                    call push_tri(v0, v1, v2, itri, a, d, b)
                    call push_tri(v0, v1, v2, itri, a, c0, d)
                end if
            end do
        end do

        call init_mesh(mesh, v0, v1, v2)
    end subroutine make_sphere

    !> XY平面上で、中心からのレイと長方形境界の最短交点を返す。
    subroutine ray_to_rectangle(hx, hy, dx, dy, x_min, x_max, y_min, y_max, z_plane, p, edge_id)
        real(dp), intent(in) :: hx, hy, dx, dy, x_min, x_max, y_min, y_max, z_plane
        real(dp), intent(out) :: p(3)
        integer(i32), intent(out) :: edge_id
        real(dp) :: s_best, s, x, y
        real(dp), parameter :: eps = 1.0d-12

        s_best = huge(1.0d0)
        edge_id = 0_i32
        p = 0.0d0

        if (dx > eps) then
            s = (x_max - hx)/dx
            y = hy + s*dy
            if (s > eps .and. y >= y_min - eps .and. y <= y_max + eps .and. s < s_best) then
                s_best = s
                edge_id = 2_i32
                p = [x_max, y, z_plane]
            end if
        else if (dx < -eps) then
            s = (x_min - hx)/dx
            y = hy + s*dy
            if (s > eps .and. y >= y_min - eps .and. y <= y_max + eps .and. s < s_best) then
                s_best = s
                edge_id = 1_i32
                p = [x_min, y, z_plane]
            end if
        end if
        if (dy > eps) then
            s = (y_max - hy)/dy
            x = hx + s*dx
            if (s > eps .and. x >= x_min - eps .and. x <= x_max + eps .and. s < s_best) then
                s_best = s
                edge_id = 4_i32
                p = [x, y_max, z_plane]
            end if
        else if (dy < -eps) then
            s = (y_min - hy)/dy
            x = hx + s*dx
            if (s > eps .and. x >= x_min - eps .and. x <= x_max + eps .and. s < s_best) then
                edge_id = 3_i32
                p = [x, y_min, z_plane]
            end if
        end if

        if (edge_id == 0_i32) error stop "failed to intersect ray with rectangle boundary"
    end subroutine ray_to_rectangle

    !> 長方形外周を反時計回りに見たとき、辺遷移で通過するコーナー数を返す。
    integer(i32) function transition_corner_count(edge_from, edge_to) result(n_corner)
        integer(i32), intent(in) :: edge_from, edge_to
        integer(i32) :: i_from, i_to

        i_from = edge_order_index(edge_from)
        i_to = edge_order_index(edge_to)
        n_corner = modulo(i_to - i_from, 4_i32)
    end function transition_corner_count

    !> 辺遷移時に通過するコーナー列（最大3点）を返す。
    subroutine transition_corners(edge_from, edge_to, x_min, x_max, y_min, y_max, z_plane, corners, n_corner)
        integer(i32), intent(in) :: edge_from, edge_to
        real(dp), intent(in) :: x_min, x_max, y_min, y_max, z_plane
        real(dp), intent(out) :: corners(3, 3)
        integer(i32), intent(out) :: n_corner
        integer(i32) :: k, edge_curr, edge_next

        corners = 0.0d0
        n_corner = transition_corner_count(edge_from, edge_to)
        edge_curr = edge_from
        do k = 1, n_corner
            edge_next = edge_next_ccw(edge_curr)
            call corner_between(edge_curr, edge_next, x_min, x_max, y_min, y_max, z_plane, corners(:, k))
            edge_curr = edge_next
        end do
    end subroutine transition_corners

    !> 境界辺の反時計回り順序インデックスを返す。
    integer(i32) function edge_order_index(edge_id) result(idx)
        integer(i32), intent(in) :: edge_id

        select case (edge_id)
        case (2_i32)
            idx = 0_i32
        case (4_i32)
            idx = 1_i32
        case (1_i32)
            idx = 2_i32
        case (3_i32)
            idx = 3_i32
        case default
            error stop "unknown edge id"
        end select
    end function edge_order_index

    !> 境界辺の次（反時計回り）を返す。
    integer(i32) function edge_next_ccw(edge_id) result(next_id)
        integer(i32), intent(in) :: edge_id

        select case (edge_id)
        case (2_i32)
            next_id = 4_i32
        case (4_i32)
            next_id = 1_i32
        case (1_i32)
            next_id = 3_i32
        case (3_i32)
            next_id = 2_i32
        case default
            error stop "unknown edge id"
        end select
    end function edge_next_ccw

    !> 反時計回りに隣接する辺ペアに対し、対応する長方形コーナー座標を返す。
    subroutine corner_between(edge_from, edge_to, x_min, x_max, y_min, y_max, z_plane, corner)
        integer(i32), intent(in) :: edge_from, edge_to
        real(dp), intent(in) :: x_min, x_max, y_min, y_max, z_plane
        real(dp), intent(out) :: corner(3)

        select case (edge_from)
        case (2_i32)
            if (edge_to /= 4_i32) error stop "invalid edge transition"
            corner = [x_max, y_max, z_plane]
        case (4_i32)
            if (edge_to /= 1_i32) error stop "invalid edge transition"
            corner = [x_min, y_max, z_plane]
        case (1_i32)
            if (edge_to /= 3_i32) error stop "invalid edge transition"
            corner = [x_min, y_min, z_plane]
        case (3_i32)
            if (edge_to /= 2_i32) error stop "invalid edge transition"
            corner = [x_max, y_min, z_plane]
        case default
            error stop "unknown edge id"
        end select
    end subroutine corner_between

    !> 球座標 `(theta, phi)` を中心 `c`・半径 `r` の直交座標へ変換する。
  !! @param[in] r 球半径 [m]。
  !! @param[in] c 球中心座標 `(x,y,z)` [m]。
  !! @param[in] theta 方位角 [rad]。
  !! @param[in] phi 極角（+Z軸基準） [rad]。
  !! @param[out] p 変換後の直交座標 `(x,y,z)` [m]。
    pure subroutine sph(r, c, theta, phi, p)
        real(dp), intent(in) :: r, c(3), theta, phi
        real(dp), intent(out) :: p(3)
        p(1) = c(1) + r*sin(phi)*cos(theta)
        p(2) = c(2) + r*sin(phi)*sin(theta)
        p(3) = c(3) + r*cos(phi)
    end subroutine sph

    !> 三角形頂点 `a,b,c` を出力配列の次インデックスへ書き込む。
  !! @param[inout] v0 三角形頂点0を保持する配列 `v0(3,nelem)`。
  !! @param[inout] v1 三角形頂点1を保持する配列 `v1(3,nelem)`。
  !! @param[inout] v2 三角形頂点2を保持する配列 `v2(3,nelem)`。
  !! @param[inout] itri 現在までに書き込んだ三角形数（呼び出し内で1増加）。
  !! @param[in] a 追加する三角形の頂点A座標。
  !! @param[in] b 追加する三角形の頂点B座標。
  !! @param[in] c 追加する三角形の頂点C座標。
    pure subroutine push_tri(v0, v1, v2, itri, a, b, c)
        real(dp), intent(inout) :: v0(:, :), v1(:, :), v2(:, :)
        integer(i32), intent(inout) :: itri
        real(dp), intent(in) :: a(3), b(3), c(3)
        itri = itri + 1
        v0(:, itri) = a
        v1(:, itri) = b
        v2(:, itri) = c
    end subroutine push_tri

end module bem_templates
