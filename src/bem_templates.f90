!> 平面/箱/円柱/球テンプレートから三角形メッシュを生成するユーティリティ。
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

    !> 円柱側面を分割生成し、必要に応じて上下キャップを追加したメッシュを生成する。
  !! @param[out] mesh 生成した円柱メッシュ。
  !! @param[in] radius 円柱半径 [m]（省略時 0.5）。
  !! @param[in] height 円柱高さ [m]（省略時 1.0）。
  !! @param[in] n_theta 周方向分割数（省略時 24）。
  !! @param[in] n_z 軸方向分割数（省略時 1）。
  !! @param[in] cap 上下キャップを生成するか（省略時 `.true.`）。
  !! @param[in] center 円柱中心座標 `(x,y,z)` [m]（省略時原点）。
    subroutine make_cylinder(mesh, radius, height, n_theta, n_z, cap, center)
        type(mesh_type), intent(out) :: mesh
        real(dp), intent(in), optional :: radius, height, center(3)
        integer(i32), intent(in), optional :: n_theta, n_z
        logical, intent(in), optional :: cap
        real(dp) :: r, h, c(3), t0, t1, z0, z1, pi
        integer(i32) :: nt, nz0, iz, it, nelem, itri
        logical :: cap0
        real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)

        pi = acos(-1.0d0)
        r = 0.5d0; h = 1.0d0; nt = 24; nz0 = 1; cap0 = .true.; c = 0.0d0
        if (present(radius)) r = radius
        if (present(height)) h = height
        if (present(n_theta)) nt = n_theta
        if (present(n_z)) nz0 = n_z
        if (present(cap)) cap0 = cap
        if (present(center)) c = center
        if (nt < 3 .or. nz0 <= 0) error stop "n_theta >= 3 and n_z > 0"

        nelem = 2*nt*nz0
        if (cap0) nelem = nelem + 2*nt
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

        if (cap0) then
            do it = 0, nt - 1
                t0 = 2.0d0*pi*real(it, dp)/real(nt, dp)
                t1 = 2.0d0*pi*real(mod(it + 1, nt), dp)/real(nt, dp)
                call push_tri(v0, v1, v2, itri, &
                              [c(1), c(2), c(3) - 0.5d0*h], &
                              [c(1) + r*cos(t1), c(2) + r*sin(t1), c(3) - 0.5d0*h], &
                              [c(1) + r*cos(t0), c(2) + r*sin(t0), c(3) - 0.5d0*h])
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
