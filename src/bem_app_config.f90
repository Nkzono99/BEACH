!> 実行設定(TOML)を既定値付きで保持し、メッシュ/粒子初期化へ変換する構成管理モジュール。
module bem_app_config
    use bem_kinds, only: dp, i32
    use bem_types, only: sim_config, mesh_type, particles_soa, bc_open, bc_reflect, bc_periodic
    use bem_templates, only: make_plane, make_box, make_cylinder, make_sphere
    use bem_mesh, only: init_mesh
    use bem_importers, only: load_obj_mesh
    use bem_injection, only: seed_rng, init_random_beam_particles, &
                             sample_uniform_positions, sample_shifted_maxwell_velocities
    use bem_particles, only: init_particles
    implicit none

    integer, parameter :: max_templates = 8
    integer, parameter :: max_particle_species = 8

    type :: particle_species_spec
        logical :: enabled = .false.
        integer(i32) :: n_particles = 0_i32
        real(dp) :: q_particle = -1.602176634d-19
        real(dp) :: m_particle = 9.10938356d-31
        real(dp) :: w_particle = 1.0d0
        real(dp) :: pos_low(3) = [-0.4d0, -0.4d0, 0.2d0]
        real(dp) :: pos_high(3) = [0.4d0, 0.4d0, 0.5d0]
        real(dp) :: drift_velocity(3) = [0.0d0, 0.0d0, -8.0d5]
        real(dp) :: temperature_k = 2.0d4
    end type particle_species_spec

    !> 1つのテンプレート形状(plane/box/cylinder/sphere)の有効化フラグと幾何パラメータを保持する設定型。
    type :: template_spec
        logical :: enabled = .false.
        character(len=16) :: kind = 'plane'
        real(dp) :: center(3) = 0.0d0
        real(dp) :: size_x = 1.0d0
        real(dp) :: size_y = 1.0d0
        real(dp) :: size(3) = [1.0d0, 1.0d0, 1.0d0]
        integer(i32) :: nx = 1
        integer(i32) :: ny = 1
        integer(i32) :: nz = 1
        real(dp) :: radius = 0.5d0
        real(dp) :: height = 1.0d0
        integer(i32) :: n_theta = 24
        integer(i32) :: n_z = 1
        logical :: cap = .true.
        integer(i32) :: n_lon = 24
        integer(i32) :: n_lat = 12
    end type template_spec

    !> シミュレーション条件・メッシュ入力・粒子注入条件・出力先を一元管理するアプリ設定型。
    type :: app_config
        character(len=16) :: mesh_mode = 'auto'  ! auto / obj / template
        character(len=256) :: obj_path = 'examples/simple_plate.obj'
        integer(i32) :: n_templates = 1
        type(template_spec) :: templates(max_templates)

        integer(i32) :: rng_seed = 12345_i32
        integer(i32) :: n_particles = 256_i32
        real(dp) :: q_particle = -1.602176634d-19
        real(dp) :: m_particle = 9.10938356d-31
        real(dp) :: w_particle = 1.0d0
        real(dp) :: pos_low(3) = [-0.4d0, -0.4d0, 0.2d0]
        real(dp) :: pos_high(3) = [0.4d0, 0.4d0, 0.5d0]
        real(dp) :: drift_velocity(3) = [0.0d0, 0.0d0, -8.0d5]
        real(dp) :: temperature_k = 2.0d4
        integer(i32) :: n_particle_species = 0_i32
        type(particle_species_spec) :: particle_species(max_particle_species)

        logical :: write_output = .true.
        character(len=256) :: output_dir = 'outputs/latest'
        integer(i32) :: history_stride = 1

        type(sim_config) :: sim
    end type app_config

contains

    !> 設定から総粒子数を返す。
  !! @param[in] cfg 入力引数。
  !! @return 総粒子数。
    integer(i32) function total_particles_from_config(cfg) result(total_n)
        type(app_config), intent(in) :: cfg
        integer(i32) :: s

        if (cfg%n_particle_species <= 0) then
            total_n = max(0_i32, cfg%n_particles)
            return
        end if

        total_n = 0_i32
        do s = 1, cfg%n_particle_species
            if (cfg%particle_species(s)%enabled) total_n = total_n + max(0_i32, cfg%particle_species(s)%n_particles)
        end do
    end function total_particles_from_config

    !> `app_config` に実用的な既定値を設定し、設定ファイル未指定でも実行可能な状態を作る。
  !! @param[out] cfg 出力引数。
    subroutine default_app_config(cfg)
        type(app_config), intent(out) :: cfg
        cfg%sim%dt = 1.0d-9
        cfg%sim%npcls_per_step = 64
        cfg%sim%max_step = 400
        cfg%sim%tol_rel = 1.0d-8
        cfg%sim%softening = 1.0d-6
        cfg%sim%b0 = [0.0d0, 0.0d0, 0.0d0]

        cfg%templates(1)%enabled = .true.
        cfg%templates(1)%kind = 'plane'
        cfg%templates(1)%size_x = 1.0d0
        cfg%templates(1)%size_y = 1.0d0
        cfg%templates(1)%nx = 1
        cfg%templates(1)%ny = 1
        cfg%templates(1)%center = [0.0d0, 0.0d0, 0.0d0]
    end subroutine default_app_config

    !> TOML設定ファイルを読み込み、既定値に対して上書き適用する。
  !! @param[in] path 入力引数。
  !! @param[inout] cfg 入出力引数。
    subroutine load_app_config(path, cfg)
        character(len=*), intent(in) :: path
        type(app_config), intent(inout) :: cfg
        if (.not. ends_with(lower(trim(path)), '.toml')) then
            error stop 'Only TOML config is supported. Please pass a .toml file.'
        end if
        call load_toml_config(path, cfg)
    end subroutine load_app_config

    !> `mesh_mode` とファイル存在有無に応じて OBJ 読込またはテンプレート生成を選択してメッシュを構築する。
  !! @param[in] cfg 入力引数。
  !! @param[out] mesh 出力引数。
    subroutine build_mesh_from_config(cfg, mesh)
        type(app_config), intent(in) :: cfg
        type(mesh_type), intent(out) :: mesh
        logical :: has_obj

        select case (trim(lower(cfg%mesh_mode)))
        case ('obj')
            call load_obj_mesh(trim(cfg%obj_path), mesh)
        case ('template')
            call build_template_mesh(cfg, mesh)
        case default
            inquire (file=trim(cfg%obj_path), exist=has_obj)
            if (has_obj) then
                call load_obj_mesh(trim(cfg%obj_path), mesh)
            else
                call build_template_mesh(cfg, mesh)
            end if
        end select
    end subroutine build_mesh_from_config

    !> 設定値から乱数シードとビーム注入条件を適用し、粒子SoAを初期化する。
  !! @param[in] cfg 入力引数。
  !! @param[out] pcls 出力引数。
    subroutine init_particles_from_config(cfg, pcls)
        type(app_config), intent(in) :: cfg
        type(particles_soa), intent(out) :: pcls

        integer(i32) :: s, total_n, max_n, out_idx, rank_idx
        integer(i32), allocatable :: counts(:)
        real(dp), allocatable :: x_species(:, :, :), v_species(:, :, :)
        real(dp), allocatable :: q(:), m(:), w(:), x(:, :), v(:, :)

        if (cfg%n_particle_species <= 0) then
            call seed_rng([cfg%rng_seed])
            call init_random_beam_particles( &
                pcls=pcls, &
                n=cfg%n_particles, &
                q_particle=cfg%q_particle, &
                m_particle=cfg%m_particle, &
                w_particle=cfg%w_particle, &
                pos_low=cfg%pos_low, &
                pos_high=cfg%pos_high, &
                drift_velocity=cfg%drift_velocity, &
                temperature_k=cfg%temperature_k &
                )
            return
        end if

        call seed_rng([cfg%rng_seed])

        allocate (counts(cfg%n_particle_species))
        counts = 0_i32
        do s = 1, cfg%n_particle_species
            if (cfg%particle_species(s)%enabled) then
                counts(s) = max(0_i32, cfg%particle_species(s)%n_particles)
            else
                counts(s) = 0_i32
            end if
        end do
        total_n = sum(counts)
        max_n = max(1_i32, maxval(counts))

        allocate (x_species(3, max_n, cfg%n_particle_species))
        allocate (v_species(3, max_n, cfg%n_particle_species))
        x_species = 0.0d0
        v_species = 0.0d0

        do s = 1, cfg%n_particle_species
            if (counts(s) <= 0) cycle
            call sample_species_state( &
                cfg%particle_species(s), counts(s), x_species(:, 1:counts(s), s), v_species(:, 1:counts(s), s) &
                )
        end do

        allocate (x(3, total_n), v(3, total_n), q(total_n), m(total_n), w(total_n))
        out_idx = 0
        do rank_idx = 1, maxval(counts)
            do s = 1, cfg%n_particle_species
                if (rank_idx > counts(s)) cycle
                out_idx = out_idx + 1
                x(:, out_idx) = x_species(:, rank_idx, s)
                v(:, out_idx) = v_species(:, rank_idx, s)
                q(out_idx) = cfg%particle_species(s)%q_particle
                m(out_idx) = cfg%particle_species(s)%m_particle
                w(out_idx) = cfg%particle_species(s)%w_particle
            end do
        end do

        call init_particles(pcls, x, v, q, m, w)
    end subroutine init_particles_from_config

    !> 粒子注入用乱数のシードを設定する。
  !! @param[in] cfg 入力引数。
    subroutine seed_particles_from_config(cfg)
        type(app_config), intent(in) :: cfg
        call seed_rng([cfg%rng_seed])
    end subroutine seed_particles_from_config

    !> 指定区間（1始まり）の粒子インデックスに対応する粒子バッチを生成する。
  !! @param[in] cfg 入力引数。
  !! @param[in] start_idx 生成開始インデックス（1始まり）。
  !! @param[in] batch_n 生成粒子数。
  !! @param[out] pcls 出力引数。
    subroutine init_particle_batch_from_config(cfg, start_idx, batch_n, pcls)
        type(app_config), intent(in) :: cfg
        integer(i32), intent(in) :: start_idx, batch_n
        type(particles_soa), intent(out) :: pcls

        integer(i32) :: s, i, global_idx, end_idx, total_n, max_rank, out_idx
        integer(i32), allocatable :: counts(:), batch_counts(:), species_cursor(:), species_id(:)
        real(dp), allocatable :: x_species(:, :, :), v_species(:, :, :), x(:, :), v(:, :), q(:), m(:), w(:)

        if (batch_n < 0_i32) error stop 'batch_n must be >= 0.'
        if (batch_n == 0_i32) then
            allocate (x(3, 0), v(3, 0), q(0), m(0), w(0))
            call init_particles(pcls, x, v, q, m, w)
            return
        end if

        if (cfg%n_particle_species <= 0) then
            if (start_idx < 1_i32 .or. start_idx + batch_n - 1_i32 > cfg%n_particles) then
                error stop 'Requested particle batch is out of range.'
            end if
            call init_random_beam_particles( &
                pcls=pcls, &
                n=batch_n, &
                q_particle=cfg%q_particle, &
                m_particle=cfg%m_particle, &
                w_particle=cfg%w_particle, &
                pos_low=cfg%pos_low, &
                pos_high=cfg%pos_high, &
                drift_velocity=cfg%drift_velocity, &
                temperature_k=cfg%temperature_k &
                )
            return
        end if

        allocate (counts(cfg%n_particle_species))
        counts = 0_i32
        do s = 1, cfg%n_particle_species
            if (cfg%particle_species(s)%enabled) counts(s) = max(0_i32, cfg%particle_species(s)%n_particles)
        end do
        total_n = sum(counts)
        if (start_idx < 1_i32 .or. start_idx + batch_n - 1_i32 > total_n) then
            error stop 'Requested particle batch is out of range.'
        end if

        allocate (species_id(batch_n), batch_counts(cfg%n_particle_species))
        batch_counts = 0_i32
        end_idx = start_idx + batch_n - 1_i32
        global_idx = 0_i32
        max_rank = max(1_i32, maxval(counts))
        out_idx = 0_i32
        do i = 1, max_rank
            do s = 1, cfg%n_particle_species
                if (i > counts(s)) cycle
                global_idx = global_idx + 1_i32
                if (global_idx < start_idx) cycle
                if (global_idx > end_idx) exit
                out_idx = out_idx + 1_i32
                species_id(out_idx) = s
                batch_counts(s) = batch_counts(s) + 1_i32
            end do
            if (global_idx > end_idx) exit
        end do

        allocate (x_species(3, max(1_i32, maxval(batch_counts)), cfg%n_particle_species))
        allocate (v_species(3, max(1_i32, maxval(batch_counts)), cfg%n_particle_species))
        x_species = 0.0d0
        v_species = 0.0d0
        do s = 1, cfg%n_particle_species
            if (batch_counts(s) <= 0_i32) cycle
            call sample_species_state( &
                cfg%particle_species(s), batch_counts(s), &
                x_species(:, 1:batch_counts(s), s), &
                v_species(:, 1:batch_counts(s), s) &
            )
        end do

        allocate (x(3, batch_n), v(3, batch_n), q(batch_n), m(batch_n), w(batch_n), species_cursor(cfg%n_particle_species))
        species_cursor = 0_i32
        do i = 1, batch_n
            s = species_id(i)
            species_cursor(s) = species_cursor(s) + 1_i32
            x(:, i) = x_species(:, species_cursor(s), s)
            v(:, i) = v_species(:, species_cursor(s), s)
            q(i) = cfg%particle_species(s)%q_particle
            m(i) = cfg%particle_species(s)%m_particle
            w(i) = cfg%particle_species(s)%w_particle
        end do

        call init_particles(pcls, x, v, q, m, w)
    end subroutine init_particle_batch_from_config

    subroutine sample_species_state(spec, n, x, v)
        type(particle_species_spec), intent(in) :: spec
        integer(i32), intent(in) :: n
        real(dp), intent(out) :: x(:, :), v(:, :)

        if (n <= 0) return
        call sample_uniform_positions(spec%pos_low, spec%pos_high, x)
        call sample_shifted_maxwell_velocities( &
            spec%drift_velocity, spec%m_particle, v, temperature_k=spec%temperature_k &
            )
    end subroutine sample_species_state

    !> 有効化された複数テンプレートを結合し、1つの三角形メッシュとして初期化する。
  !! @param[in] cfg 入力引数。
  !! @param[out] mesh 出力引数。
    subroutine build_template_mesh(cfg, mesh)
        type(app_config), intent(in) :: cfg
        type(mesh_type), intent(out) :: mesh
        type(mesh_type) :: part
        real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)
        integer(i32) :: i

        allocate (v0(3, 0), v1(3, 0), v2(3, 0))
        do i = 1, min(cfg%n_templates, int(max_templates, i32))
            if (.not. cfg%templates(i)%enabled) cycle
            call build_one_template(cfg%templates(i), part)
            call append_triangles(v0, v1, v2, part%v0, part%v1, part%v2)
        end do

        if (size(v0, 2) == 0) then
            error stop 'No enabled template found in configuration.'
        end if
        call init_mesh(mesh, v0, v1, v2)
    end subroutine build_template_mesh

    !> 単一テンプレート設定を形状別生成ルーチンへディスパッチしてメッシュ化する。
  !! @param[in] spec 入力引数。
  !! @param[out] mesh 出力引数。
    subroutine build_one_template(spec, mesh)
        type(template_spec), intent(in) :: spec
        type(mesh_type), intent(out) :: mesh

        select case (trim(lower(spec%kind)))
        case ('plane')
            call make_plane(mesh, size_x=spec%size_x, size_y=spec%size_y, nx=spec%nx, ny=spec%ny, center=spec%center)
        case ('box')
            call make_box(mesh, size=spec%size, center=spec%center, nx=spec%nx, ny=spec%ny, nz=spec%nz)
        case ('cylinder')
            call make_cylinder( &
                mesh, radius=spec%radius, height=spec%height, &
                n_theta=spec%n_theta, n_z=spec%n_z, &
                cap=spec%cap, center=spec%center &
            )
        case ('sphere')
            call make_sphere(mesh, radius=spec%radius, n_lon=spec%n_lon, n_lat=spec%n_lat, center=spec%center)
        case default
            error stop 'Unknown template kind in config.'
        end select
    end subroutine build_one_template

    !> 既存三角形配列へ追加三角形群を連結し、結合済み配列へ再確保して詰め替える。
  !! @param[inout] v0 入出力引数。
  !! @param[inout] v1 入出力引数。
  !! @param[inout] v2 入出力引数。
  !! @param[in] add_v0 入力引数。
  !! @param[in] add_v1 入力引数。
  !! @param[in] add_v2 入力引数。
    subroutine append_triangles(v0, v1, v2, add_v0, add_v1, add_v2)
        real(dp), allocatable, intent(inout) :: v0(:, :), v1(:, :), v2(:, :)
        real(dp), intent(in) :: add_v0(:, :), add_v1(:, :), add_v2(:, :)
        real(dp), allocatable :: t0(:, :), t1(:, :), t2(:, :)
        integer :: n0, n1

        n0 = size(v0, 2)
        n1 = size(add_v0, 2)
        allocate (t0(3, n0 + n1), t1(3, n0 + n1), t2(3, n0 + n1))
        if (n0 > 0) then
            t0(:, 1:n0) = v0
            t1(:, 1:n0) = v1
            t2(:, 1:n0) = v2
        end if
        t0(:, n0 + 1:n0 + n1) = add_v0
        t1(:, n0 + 1:n0 + n1) = add_v1
        t2(:, n0 + 1:n0 + n1) = add_v2
        call move_alloc(t0, v0)
        call move_alloc(t1, v1)
        call move_alloc(t2, v2)
    end subroutine append_triangles

    !> 最小TOMLパーサでセクション/キー値を走査し、`app_config` 各項目へ反映する。
  !! @param[in] path 入力引数。
  !! @param[inout] cfg 入出力引数。
    subroutine load_toml_config(path, cfg)
        character(len=*), intent(in) :: path
        type(app_config), intent(inout) :: cfg
        integer :: u, ios, i, t_idx, s_idx
        character(len=512) :: raw, line, section

        t_idx = 0
        section = ''
        s_idx = 0
        open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'Could not open TOML file.'

        do
            read (u, '(A)', iostat=ios) raw
            if (ios /= 0) exit
            line = strip_comment(trim(raw))
            if (len_trim(line) == 0) cycle

            if (line(1:1) == '[') then
                if (trim(line) == '[[mesh.templates]]') then
                    t_idx = t_idx + 1
                    if (t_idx > max_templates) error stop 'Too many mesh.templates entries.'
                    cfg%templates(t_idx)%enabled = .true.
                    section = 'mesh.template'
                else if (trim(line) == '[[particles.species]]') then
                    s_idx = s_idx + 1
                    if (s_idx > max_particle_species) error stop 'Too many particles.species entries.'
                    cfg%particle_species(s_idx) = species_from_legacy(cfg)
                    cfg%particle_species(s_idx)%enabled = .true.
                    section = 'particles.species'
                else
                    section = lower(trim(adjustl(line(2:len_trim(line) - 1))))
                end if
                cycle
            end if

            select case (trim(section))
            case ('sim')
                call apply_sim_kv(cfg, line)
            case ('particles')
                call apply_particles_kv(cfg, line)
            case ('particles.species')
                call apply_particles_species_kv(cfg%particle_species(s_idx), line)
            case ('mesh')
                call apply_mesh_kv(cfg, line)
            case ('mesh.template')
                call apply_template_kv(cfg%templates(t_idx), line)
            case ('output')
                call apply_output_kv(cfg, line)
            end select
        end do
        close (u)
        if (t_idx > 0) cfg%n_templates = t_idx
        if (s_idx > 0) then
            cfg%n_particle_species = s_idx
            cfg%n_particles = sum(max(0_i32, cfg%particle_species(1:s_idx)%n_particles))
        end if
    end subroutine load_toml_config

    !> `[sim]` セクションのキー値を `sim_config` へ変換して適用する。
  !! @param[inout] cfg 入出力引数。
  !! @param[in] line 入力引数。
    subroutine apply_sim_kv(cfg, line)
        type(app_config), intent(inout) :: cfg
        character(len=*), intent(in) :: line
        character(len=64) :: k
        character(len=256) :: v
        call split_key_value(line, k, v)
        select case (trim(k))
        case ('dt'); call parse_real(v, cfg%sim%dt)
        case ('npcls_per_step'); call parse_int(v, cfg%sim%npcls_per_step)
        case ('max_step'); call parse_int(v, cfg%sim%max_step)
        case ('tol_rel'); call parse_real(v, cfg%sim%tol_rel)
        case ('q_floor'); call parse_real(v, cfg%sim%q_floor)
        case ('softening'); call parse_real(v, cfg%sim%softening)
        case ('use_hybrid'); call parse_logical(v, cfg%sim%use_hybrid)
        case ('r_switch_factor'); call parse_real(v, cfg%sim%r_switch_factor)
        case ('n_sub'); call parse_int(v, cfg%sim%n_sub)
        case ('softening_factor'); call parse_real(v, cfg%sim%softening_factor)
        case ('b0'); call parse_real3(v, cfg%sim%b0)
        case ('use_box'); call parse_logical(v, cfg%sim%use_box)
        case ('box_min'); call parse_real3(v, cfg%sim%box_min)
        case ('box_max'); call parse_real3(v, cfg%sim%box_max)
        case ('bc_x_low'); call parse_boundary_mode(v, cfg%sim%bc_low(1))
        case ('bc_x_high'); call parse_boundary_mode(v, cfg%sim%bc_high(1))
        case ('bc_y_low'); call parse_boundary_mode(v, cfg%sim%bc_low(2))
        case ('bc_y_high'); call parse_boundary_mode(v, cfg%sim%bc_high(2))
        case ('bc_z_low'); call parse_boundary_mode(v, cfg%sim%bc_low(3))
        case ('bc_z_high'); call parse_boundary_mode(v, cfg%sim%bc_high(3))
        end select
    end subroutine apply_sim_kv

    !> `[particles]` セクションのキー値を粒子注入パラメータへ適用する。
  !! @param[inout] cfg 入出力引数。
  !! @param[in] line 入力引数。
    subroutine apply_particles_kv(cfg, line)
        type(app_config), intent(inout) :: cfg
        character(len=*), intent(in) :: line
        character(len=64) :: k
        character(len=256) :: v
        call split_key_value(line, k, v)
        select case (trim(k))
        case ('rng_seed'); call parse_int(v, cfg%rng_seed)
        case ('n_particles'); call parse_int(v, cfg%n_particles)
        case ('q_particle'); call parse_real(v, cfg%q_particle)
        case ('m_particle'); call parse_real(v, cfg%m_particle)
        case ('w_particle'); call parse_real(v, cfg%w_particle)
        case ('pos_low'); call parse_real3(v, cfg%pos_low)
        case ('pos_high'); call parse_real3(v, cfg%pos_high)
        case ('drift_velocity'); call parse_real3(v, cfg%drift_velocity)
        case ('temperature_k'); call parse_real(v, cfg%temperature_k)
        end select
    end subroutine apply_particles_kv

    subroutine apply_particles_species_kv(spec, line)
        type(particle_species_spec), intent(inout) :: spec
        character(len=*), intent(in) :: line
        character(len=64) :: k
        character(len=256) :: v
        call split_key_value(line, k, v)
        select case (trim(k))
        case ('enabled'); call parse_logical(v, spec%enabled)
        case ('n_particles'); call parse_int(v, spec%n_particles)
        case ('q_particle'); call parse_real(v, spec%q_particle)
        case ('m_particle'); call parse_real(v, spec%m_particle)
        case ('w_particle'); call parse_real(v, spec%w_particle)
        case ('pos_low'); call parse_real3(v, spec%pos_low)
        case ('pos_high'); call parse_real3(v, spec%pos_high)
        case ('drift_velocity'); call parse_real3(v, spec%drift_velocity)
        case ('temperature_k'); call parse_real(v, spec%temperature_k)
        end select
    end subroutine apply_particles_species_kv

    pure function species_from_legacy(cfg) result(spec)
        type(app_config), intent(in) :: cfg
        type(particle_species_spec) :: spec
        spec%enabled = .true.
        spec%n_particles = cfg%n_particles
        spec%q_particle = cfg%q_particle
        spec%m_particle = cfg%m_particle
        spec%w_particle = cfg%w_particle
        spec%pos_low = cfg%pos_low
        spec%pos_high = cfg%pos_high
        spec%drift_velocity = cfg%drift_velocity
        spec%temperature_k = cfg%temperature_k
    end function species_from_legacy

    !> `[mesh]` セクションのキー値をメッシュ入力モード/OBJパスへ適用する。
  !! @param[inout] cfg 入出力引数。
  !! @param[in] line 入力引数。
    subroutine apply_mesh_kv(cfg, line)
        type(app_config), intent(inout) :: cfg
        character(len=*), intent(in) :: line
        character(len=64) :: k
        character(len=256) :: v
        call split_key_value(line, k, v)
        select case (trim(k))
        case ('mode'); call parse_string(v, cfg%mesh_mode)
        case ('obj_path'); call parse_string(v, cfg%obj_path)
        case ('n_templates'); call parse_int(v, cfg%n_templates)
        end select
    end subroutine apply_mesh_kv

    !> `[template]` セクションのキー値を単一テンプレート設定へ適用する。
  !! @param[inout] spec 入出力引数。
  !! @param[in] line 入力引数。
    subroutine apply_template_kv(spec, line)
        type(template_spec), intent(inout) :: spec
        character(len=*), intent(in) :: line
        character(len=64) :: k
        character(len=256) :: v
        call split_key_value(line, k, v)
        select case (trim(k))
        case ('enabled'); call parse_logical(v, spec%enabled)
        case ('kind'); call parse_string(v, spec%kind)
        case ('center'); call parse_real3(v, spec%center)
        case ('size_x'); call parse_real(v, spec%size_x)
        case ('size_y'); call parse_real(v, spec%size_y)
        case ('size'); call parse_real3(v, spec%size)
        case ('nx'); call parse_int(v, spec%nx)
        case ('ny'); call parse_int(v, spec%ny)
        case ('nz'); call parse_int(v, spec%nz)
        case ('radius'); call parse_real(v, spec%radius)
        case ('height'); call parse_real(v, spec%height)
        case ('n_theta'); call parse_int(v, spec%n_theta)
        case ('n_z'); call parse_int(v, spec%n_z)
        case ('cap'); call parse_logical(v, spec%cap)
        case ('n_lon'); call parse_int(v, spec%n_lon)
        case ('n_lat'); call parse_int(v, spec%n_lat)
        end select
    end subroutine apply_template_kv

    !> `[output]` セクションのキー値を出力有効化フラグと出力先ディレクトリへ適用する。
  !! @param[inout] cfg 入出力引数。
  !! @param[in] line 入力引数。
    subroutine apply_output_kv(cfg, line)
        type(app_config), intent(inout) :: cfg
        character(len=*), intent(in) :: line
        character(len=64) :: k
        character(len=256) :: v
        call split_key_value(line, k, v)
        select case (trim(k))
        case ('write_files'); call parse_logical(v, cfg%write_output)
        case ('dir'); call parse_string(v, cfg%output_dir)
        case ('history_stride'); call parse_int(v, cfg%history_stride)
        end select
    end subroutine apply_output_kv

    !> `key = value` 形式の1行を分割し、前後空白を除去して返す。
  !! @param[in] line 入力引数。
  !! @param[out] key 出力引数。
  !! @param[out] value 出力引数。
    subroutine split_key_value(line, key, value)
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: key
        character(len=*), intent(out) :: value
        integer :: p
        p = index(line, '=')
        if (p <= 0) then
            key = ''
            value = ''
            return
        end if
        key = lower(trim(adjustl(line(:p - 1))))
        value = trim(adjustl(line(p + 1:)))
    end subroutine split_key_value

    !> 文字列を倍精度実数として読み取り、設定値へ変換する。
  !! @param[in] text 入力引数。
  !! @param[out] out 出力引数。
    subroutine parse_real(text, out)
        character(len=*), intent(in) :: text
        real(dp), intent(out) :: out
        read (text, *) out
    end subroutine parse_real

    !> 文字列を32bit整数として読み取り、設定値へ変換する。
  !! @param[in] text 入力引数。
  !! @param[out] out 出力引数。
    subroutine parse_int(text, out)
        character(len=*), intent(in) :: text
        integer(i32), intent(out) :: out
        read (text, *) out
    end subroutine parse_int

    !> `true/.true.` を真として解釈し、論理値へ変換する。
  !! @param[in] text 入力引数。
  !! @param[out] out 出力引数。
    subroutine parse_logical(text, out)
        character(len=*), intent(in) :: text
        logical, intent(out) :: out
        character(len=32) :: t
        t = lower(trim(adjustl(text)))
        out = (t == 'true' .or. t == '.true.')
    end subroutine parse_logical

    !> 文字列リテラルから前後の引用符を除去して格納する。
  !! @param[in] text 入力引数。
  !! @param[out] out 出力引数。
    subroutine parse_string(text, out)
        character(len=*), intent(in) :: text
        character(len=*), intent(out) :: out
        character(len=:), allocatable :: tmp
        tmp = trim(adjustl(text))
        if (len(tmp) >= 2 .and. tmp(1:1) == '"' .and. tmp(len(tmp):len(tmp)) == '"') then
            tmp = tmp(2:len(tmp) - 1)
        end if
        out = trim(tmp)
    end subroutine parse_string

    !> `[x, y, z]` 形式の3要素実数ベクトル文字列を配列へ変換する。
  !! @param[in] text 入力引数。
  !! @param[out] out 出力引数。
    subroutine parse_real3(text, out)
        character(len=*), intent(in) :: text
        real(dp), intent(out) :: out(3)
        character(len=256) :: t
        t = trim(adjustl(text))
        if (t(1:1) == '[') t = t(2:)
        if (t(len_trim(t):len_trim(t)) == ']') t = t(:len_trim(t) - 1)
        read (t, *) out(1), out(2), out(3)
    end subroutine parse_real3

    subroutine parse_boundary_mode(text, out)
        character(len=*), intent(in) :: text
        integer(i32), intent(out) :: out
        character(len=64) :: mode

        call parse_string(text, mode)
        select case (trim(lower(mode)))
        case ('open', 'outflow', 'escape')
            out = bc_open
        case ('reflect', 'reflection')
            out = bc_reflect
        case ('periodic')
            out = bc_periodic
        case default
            error stop 'Unknown boundary condition mode in [sim].'
        end select
    end subroutine parse_boundary_mode

    !> 行内コメント(`#`以降)を除去して、値パース対象の文字列を返す。
  !! @param[in] line 入力引数。
  !! @return out 関数の戻り値。
    pure function strip_comment(line) result(out)
        character(len=*), intent(in) :: line
        character(len=len(line)) :: out
        integer :: p
        p = index(line, '#')
        if (p > 0) then
            out = trim(line(:p - 1))
        else
            out = trim(line)
        end if
    end function strip_comment

    !> ASCII英字を小文字化し、大文字小文字非依存のキー比較を可能にする。
  !! @param[in] s 入力引数。
  !! @return o 関数の戻り値。
    pure function lower(s) result(o)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: o
        integer :: i, c
        o = s
        do i = 1, len(s)
            c = iachar(s(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) o(i:i) = achar(c + 32)
        end do
    end function lower

    !> 文字列 `s` が指定サフィックスで終わるかを判定する。
  !! @param[in] s 入力引数。
  !! @param[in] suffix 入力引数。
  !! @return ends_with 関数の戻り値。
    pure logical function ends_with(s, suffix)
        character(len=*), intent(in) :: s, suffix
        integer :: ls, lf
        ls = len_trim(s)
        lf = len_trim(suffix)
        if (lf > ls) then
            ends_with = .false.
        else
            ends_with = (s(ls - lf + 1:ls) == suffix(1:lf))
        end if
    end function ends_with

end module bem_app_config
