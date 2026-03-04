!> `app_config` からメッシュ・粒子群を構築する実行時変換モジュール。
module bem_app_config_runtime
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, particles_soa
  use bem_templates, only: make_plane, make_box, make_cylinder, make_sphere
  use bem_mesh, only: init_mesh
  use bem_importers, only: load_obj_mesh
  use bem_injection, only: seed_rng, sample_uniform_positions, sample_shifted_maxwell_velocities
  use bem_particles, only: init_particles
  use bem_app_config_types, only: &
    app_config, particle_species_spec, template_spec, max_templates, particles_per_batch_from_config, &
    total_particles_from_config
  use bem_app_config_parser, only: lower
  implicit none

contains

  !> `mesh_mode` と OBJ ファイル有無に応じてメッシュ生成方法を選ぶ。
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

  !> 設定全体ぶんの粒子群を生成し、SoA へ詰める。
  !! 粒子種ごとに乱数サンプルした後、種ごとに rank を揃えて interleave する。
  subroutine init_particles_from_config(cfg, pcls)
    type(app_config), intent(in) :: cfg
    type(particles_soa), intent(out) :: pcls

    integer(i32) :: s, total_n, max_n, out_idx, rank_idx
    integer(i32), allocatable :: counts(:)
    real(dp), allocatable :: x_species(:, :, :), v_species(:, :, :)
    real(dp), allocatable :: q(:), m(:), w(:), x(:, :), v(:, :)

    if (cfg%n_particle_species <= 0) error stop 'At least one [[particles.species]] entry is required.'

    call seed_rng([cfg%sim%rng_seed])

    allocate (counts(cfg%n_particle_species))
    counts = 0_i32
    do s = 1, cfg%n_particle_species
      if (cfg%particle_species(s)%enabled) then
        if (cfg%particle_species(s)%npcls_per_step < 0_i32) then
          error stop 'particles.species.npcls_per_step must be >= 0.'
        end if
        counts(s) = cfg%sim%batch_count * cfg%particle_species(s)%npcls_per_step
      end if
    end do

    total_n = total_particles_from_config(cfg)
    max_n = max(1_i32, maxval(counts))

    allocate (x_species(3, max_n, cfg%n_particle_species))
    allocate (v_species(3, max_n, cfg%n_particle_species))
    x_species = 0.0d0
    v_species = 0.0d0

    do s = 1, cfg%n_particle_species
      if (counts(s) <= 0_i32) cycle
      call sample_species_state(cfg%particle_species(s), counts(s), x_species(:, 1:counts(s), s), v_species(:, 1:counts(s), s))
    end do

    allocate (x(3, total_n), v(3, total_n), q(total_n), m(total_n), w(total_n))
    out_idx = 0_i32
    do rank_idx = 1, maxval(counts)
      do s = 1, cfg%n_particle_species
        if (rank_idx > counts(s)) cycle
        out_idx = out_idx + 1_i32
        x(:, out_idx) = x_species(:, rank_idx, s)
        v(:, out_idx) = v_species(:, rank_idx, s)
        q(out_idx) = cfg%particle_species(s)%q_particle
        m(out_idx) = cfg%particle_species(s)%m_particle
        w(out_idx) = cfg%particle_species(s)%w_particle
      end do
    end do

    call init_particles(pcls, x, v, q, m, w)
  end subroutine init_particles_from_config

  !> バッチ生成前に乱数シードだけを初期化する。
  subroutine seed_particles_from_config(cfg)
    type(app_config), intent(in) :: cfg

    call seed_rng([cfg%sim%rng_seed])
  end subroutine seed_particles_from_config

  !> 指定バッチ番号に対応する粒子バッチを生成する。
  subroutine init_particle_batch_from_config(cfg, batch_idx, pcls)
    type(app_config), intent(in) :: cfg
    integer(i32), intent(in) :: batch_idx
    type(particles_soa), intent(out) :: pcls

    integer(i32) :: s, i, batch_n, max_rank, out_idx
    integer(i32), allocatable :: counts(:), species_cursor(:), species_id(:)
    real(dp), allocatable :: x_species(:, :, :), v_species(:, :, :), x(:, :), v(:, :), q(:), m(:), w(:)

    if (cfg%sim%batch_count <= 0_i32) error stop 'sim.batch_count must be > 0.'
    if (batch_idx < 1_i32 .or. batch_idx > cfg%sim%batch_count) then
      error stop 'Requested batch index is out of range.'
    end if

    batch_n = particles_per_batch_from_config(cfg)

    allocate (counts(cfg%n_particle_species))
    counts = 0_i32
    do s = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      counts(s) = cfg%particle_species(s)%npcls_per_step
    end do

    allocate (species_id(batch_n))
    max_rank = max(1_i32, maxval(counts))
    out_idx = 0_i32
    do i = 1, max_rank
      do s = 1, cfg%n_particle_species
        if (i > counts(s)) cycle
        out_idx = out_idx + 1_i32
        species_id(out_idx) = s
      end do
    end do

    allocate (x_species(3, max_rank, cfg%n_particle_species))
    allocate (v_species(3, max_rank, cfg%n_particle_species))
    x_species = 0.0d0
    v_species = 0.0d0
    do s = 1, cfg%n_particle_species
      if (counts(s) <= 0_i32) cycle
      call sample_species_state(cfg%particle_species(s), counts(s), x_species(:, 1:counts(s), s), v_species(:, 1:counts(s), s))
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

  !> 1粒子種ぶんの位置・速度サンプルをまとめて生成する。
  subroutine sample_species_state(spec, n, x, v)
    type(particle_species_spec), intent(in) :: spec
    integer(i32), intent(in) :: n
    real(dp), intent(out) :: x(:, :)
    real(dp), intent(out) :: v(:, :)

    if (n <= 0_i32) return
    call sample_uniform_positions(spec%pos_low, spec%pos_high, x)
    call sample_shifted_maxwell_velocities(spec%drift_velocity, spec%m_particle, v, temperature_k=spec%temperature_k)
  end subroutine sample_species_state

  !> 有効なテンプレートを連結し、1つのメッシュへまとめる。
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

  !> テンプレート種別に応じて形状生成ルーチンへディスパッチする。
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
        mesh, radius=spec%radius, height=spec%height, n_theta=spec%n_theta, n_z=spec%n_z, &
        cap=spec%cap, center=spec%center &
      )
    case ('sphere')
      call make_sphere(mesh, radius=spec%radius, n_lon=spec%n_lon, n_lat=spec%n_lat, center=spec%center)
    case default
      error stop 'Unknown template kind in config.'
    end select
  end subroutine build_one_template

  !> 既存三角形配列へ追加分を連結し、再確保後の配列へ差し替える。
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

end module bem_app_config_runtime
