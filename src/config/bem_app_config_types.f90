!> アプリ設定の型定義と、設定由来の粒子数計算をまとめるモジュール。
module bem_app_config_types
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config
  implicit none

  integer, parameter :: max_templates = 8
  integer, parameter :: max_particle_species = 8

  !> 1粒子種の注入設定を表す。
  type :: particle_species_spec
    logical :: enabled = .false.
    integer(i32) :: npcls_per_step = 0_i32
    logical :: has_npcls_per_step = .false.
    character(len=16) :: source_mode = 'volume_seed'
    real(dp) :: number_density_cm3 = 0.0d0
    real(dp) :: number_density_m3 = 0.0d0
    logical :: has_number_density_cm3 = .false.
    logical :: has_number_density_m3 = .false.
    real(dp) :: q_particle = -1.602176634d-19
    real(dp) :: m_particle = 9.10938356d-31
    real(dp) :: w_particle = 1.0d0
    logical :: has_w_particle = .false.
    integer(i32) :: target_macro_particles_per_batch = 0_i32
    logical :: has_target_macro_particles_per_batch = .false.
    real(dp) :: pos_low(3) = [-0.4d0, -0.4d0, 0.2d0]
    real(dp) :: pos_high(3) = [0.4d0, 0.4d0, 0.5d0]
    real(dp) :: drift_velocity(3) = [0.0d0, 0.0d0, -8.0d5]
    real(dp) :: temperature_k = 2.0d4
    real(dp) :: temperature_ev = -1.0d0
    logical :: has_temperature_k = .false.
    logical :: has_temperature_ev = .false.
    real(dp) :: emit_current_density_a_m2 = 0.0d0
    integer(i32) :: rays_per_batch = 0_i32
    real(dp) :: normal_drift_speed = 0.0d0
    real(dp) :: ray_direction(3) = [0.0d0, 0.0d0, 0.0d0]
    logical :: has_ray_direction = .false.
    character(len=16) :: inject_face = ''
  end type particle_species_spec

  !> 1つのテンプレート形状の有効化フラグと幾何パラメータを保持する。
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

  !> 実行条件・メッシュ入力・粒子注入・出力設定を一元管理する。
  type :: app_config
    character(len=16) :: mesh_mode = 'auto'
    character(len=256) :: obj_path = 'examples/simple_plate.obj'
    integer(i32) :: n_templates = 1
    type(template_spec) :: templates(max_templates)

    integer(i32) :: n_particles = 0_i32
    integer(i32) :: n_particle_species = 0_i32
    type(particle_species_spec) :: particle_species(max_particle_species)

    logical :: write_output = .true.
    character(len=256) :: output_dir = 'outputs/latest'
    integer(i32) :: history_stride = 1
    logical :: resume_output = .false.

    type(sim_config) :: sim
  end type app_config

contains

  !> 有効な粒子種の `npcls_per_step` を合計し、1バッチあたりの粒子数を返す。
  !! 1つ以上の粒子種が有効で、かつ合計が正でない場合は停止する。
  !! @param[in] cfg 粒子種設定を含むアプリ設定。
  !! @return batch_n 1バッチあたりの総粒子数。
  integer(i32) function particles_per_batch_from_config(cfg) result(batch_n)
    type(app_config), intent(in) :: cfg
    integer(i32) :: s
    logical :: has_dynamic_source

    if (cfg%n_particle_species <= 0) then
      error stop 'At least one [[particles.species]] entry is required.'
    end if

    batch_n = 0_i32
    has_dynamic_source = .false.
    do s = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      select case (trim(cfg%particle_species(s)%source_mode))
      case ('volume_seed')
        if (cfg%particle_species(s)%npcls_per_step < 0_i32) then
          error stop 'particles.species.npcls_per_step must be >= 0.'
        end if
        batch_n = batch_n + cfg%particle_species(s)%npcls_per_step
      case ('reservoir_face')
        has_dynamic_source = .true.
      case ('photo_raycast')
        has_dynamic_source = .true.
      case default
        error stop 'Unknown particles.species.source_mode.'
      end select
    end do

    if (batch_n <= 0_i32 .and. .not. has_dynamic_source) then
      error stop 'At least one enabled [[particles.species]] entry must have npcls_per_step > 0.'
    end if
  end function particles_per_batch_from_config

  !> バッチ数と1バッチ粒子数から総粒子数を返す。
  !! @param[in] cfg バッチ数と粒子種設定を含むアプリ設定。
  !! @return total_n 総粒子数。
  integer(i32) function total_particles_from_config(cfg) result(total_n)
    type(app_config), intent(in) :: cfg

    if (cfg%sim%batch_count <= 0_i32) then
      error stop 'sim.batch_count must be > 0.'
    end if
    total_n = cfg%sim%batch_count * particles_per_batch_from_config(cfg)
  end function total_particles_from_config

  !> `app_config` を既定値で初期化し、TOML 上書き前の状態を作る。
  !! @param[out] cfg 既定値で初期化したアプリ設定。
  subroutine default_app_config(cfg)
    type(app_config), intent(out) :: cfg

    cfg = app_config()
    cfg%sim%dt = 1.0d-9
    cfg%sim%rng_seed = 12345_i32
    cfg%sim%batch_count = 1_i32
    cfg%sim%max_step = 400
    cfg%sim%tol_rel = 1.0d-8
    cfg%sim%softening = 1.0d-6
    cfg%sim%b0 = [0.0d0, 0.0d0, 0.0d0]
    cfg%sim%reservoir_potential_model = 'none'
    cfg%sim%phi_infty = 0.0d0
    cfg%sim%injection_face_phi_grid_n = 3_i32
    cfg%sim%raycast_max_bounce = 16_i32
    cfg%n_particles = 0_i32

    cfg%templates(1)%enabled = .true.
    cfg%templates(1)%kind = 'plane'
    cfg%templates(1)%size_x = 1.0d0
    cfg%templates(1)%size_y = 1.0d0
    cfg%templates(1)%nx = 1
    cfg%templates(1)%ny = 1
    cfg%templates(1)%center = [0.0d0, 0.0d0, 0.0d0]
  end subroutine default_app_config

  !> `[[particles.species]]` の既定値を返す。
  !! 現行仕様では、列挙された粒子種は既定で有効とみなす。
  !! @return spec 粒子種の既定設定。
  pure function species_from_defaults() result(spec)
    type(particle_species_spec) :: spec

    spec = particle_species_spec()
    spec%enabled = .true.
    spec%npcls_per_step = 0_i32
  end function species_from_defaults

end module bem_app_config_types
