!> アプリ設定の型定義と、設定由来の粒子数計算をまとめるモジュール。
module bem_app_config_types
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config
  use bem_filesystem, only: default_periodic_cache_directory
  use bem_physics_config_types, only: &
    field_physics_config, periodic2_physics_config, panel_kernel_config, &
    normalize_legacy_physics_config
  implicit none

  !> 設定配列の初期確保サイズ（上限ではない）。
  integer, parameter :: max_templates = 8
  integer, parameter :: max_particle_species = 8
  integer(i32), parameter :: particle_bc_inherit = -1_i32
  integer(i32), parameter :: particle_inflow_none = 0_i32
  integer(i32), parameter :: particle_inflow_reservoir = 1_i32

  !> 外部シースの零電流根から species 別の表面電流と境界写像を与えるモデルの設定。
  type :: surface_current_model_config
    character(len=32) :: model = 'none'
    character(len=16) :: zhao_branch = 'auto'
    character(len=64) :: electron_species = ''
    character(len=64) :: ion_species = ''
    character(len=64) :: photoelectron_species = ''
    character(len=64) :: secondary_electron_species = ''
    real(dp) :: solar_elevation_deg = 0.0_dp
    real(dp) :: photoelectron_ref_density_m3 = 0.0_dp
    real(dp) :: photoelectron_source_scale = 1.0_dp
    real(dp) :: reference_area_m2 = 0.0_dp
    logical :: has_reference_area_m2 = .false.
    !> zhao_stationary の外部シース源を観測PE流出で解き直す accepted batch 間隔。0 は初期根を固定する。
    integer(i32) :: outflow_refresh_batches = 0_i32
    !> 電子drift時のA/C根で、上流に接するE^2<0の幅を |phi_m|(A) / |phi_0|(C) の何割まで許すか。0 は厳密判定。
    real(dp) :: zhao_upstream_band_tolerance = 0.0_dp
    !> 上流境界を「場がこの外部レベル [V/m] まで落ちた点」に置く。電子drift時のA/C根で、帯の内側の場の不足
    !! sqrt(-min E^2) がこれ以下なら受理する。0 はこの判定を使わない。
    real(dp) :: zhao_upstream_field_floor_v_m = 0.0_dp
    !> 外部根のイオン流体の圧力係数 (Sana & Mishra 2026, Eq. 17)。1 は等温。イオン温度は ion species の温度を使う。
    real(dp) :: zhao_ion_pressure_factor = 1.0_dp
  end type surface_current_model_config

  !> 1粒子種の注入設定を表す。
  type :: particle_species_spec
    character(len=64) :: species_key = ''
    logical :: enabled = .false.
    integer(i32) :: npcls_per_step = 0_i32
    character(len=32) :: source_mode = 'volume_seed'
    logical :: secondary_emission_enabled = .false.
    character(len=16) :: secondary_emission_model = 'whipple'
    real(dp) :: secondary_emission_yield_max = 1.5_dp
    real(dp) :: secondary_emission_energy_max_ev = 500.0_dp
    real(dp) :: secondary_emission_const_yield = 0.0_dp
    character(len=64) :: secondary_emission_target_species = ''
    integer(i32) :: secondary_emission_target_idx = 0_i32
    logical :: has_secondary_forbidden_source_option = .false.
    logical :: has_surface_charge_closure = .false.
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
    logical :: has_temperature_ev = .false.
    real(dp) :: emit_current_density_a_m2 = 0.0d0
    integer(i32) :: rays_per_batch = 0_i32
    logical :: deposit_opposite_charge_on_emit = .false.
    real(dp) :: normal_drift_speed = 0.0d0
    real(dp) :: ray_direction(3) = [0.0d0, 0.0d0, 0.0d0]
    logical :: has_ray_direction = .false.
    character(len=16) :: inject_face = ''
    integer(i32) :: boundary_low(3) = [particle_bc_inherit, particle_bc_inherit, particle_bc_inherit]
    integer(i32) :: boundary_high(3) = [particle_bc_inherit, particle_bc_inherit, particle_bc_inherit]
    integer(i32) :: boundary_inflow_low(3) = particle_inflow_none
    integer(i32) :: boundary_inflow_high(3) = particle_inflow_none
    character(len=16) :: surface_charge_closure = 'explicit'
    real(dp) :: target_absorbed_current_a = 0.0_dp
    logical :: has_target_absorbed_current_a = .false.
    real(dp) :: target_emission_current_a = 0.0_dp
    logical :: has_target_emission_current_a = .false.
  end type particle_species_spec

  !> OBJ または template 1 要素の表面と幾何パラメータを保持する。
  type :: template_spec
    logical :: enabled = .false.
    character(len=16) :: kind = ''
    character(len=256) :: path = ''
    real(dp) :: scale = 1.0_dp
    real(dp) :: rotation(3) = 0.0_dp
    real(dp) :: offset(3) = 0.0_dp
    character(len=16) :: surface_model = 'insulator'
    character(len=32) :: surface_side_policy = 'not_specified'
    real(dp) :: epsilon_r = 1.0d0
    real(dp) :: center(3) = 0.0d0
    real(dp) :: size_x = 1.0d0
    real(dp) :: size_y = 1.0d0
    real(dp) :: size(3) = [1.0d0, 1.0d0, 1.0d0]
    integer(i32) :: nx = 1
    integer(i32) :: ny = 1
    integer(i32) :: nz = 1
    real(dp) :: radius = 0.5d0
    real(dp) :: inner_radius = 0.25d0
    real(dp) :: height = 1.0d0
    integer(i32) :: n_theta = 24
    integer(i32) :: n_r = 4
    integer(i32) :: n_z = 1
    logical :: cap = .true.
    logical :: cap_top = .true.
    logical :: cap_bottom = .true.
    logical :: has_cap_top = .false.
    logical :: has_cap_bottom = .false.
    integer(i32) :: n_lon = 24
    integer(i32) :: n_lat = 12
  end type template_spec

  !> 実行条件・メッシュ入力・粒子注入・出力設定を一元管理する。
  type :: app_config
    integer(i32) :: n_mesh_parts = 0_i32
    type(template_spec), allocatable :: mesh_parts(:)

    integer(i32) :: n_particle_species = 0_i32
    type(particle_species_spec), allocatable :: particle_species(:)
    integer(i32) :: secondary_emission_max_generations = 3_i32
    integer(i32) :: particle_boundary_low(3) = &
                    [particle_bc_inherit, particle_bc_inherit, particle_bc_inherit]
    integer(i32) :: particle_boundary_high(3) = &
                    [particle_bc_inherit, particle_bc_inherit, particle_bc_inherit]

    logical :: write_output = .true.
    logical :: write_mesh_potential = .false.
    logical :: write_potential_history = .false.
    character(len=256) :: output_dir = 'outputs/latest'
    character(len=256) :: output_restart_from = ''
    integer(i32) :: history_stride = 1
    integer(i32) :: checkpoint_stride = 0
    logical :: resume_output = .false.

    type(sim_config) :: sim
    type(field_physics_config) :: field
    type(periodic2_physics_config) :: periodic2
    type(panel_kernel_config) :: panel
    type(surface_current_model_config) :: surface_current
  end type app_config

contains

  !> 診断に使う公開 mesh 配列の添字付きパスを返す。
  function mesh_part_config_path(part_idx) result(path)
    integer, intent(in) :: part_idx
    character(len=:), allocatable :: path
    character(len=32) :: index_text

    write (index_text, '(i0)') part_idx
    path = 'mesh.parts['//trim(index_text)//']'
  end function mesh_part_config_path

  !> 診断に使う公開 species 配列の添字付きパスを返す。
  function species_config_path(species_idx) result(path)
    integer, intent(in) :: species_idx
    character(len=:), allocatable :: path
    character(len=32) :: index_text

    write (index_text, '(i0)') species_idx
    path = 'particles.species['//trim(index_text)//']'
  end function species_config_path

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
          error stop species_config_path(int(s))//'.sampling.volume_macro_particles_per_batch must be >= 0.'
        end if
        batch_n = batch_n + cfg%particle_species(s)%npcls_per_step
      case ('photo_raycast')
        has_dynamic_source = .true.
      case ('secondary_emission')
        cycle
      case default
        error stop 'Unknown '//species_config_path(int(s))//'.source.mode.'
      end select
      has_dynamic_source = has_dynamic_source .or. &
                           any(cfg%particle_species(s)%boundary_inflow_low == particle_inflow_reservoir) .or. &
                           any(cfg%particle_species(s)%boundary_inflow_high == particle_inflow_reservoir)
    end do

    if (batch_n <= 0_i32 .and. .not. has_dynamic_source) then
      error stop 'At least one enabled [[particles.species]] entry must have sampling.volume_macro_particles_per_batch > 0.'
    end if
  end function particles_per_batch_from_config

  !> `app_config` を既定値で初期化し、TOML 上書き前の状態を作る。
  !! @param[out] cfg 既定値で初期化したアプリ設定。
  subroutine default_app_config(cfg)
    type(app_config), intent(out) :: cfg

    cfg%n_mesh_parts = 0_i32
    cfg%n_particle_species = 0_i32
    cfg%write_output = .true.
    cfg%write_mesh_potential = .false.
    cfg%write_potential_history = .false.
    cfg%output_dir = 'outputs/latest'
    cfg%output_restart_from = ''
    cfg%history_stride = 1_i32
    cfg%checkpoint_stride = 0_i32
    cfg%resume_output = .false.
    cfg%sim = sim_config()
    cfg%sim%field_periodic_cache_dir = default_periodic_cache_directory(missing_ok=.true.)
    allocate (cfg%mesh_parts(max_templates))
    allocate (cfg%particle_species(max_particle_species))
    cfg%particle_species = particle_species_spec()
    cfg%sim%dt = 1.0d-9
    cfg%sim%rng_seed = 12345_i32
    cfg%sim%batch_count = 1_i32
    cfg%sim%max_step = 400
    cfg%sim%tol_rel = 1.0d-8
    cfg%sim%field_solver = 'auto'
    cfg%sim%field_normalization = 'si'
    cfg%sim%field_length_scale = 1.0d0
    cfg%sim%field_bc_mode = 'free'
    cfg%sim%field_periodic_image_layers = 1_i32
    cfg%sim%field_periodic_far_correction = 'none'
    cfg%sim%field_periodic_ewald_alpha = 0.0d0
    cfg%sim%field_periodic_ewald_layers = 4_i32
    cfg%sim%tree_theta = 0.5d0
    cfg%sim%tree_leaf_max = 16_i32
    cfg%sim%tree_min_nelem = 256_i32
    cfg%sim%e0 = [0.0d0, 0.0d0, 0.0d0]
    cfg%sim%has_e0_vector = .false.
    cfg%sim%e0_abs = 0.0d0
    cfg%sim%has_e0_abs = .false.
    cfg%sim%e0_phi_xy_deg = 0.0d0
    cfg%sim%has_e0_phi_xy_deg = .false.
    cfg%sim%e0_phi_z_deg = 0.0d0
    cfg%sim%has_e0_phi_z_deg = .false.
    cfg%sim%b0 = [0.0d0, 0.0d0, 0.0d0]
    cfg%sim%reservoir_potential_model = 'none'
    cfg%sim%phi_infty = 0.0d0
    cfg%sim%injection_face_phi_grid_n = 3_i32
    cfg%sim%raycast_max_bounce = 16_i32
    cfg%n_mesh_parts = 1_i32

    cfg%mesh_parts(1)%enabled = .true.
    cfg%mesh_parts(1)%kind = 'plane'
    cfg%mesh_parts(1)%size_x = 1.0d0
    cfg%mesh_parts(1)%size_y = 1.0d0
    cfg%mesh_parts(1)%nx = 1
    cfg%mesh_parts(1)%ny = 1
    cfg%mesh_parts(1)%center = [0.0d0, 0.0d0, 0.0d0]
    cfg%mesh_parts(1)%surface_side_policy = 'normal_plus'
    call normalize_legacy_physics_config(cfg%sim, cfg%field, cfg%periodic2, cfg%panel)
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
