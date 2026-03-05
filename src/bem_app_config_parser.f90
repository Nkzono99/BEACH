!> TOML風設定ファイルを `app_config` へ読み込む軽量パーサ。
module bem_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_types, only: bc_open, bc_reflect, bc_periodic
  use bem_app_config_types, only: &
    app_config, particle_species_spec, template_spec, max_templates, max_particle_species, species_from_defaults
  implicit none

contains

  !> `.toml` 拡張子の設定ファイルを読み込み、既存値へ上書き適用する。
  !! @param[in] path 読み込む設定ファイルパス（`.toml` 必須）。
  !! @param[inout] cfg 読み込み結果で上書きするアプリ設定。
  subroutine load_app_config(path, cfg)
    character(len=*), intent(in) :: path
    type(app_config), intent(inout) :: cfg

    if (.not. ends_with(lower(trim(path)), '.toml')) then
      error stop 'Only TOML config is supported. Please pass a .toml file.'
    end if
    call load_toml_config(path, cfg)
  end subroutine load_app_config

  !> 最小限の TOML セクションと `key = value` を解釈して設定へ反映する。
  !! 現在は `sim` / `mesh` / `output` / `[[mesh.templates]]` / `[[particles.species]]` を扱う。
  !! @param[in] path 読み込むTOMLファイルパス。
  !! @param[inout] cfg 読み込み結果で更新するアプリ設定。
  subroutine load_toml_config(path, cfg)
    character(len=*), intent(in) :: path
    type(app_config), intent(inout) :: cfg
    integer :: u, ios, i, t_idx, s_idx
    integer(i32) :: per_batch_particles
    logical :: has_reservoir_species
    character(len=512) :: raw, line, section

    t_idx = 0
    s_idx = 0
    section = ''

    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Could not open TOML file.'

    do
      read (u, '(A)', iostat=ios) raw
      if (ios /= 0) exit
      line = strip_comment(trim(raw))
      if (len_trim(line) == 0) cycle

      if (line(1:1) == '[') then
        ! 配列テーブルは専用カウンタを進める。
        if (trim(line) == '[[mesh.templates]]') then
          t_idx = t_idx + 1
          if (t_idx > max_templates) error stop 'Too many mesh.templates entries.'
          cfg%templates(t_idx)%enabled = .true.
          section = 'mesh.template'
        else if (trim(line) == '[[particles.species]]') then
          s_idx = s_idx + 1
          if (s_idx > max_particle_species) error stop 'Too many particles.species entries.'
          cfg%particle_species(s_idx) = species_from_defaults()
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
        call apply_particles_kv(line)
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
    if (cfg%sim%batch_count <= 0_i32) error stop 'sim.batch_count must be > 0.'
    if (s_idx <= 0) error stop 'At least one [[particles.species]] entry is required.'

    cfg%n_particle_species = s_idx
    per_batch_particles = 0_i32
    has_reservoir_species = .false.
    do i = 1, s_idx
      if (.not. cfg%particle_species(i)%enabled) cycle

      cfg%particle_species(i)%source_mode = lower(trim(cfg%particle_species(i)%source_mode))
      select case (trim(cfg%particle_species(i)%source_mode))
      case ('volume_seed')
        if (cfg%particle_species(i)%npcls_per_step < 0_i32) then
          error stop 'particles.species.npcls_per_step must be >= 0.'
        end if
        per_batch_particles = per_batch_particles + cfg%particle_species(i)%npcls_per_step
      case ('reservoir_face')
        has_reservoir_species = .true.
        call validate_reservoir_species(cfg, i)
      case default
        error stop 'Unknown particles.species.source_mode.'
      end select
    end do

    if (per_batch_particles <= 0_i32 .and. .not. has_reservoir_species) then
      error stop 'At least one enabled [[particles.species]] entry must have npcls_per_step > 0.'
    end if
    cfg%n_particles = cfg%sim%batch_count * per_batch_particles
  end subroutine load_toml_config

  !> `[sim]` セクションのキーを `sim_config` へ適用する。
  !! @param[inout] cfg 更新対象のアプリ設定。
  !! @param[in] line `key = value` 形式の設定行。
  subroutine apply_sim_kv(cfg, line)
    type(app_config), intent(inout) :: cfg
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v

    call split_key_value(line, k, v)
    select case (trim(k))
    case ('dt')
      call parse_real(v, cfg%sim%dt)
    case ('rng_seed')
      call parse_int(v, cfg%sim%rng_seed)
    case ('batch_count')
      call parse_int(v, cfg%sim%batch_count)
    case ('batch_duration')
      call parse_real(v, cfg%sim%batch_duration)
    case ('npcls_per_step')
      error stop 'sim.npcls_per_step was removed. Use sim.batch_count instead.'
    case ('max_step')
      call parse_int(v, cfg%sim%max_step)
    case ('tol_rel')
      call parse_real(v, cfg%sim%tol_rel)
    case ('q_floor')
      call parse_real(v, cfg%sim%q_floor)
    case ('softening')
      call parse_real(v, cfg%sim%softening)
    case ('use_hybrid')
      call parse_logical(v, cfg%sim%use_hybrid)
    case ('r_switch_factor')
      call parse_real(v, cfg%sim%r_switch_factor)
    case ('n_sub')
      call parse_int(v, cfg%sim%n_sub)
    case ('softening_factor')
      call parse_real(v, cfg%sim%softening_factor)
    case ('b0')
      call parse_real3(v, cfg%sim%b0)
    case ('use_box')
      call parse_logical(v, cfg%sim%use_box)
    case ('box_min')
      call parse_real3(v, cfg%sim%box_min)
    case ('box_max')
      call parse_real3(v, cfg%sim%box_max)
    case ('bc_x_low')
      call parse_boundary_mode(v, cfg%sim%bc_low(1))
    case ('bc_x_high')
      call parse_boundary_mode(v, cfg%sim%bc_high(1))
    case ('bc_y_low')
      call parse_boundary_mode(v, cfg%sim%bc_low(2))
    case ('bc_y_high')
      call parse_boundary_mode(v, cfg%sim%bc_high(2))
    case ('bc_z_low')
      call parse_boundary_mode(v, cfg%sim%bc_low(3))
    case ('bc_z_high')
      call parse_boundary_mode(v, cfg%sim%bc_high(3))
    end select
  end subroutine apply_sim_kv

  !> 廃止済み `[particles]` セクションを検出し、新仕様への移行を促す。
  !! @param[in] line 廃止セクション内の `key = value` 行。
  subroutine apply_particles_kv(line)
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v

    call split_key_value(line, k, v)
    select case (trim(k))
    case ('rng_seed')
      error stop 'particles.rng_seed was moved to sim.rng_seed.'
    case ('n_particles')
      error stop 'particles.n_particles was removed. Define [[particles.species]] entries with npcls_per_step instead.'
    case ('q_particle', 'm_particle', 'w_particle', 'pos_low', 'pos_high', 'drift_velocity', 'temperature_k')
      error stop '[particles] defaults were removed. Define particle properties inside each [[particles.species]] entry.'
    case default
      error stop '[particles] defaults were removed. Use sim.rng_seed and [[particles.species]] instead.'
    end select
  end subroutine apply_particles_kv

  !> `[[particles.species]]` のキーを粒子種設定へ適用する。
  !! @param[inout] spec 更新対象の粒子種設定。
  !! @param[in] line `key = value` 形式の設定行。
  subroutine apply_particles_species_kv(spec, line)
    type(particle_species_spec), intent(inout) :: spec
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v

    call split_key_value(line, k, v)
    select case (trim(k))
    case ('enabled')
      call parse_logical(v, spec%enabled)
    case ('npcls_per_step')
      call parse_int(v, spec%npcls_per_step)
      spec%has_npcls_per_step = .true.
    case ('source_mode')
      call parse_string(v, spec%source_mode)
      spec%source_mode = lower(trim(spec%source_mode))
    case ('number_density_cm3')
      call parse_real(v, spec%number_density_cm3)
      spec%has_number_density_cm3 = .true.
    case ('number_density_m3')
      call parse_real(v, spec%number_density_m3)
      spec%has_number_density_m3 = .true.
    case ('n_particles')
      error stop 'particles.species.n_particles was removed. Use particles.species.npcls_per_step instead.'
    case ('q_particle')
      call parse_real(v, spec%q_particle)
    case ('m_particle')
      call parse_real(v, spec%m_particle)
    case ('w_particle')
      call parse_real(v, spec%w_particle)
    case ('pos_low')
      call parse_real3(v, spec%pos_low)
    case ('pos_high')
      call parse_real3(v, spec%pos_high)
    case ('drift_velocity')
      call parse_real3(v, spec%drift_velocity)
    case ('temperature_k')
      call parse_real(v, spec%temperature_k)
      spec%has_temperature_k = .true.
    case ('temperature_ev')
      call parse_real(v, spec%temperature_ev)
      spec%has_temperature_ev = .true.
    case ('inject_face')
      call parse_string(v, spec%inject_face)
      spec%inject_face = lower(trim(spec%inject_face))
    end select
  end subroutine apply_particles_species_kv

  !> `reservoir_face` 用の必須項目と整合性を検証する。
  !! @param[in] cfg 検証対象のアプリ設定。
  !! @param[in] species_idx 検証する粒子種のインデックス（1始まり）。
  subroutine validate_reservoir_species(cfg, species_idx)
    type(app_config), intent(in) :: cfg
    integer, intent(in) :: species_idx

    integer :: axis, axis_t1, axis_t2
    real(dp) :: boundary_value
    type(particle_species_spec) :: spec

    spec = cfg%particle_species(species_idx)

    if (spec%has_npcls_per_step) then
      error stop 'particles.species.npcls_per_step is auto-computed for reservoir_face.'
    end if
    if (spec%w_particle <= 0.0d0) then
      error stop 'particles.species.w_particle must be > 0 for reservoir_face.'
    end if
    if (.not. cfg%sim%use_box) then
      error stop 'particles.species.source_mode="reservoir_face" requires sim.use_box = true.'
    end if
    if (cfg%sim%batch_duration <= 0.0d0) then
      error stop 'sim.batch_duration must be > 0 for reservoir_face.'
    end if
    if (spec%has_number_density_cm3 .and. spec%has_number_density_m3) then
      error stop 'Specify either number_density_cm3 or number_density_m3, not both.'
    end if
    if (.not. spec%has_number_density_cm3 .and. .not. spec%has_number_density_m3) then
      error stop 'reservoir_face requires number_density_cm3 or number_density_m3.'
    end if
    if (spec%has_number_density_cm3) then
      if (spec%number_density_cm3 <= 0.0d0) error stop 'number_density_cm3 must be > 0.'
    else
      if (spec%number_density_m3 <= 0.0d0) error stop 'number_density_m3 must be > 0.'
    end if
    if (spec%has_temperature_ev .and. spec%has_temperature_k) then
      error stop 'Specify either temperature_ev or temperature_k, not both.'
    end if
    if (spec%has_temperature_ev) then
      if (spec%temperature_ev < 0.0d0) error stop 'temperature_ev must be >= 0.'
    else if (spec%has_temperature_k) then
      if (spec%temperature_k < 0.0d0) error stop 'temperature_k must be >= 0.'
    else
      if (spec%temperature_k < 0.0d0) error stop 'temperature_k must be >= 0.'
    end if

    call resolve_inject_face(cfg%sim%box_min, cfg%sim%box_max, spec%inject_face, axis, boundary_value)
    axis_t1 = modulo(axis, 3) + 1
    axis_t2 = modulo(axis + 1, 3) + 1
    if (abs(spec%pos_low(axis) - boundary_value) > 1.0d-12 .or. abs(spec%pos_high(axis) - boundary_value) > 1.0d-12) then
      error stop 'reservoir_face pos_low/pos_high must lie on the selected box face.'
    end if
    if (spec%pos_high(axis_t1) < spec%pos_low(axis_t1) .or. spec%pos_high(axis_t2) < spec%pos_low(axis_t2)) then
      error stop 'reservoir_face tangential bounds must satisfy pos_high >= pos_low.'
    end if
    if ((spec%pos_high(axis_t1) - spec%pos_low(axis_t1)) <= 0.0d0 .or. &
        (spec%pos_high(axis_t2) - spec%pos_low(axis_t2)) <= 0.0d0) then
      error stop 'reservoir_face opening area must be positive.'
    end if
  end subroutine validate_reservoir_species

  !> 注入面名から法線軸と対応する境界座標を返す。
  !! @param[in] box_min ボックス下限座標 `(x,y,z)` [m]。
  !! @param[in] box_max ボックス上限座標 `(x,y,z)` [m]。
  !! @param[in] inject_face 注入面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[out] axis 注入面の法線軸インデックス（1:x, 2:y, 3:z）。
  !! @param[out] boundary_value 指定面の境界座標値 [m]。
  subroutine resolve_inject_face(box_min, box_max, inject_face, axis, boundary_value)
    real(dp), intent(in) :: box_min(3), box_max(3)
    character(len=*), intent(in) :: inject_face
    integer, intent(out) :: axis
    real(dp), intent(out) :: boundary_value

    select case (trim(lower(inject_face)))
    case ('x_low')
      axis = 1
      boundary_value = box_min(1)
    case ('x_high')
      axis = 1
      boundary_value = box_max(1)
    case ('y_low')
      axis = 2
      boundary_value = box_min(2)
    case ('y_high')
      axis = 2
      boundary_value = box_max(2)
    case ('z_low')
      axis = 3
      boundary_value = box_min(3)
    case ('z_high')
      axis = 3
      boundary_value = box_max(3)
    case default
      error stop 'Unknown particles.species.inject_face.'
    end select
  end subroutine resolve_inject_face

  !> `[mesh]` セクションのキーをメッシュ入力設定へ適用する。
  !! @param[inout] cfg 更新対象のアプリ設定。
  !! @param[in] line `key = value` 形式の設定行。
  subroutine apply_mesh_kv(cfg, line)
    type(app_config), intent(inout) :: cfg
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v

    call split_key_value(line, k, v)
    select case (trim(k))
    case ('mode')
      call parse_string(v, cfg%mesh_mode)
    case ('obj_path')
      call parse_string(v, cfg%obj_path)
    case ('n_templates')
      call parse_int(v, cfg%n_templates)
    end select
  end subroutine apply_mesh_kv

  !> `[[mesh.templates]]` のキーをテンプレート設定へ適用する。
  !! @param[inout] spec 更新対象のテンプレート設定。
  !! @param[in] line `key = value` 形式の設定行。
  subroutine apply_template_kv(spec, line)
    type(template_spec), intent(inout) :: spec
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v

    call split_key_value(line, k, v)
    select case (trim(k))
    case ('enabled')
      call parse_logical(v, spec%enabled)
    case ('kind')
      call parse_string(v, spec%kind)
    case ('center')
      call parse_real3(v, spec%center)
    case ('size_x')
      call parse_real(v, spec%size_x)
    case ('size_y')
      call parse_real(v, spec%size_y)
    case ('size')
      call parse_real3(v, spec%size)
    case ('nx')
      call parse_int(v, spec%nx)
    case ('ny')
      call parse_int(v, spec%ny)
    case ('nz')
      call parse_int(v, spec%nz)
    case ('radius')
      call parse_real(v, spec%radius)
    case ('height')
      call parse_real(v, spec%height)
    case ('n_theta')
      call parse_int(v, spec%n_theta)
    case ('n_z')
      call parse_int(v, spec%n_z)
    case ('cap')
      call parse_logical(v, spec%cap)
    case ('n_lon')
      call parse_int(v, spec%n_lon)
    case ('n_lat')
      call parse_int(v, spec%n_lat)
    end select
  end subroutine apply_template_kv

  !> `[output]` セクションのキーを出力制御設定へ適用する。
  !! @param[inout] cfg 更新対象のアプリ設定。
  !! @param[in] line `key = value` 形式の設定行。
  subroutine apply_output_kv(cfg, line)
    type(app_config), intent(inout) :: cfg
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v

    call split_key_value(line, k, v)
    select case (trim(k))
    case ('write_files')
      call parse_logical(v, cfg%write_output)
    case ('dir')
      call parse_string(v, cfg%output_dir)
    case ('history_stride')
      call parse_int(v, cfg%history_stride)
    case ('resume')
      call parse_logical(v, cfg%resume_output)
    end select
  end subroutine apply_output_kv

  !> `key = value` 形式の1行を分割し、前後空白を取り除いて返す。
  !! @param[in] line 分割対象の入力行。
  !! @param[out] key 正規化したキー文字列（小文字化済み）。
  !! @param[out] value 前後空白を除去した値文字列。
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

  !> 文字列表現を倍精度実数へ変換する。
  !! @param[in] text 実数を表す文字列。
  !! @param[out] out 変換後の倍精度実数。
  subroutine parse_real(text, out)
    character(len=*), intent(in) :: text
    real(dp), intent(out) :: out

    read (text, *) out
  end subroutine parse_real

  !> 文字列表現を 32bit 整数へ変換する。
  !! @param[in] text 整数を表す文字列。
  !! @param[out] out 変換後の32bit整数。
  subroutine parse_int(text, out)
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: out

    read (text, *) out
  end subroutine parse_int

  !> `true/.true.` を真とみなして論理値へ変換する。
  !! @param[in] text 論理値文字列。
  !! @param[out] out 変換後の論理値。
  subroutine parse_logical(text, out)
    character(len=*), intent(in) :: text
    logical, intent(out) :: out
    character(len=32) :: t

    t = lower(trim(adjustl(text)))
    out = (t == 'true' .or. t == '.true.')
  end subroutine parse_logical

  !> 前後の二重引用符を外して文字列値へ変換する。
  !! @param[in] text 文字列リテラルまたは裸文字列。
  !! @param[out] out 引用符を除去した文字列値。
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

  !> `[x, y, z]` 形式の3成分ベクトルを配列へ変換する。
  !! @param[in] text 3成分ベクトル文字列。
  !! @param[out] out 変換後の3成分実数ベクトル。
  subroutine parse_real3(text, out)
    character(len=*), intent(in) :: text
    real(dp), intent(out) :: out(3)
    character(len=256) :: t

    t = trim(adjustl(text))
    if (t(1:1) == '[') t = t(2:)
    if (t(len_trim(t):len_trim(t)) == ']') t = t(:len_trim(t) - 1)
    read (t, *) out(1), out(2), out(3)
  end subroutine parse_real3

  !> 境界条件の文字列を内部定数へ変換する。
  !! 未知のモードは実行継続せず停止する。
  !! @param[in] text 境界条件モード文字列。
  !! @param[out] out 変換後の境界条件定数。
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

  !> `#` 以降の行内コメントを除去する。
  !! @param[in] line コメント除去対象の入力行。
  !! @return out `#` 以降を取り除いた行文字列。
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

  !> ASCII 英字だけを小文字化し、設定キー比較を簡単にする。
  !! @param[in] s 変換対象文字列。
  !! @return o 小文字化した文字列。
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

  !> 文字列 `s` が `suffix` で終わるかを返す。
  !! @param[in] s 判定対象文字列。
  !! @param[in] suffix 末尾一致を判定する接尾辞。
  !! @return ends_with `s` が `suffix` で終われば `.true.`。
  pure logical function ends_with(s, suffix)
    character(len=*), intent(in) :: s
    character(len=*), intent(in) :: suffix
    integer :: ls, lf

    ls = len_trim(s)
    lf = len_trim(suffix)
    if (lf > ls) then
      ends_with = .false.
    else
      ends_with = (s(ls - lf + 1:ls) == suffix(1:lf))
    end if
  end function ends_with

end module bem_app_config_parser
