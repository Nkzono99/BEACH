!> `bem_field_solver` の初期化・設定補助手続きを実装する submodule。
submodule (bem_field_solver) bem_field_solver_config
  implicit none
contains

  !> 設定値から direct/treecode/fmm の実行モードを確定し、必要なら木構造を初期化する。
  module procedure init_field_solver
    character(len=16) :: requested_mode, field_bc_mode
    integer(i32) :: axis, n_periodic
    real(dp) :: span, min_periodic_len

    self%softening = sim%softening
    self%min_nelem = sim%tree_min_nelem
    self%nelem = mesh%nelem
    self%use_periodic2 = .false.
    self%periodic_axes = 0_i32
    self%periodic_len = 0.0d0
    self%periodic_image_layers = max(0_i32, sim%field_periodic_image_layers)
    self%periodic_far_correction = lower_ascii(trim(sim%field_periodic_far_correction))
    self%periodic_ewald_alpha = sim%field_periodic_ewald_alpha
    self%periodic_ewald_layers = max(0_i32, sim%field_periodic_ewald_layers)

    requested_mode = lower_ascii(trim(sim%field_solver))
    field_bc_mode = lower_ascii(trim(sim%field_bc_mode))
    self%field_bc_mode = field_bc_mode
    select case (trim(field_bc_mode))
    case ('free')
      self%periodic_far_correction = 'none'
      continue
    case ('periodic2')
      if (trim(requested_mode) /= 'fmm') then
        error stop 'sim.field_bc_mode must be "free" unless sim.field_solver="fmm".'
      end if
      n_periodic = 0_i32
      do axis = 1_i32, 3_i32
        if ((sim%bc_low(axis) == bc_periodic) .neqv. (sim%bc_high(axis) == bc_periodic)) then
          error stop 'periodic2 requires bc_low(axis)=bc_high(axis)=periodic for periodic axes.'
        end if
        if (sim%bc_low(axis) == bc_periodic) then
          n_periodic = n_periodic + 1_i32
          if (n_periodic <= 2_i32) self%periodic_axes(n_periodic) = axis
        end if
      end do
      if (n_periodic /= 2_i32) then
        error stop 'sim.field_bc_mode="periodic2" requires exactly two periodic axes.'
      end if
      do axis = 1_i32, 2_i32
        span = sim%box_max(self%periodic_axes(axis)) - sim%box_min(self%periodic_axes(axis))
        if (span <= 0.0d0) error stop 'periodic2 requires positive box length on periodic axes.'
        self%periodic_len(axis) = span
      end do
      min_periodic_len = min(self%periodic_len(1), self%periodic_len(2))
      if (self%periodic_ewald_alpha <= 0.0d0) then
        ! 最初の省略画像殻で erfc(alpha*r) ~ 0.1 程度になるように自動設定する。
        self%periodic_ewald_alpha = 1.2d0 / (real(self%periodic_image_layers + 1_i32, dp) * min_periodic_len)
      end if
      self%use_periodic2 = .true.
    case default
      error stop 'Unknown sim.field_bc_mode in field solver init.'
    end select

    select case (trim(requested_mode))
    case ('direct')
      self%mode = 'direct'
    case ('treecode')
      self%mode = 'treecode'
    case ('fmm')
      self%mode = 'fmm'
    case ('auto')
      if (mesh%nelem >= self%min_nelem) then
        self%mode = 'treecode'
      else
        self%mode = 'direct'
      end if
    case default
      error stop 'Unknown sim.field_solver in field solver init.'
    end select

    if (trim(self%mode) == 'treecode' .or. trim(self%mode) == 'fmm') then
      call estimate_auto_tree_params(mesh%nelem, self%theta, self%leaf_max)
      if (sim%has_tree_theta) self%theta = sim%tree_theta
      if (sim%has_tree_leaf_max) self%leaf_max = sim%tree_leaf_max
    else
      self%theta = sim%tree_theta
      self%leaf_max = sim%tree_leaf_max
    end if

    if ((trim(self%mode) == 'treecode' .or. trim(self%mode) == 'fmm') .and. mesh%nelem > 0_i32) then
      call build_tree_topology(self, mesh)
      call refresh_field_solver(self, mesh)
    else
      call reset_tree_storage(self)
      self%tree_ready = .false.
    end if
  end procedure init_field_solver

  !> 要素数レンジに応じた treecode/FMM の `theta` と `leaf_max` 推奨値を返す。
  module procedure estimate_auto_tree_params
    if (nelem < 1500_i32) then
      theta = 0.40d0
      leaf_max = 12_i32
    else if (nelem < 10000_i32) then
      theta = 0.50d0
      leaf_max = 16_i32
    else if (nelem < 50000_i32) then
      theta = 0.58d0
      leaf_max = 20_i32
    else
      theta = 0.65d0
      leaf_max = 24_i32
    end if
  end procedure estimate_auto_tree_params

  !> ASCII 英字を小文字化し、設定キー比較に使う正規化文字列を返す。
  module procedure lower_ascii
    integer :: i, code

    out = s
    do i = 1, len(s)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        out(i:i) = achar(code + 32)
      end if
    end do
  end procedure lower_ascii

end submodule bem_field_solver_config
