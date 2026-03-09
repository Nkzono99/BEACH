!> `bem_field_solver` の初期化・設定補助手続きを実装する submodule。
submodule (bem_field_solver) bem_field_solver_config
  implicit none
contains

  !> 設定値から direct/treecode の実行モードを確定し、必要なら木構造を初期化する。
  module procedure init_field_solver
    character(len=16) :: requested_mode

    self%softening = sim%softening
    self%theta = sim%tree_theta
    self%leaf_max = sim%tree_leaf_max
    self%min_nelem = sim%tree_min_nelem
    self%nelem = mesh%nelem

    requested_mode = lower_ascii(trim(sim%field_solver))
    select case (trim(requested_mode))
    case ('direct')
      self%mode = 'direct'
    case ('treecode')
      self%mode = 'treecode'
    case ('auto')
      if (mesh%nelem >= self%min_nelem) then
        self%mode = 'treecode'
        call estimate_auto_tree_params(mesh%nelem, self%theta, self%leaf_max)
      else
        self%mode = 'direct'
      end if
    case default
      error stop 'Unknown sim.field_solver in field solver init.'
    end select

    if (trim(self%mode) == 'treecode' .and. mesh%nelem > 0_i32) then
      call build_tree_topology(self, mesh)
      call refresh_field_solver(self, mesh)
    else
      call reset_tree_storage(self)
      self%tree_ready = .false.
    end if
  end procedure init_field_solver

  !> 要素数レンジに応じた treecode の `theta` と `leaf_max` の推奨値を返す。
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
