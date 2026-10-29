!> 粒子SoAデータ構造の初期化を提供するモジュール。
module bem_particles
  use bem_kinds, only: dp, i32
  use bem_types, only: particles_soa
  implicit none
contains

  !> 粒子数に一致するSoA領域を確保する。位置・速度・電荷・質量は呼出し元で埋める。
  !! SoA の配列長は容量で、有効な粒子は常に先頭 `pcls%n` 個（`extend_particles` 後は配列長の方が長くなり得る）。
  !! 重みは1、粒子種IDは0、放出元は-1、aliveはtrueで初期化する。
  subroutine allocate_particles(pcls, n)
    type(particles_soa), intent(out) :: pcls
    integer(i32), intent(in) :: n

    if (n < 0_i32) error stop "n must be non-negative"
    pcls%n = n
    allocate (pcls%x(3, n), pcls%v(3, n), pcls%q(n), pcls%m(n), pcls%w(n), &
              pcls%species_id(n), pcls%source_element(n), pcls%alive(n))
    pcls%w = 1.0_dp
    pcls%species_id = 0_i32
    pcls%source_element = -1_i32
    pcls%alive = .true.
  end subroutine allocate_particles

  !> 位置・速度・電荷・質量(と任意重み)配列から `particles_soa` を検証付きで構築する。
  !! @param[out] pcls 検証済み配列を内部に保持した `particles_soa` 構造体。
  !! @param[in] x 粒子位置配列 `x(3,n)` [m]。
  !! @param[in] v 粒子速度配列 `v(3,n)` [m/s]。
  !! @param[in] q 粒子電荷配列 `q(n)` [C]。
  !! @param[in] m 粒子質量配列 `m(n)` [kg]。
  !! @param[in] w マクロ粒子重み配列 `w(n)`（省略時は1）。
  !! @param[in] species_id 粒子種ID `species_id(n)`（省略時は0）。
  !! @param[in] source_element 表面放出元要素 `source_element(n)`（省略時は-1）。
  subroutine init_particles(pcls, x, v, q, m, w, species_id, source_element)
    type(particles_soa), intent(out) :: pcls
    real(dp), intent(in) :: x(:, :), v(:, :), q(:), m(:)
    real(dp), intent(in), optional :: w(:)
    integer(i32), intent(in), optional :: species_id(:)
    integer(i32), intent(in), optional :: source_element(:)
    integer(i32) :: n

    n = size(q)
    if (size(x, 1) /= 3 .or. size(v, 1) /= 3) then
      error stop "particle input first dimension must be 3"
    end if
    if (size(x, 2) /= n .or. size(v, 2) /= n .or. size(m) /= n) then
      error stop "particle input size mismatch"
    end if

    call allocate_particles(pcls, n)
    pcls%x = x
    pcls%v = v
    pcls%q = q
    pcls%m = m
    if (present(w)) then
      if (size(w) /= n) error stop "w size mismatch"
      pcls%w = w
    end if
    if (present(species_id)) then
      if (size(species_id) /= n) error stop "species_id size mismatch"
      pcls%species_id = species_id
    end if
    if (present(source_element)) then
      if (size(source_element) /= n) error stop "source_element size mismatch"
      pcls%source_element = source_element
    end if
  end subroutine init_particles

  !> 既存 SoA の末尾へ粒子群を追加し、既存の alive 状態を保持する。
  !! 容量が足りなければ `extend_particles` が 1.5 倍以上に広げるので、配列長は `pcls%n` より長くなり得る。
  subroutine append_particles(pcls, x, v, q, m, w, species_id, source_element)
    type(particles_soa), intent(inout) :: pcls
    real(dp), intent(in) :: x(:, :), v(:, :), q(:), m(:), w(:)
    integer(i32), intent(in) :: species_id(:)
    integer(i32), intent(in), optional :: source_element(:)
    integer(i32) :: old_count, added_count

    old_count = pcls%n
    added_count = int(size(q), i32)
    if (size(x, 1) /= 3 .or. size(v, 1) /= 3) error stop 'particle append first dimension must be 3'
    if (size(x, 2) /= added_count .or. size(v, 2) /= added_count .or. &
        size(m) /= added_count .or. size(w) /= added_count .or. size(species_id) /= added_count) then
      error stop 'particle append size mismatch'
    end if
    if (present(source_element)) then
      if (size(source_element) /= added_count) error stop 'particle append source_element size mismatch'
    end if
    if (added_count == 0_i32) return
    call extend_particles(pcls, added_count)
    pcls%x(:, old_count + 1:pcls%n) = x
    pcls%v(:, old_count + 1:pcls%n) = v
    pcls%q(old_count + 1:pcls%n) = q
    pcls%m(old_count + 1:pcls%n) = m
    pcls%w(old_count + 1:pcls%n) = w
    pcls%species_id(old_count + 1:pcls%n) = species_id
    if (present(source_element)) pcls%source_element(old_count + 1:pcls%n) = source_element
  end subroutine append_particles

  !> 末尾に `added` 個の枠を足して `pcls%n` を進める。新しい枠は重み 1、粒子種 0、放出元 -1、alive で、
  !! 位置・速度・電荷・質量は呼出し元が埋める。容量が足りなければ必要数と現容量の 1.5 倍の大きい方へ広げ、
  !! 二次電子の世代ごとの追加で毎回全粒子を写し直さない。
  subroutine extend_particles(pcls, added)
    type(particles_soa), intent(inout) :: pcls
    integer(i32), intent(in) :: added
    integer(i32) :: old_count, new_count, capacity

    old_count = pcls%n
    if (added < 0_i32) error stop 'particle append count must be >= 0'
    if (old_count < 0_i32 .or. old_count > huge(0_i32) - added) error stop 'particle append count overflow'
    if (old_count > 0_i32) then
      if (.not. allocated(pcls%x) .or. .not. allocated(pcls%v) .or. .not. allocated(pcls%q) .or. &
          .not. allocated(pcls%m) .or. .not. allocated(pcls%w) .or. .not. allocated(pcls%species_id) .or. &
          .not. allocated(pcls%source_element) .or. .not. allocated(pcls%alive)) then
        error stop 'particle append requires a complete existing SoA'
      end if
    end if
    new_count = old_count + added
    capacity = particle_capacity(pcls)
    if (new_count > capacity) then
      call reserve_particles(pcls, max(new_count, capacity + min(capacity/2_i32, huge(0_i32) - capacity)))
    end if
    pcls%w(old_count + 1:new_count) = 1.0_dp
    pcls%species_id(old_count + 1:new_count) = 0_i32
    pcls%source_element(old_count + 1:new_count) = -1_i32
    pcls%alive(old_count + 1:new_count) = .true.
    pcls%n = new_count
  end subroutine extend_particles

  !> 容量を `capacity` 以上にする。先頭 `pcls%n` 個だけを写し、それより後ろの値は未定義のまま残す。
  subroutine reserve_particles(pcls, capacity)
    type(particles_soa), intent(inout) :: pcls
    integer(i32), intent(in) :: capacity
    real(dp), allocatable :: grown_x(:, :), grown_v(:, :), grown_q(:), grown_m(:), grown_w(:)
    integer(i32), allocatable :: grown_species_id(:), grown_source_element(:)
    logical, allocatable :: grown_alive(:)
    integer(i32) :: n

    if (capacity <= particle_capacity(pcls)) return
    n = pcls%n
    allocate (grown_x(3, capacity), grown_v(3, capacity), grown_q(capacity), grown_m(capacity), &
              grown_w(capacity), grown_species_id(capacity), grown_source_element(capacity), grown_alive(capacity))
    if (n > 0_i32) then
      grown_x(:, :n) = pcls%x(:, :n)
      grown_v(:, :n) = pcls%v(:, :n)
      grown_q(:n) = pcls%q(:n)
      grown_m(:n) = pcls%m(:n)
      grown_w(:n) = pcls%w(:n)
      grown_species_id(:n) = pcls%species_id(:n)
      grown_source_element(:n) = pcls%source_element(:n)
      grown_alive(:n) = pcls%alive(:n)
    end if
    call move_alloc(grown_x, pcls%x)
    call move_alloc(grown_v, pcls%v)
    call move_alloc(grown_q, pcls%q)
    call move_alloc(grown_m, pcls%m)
    call move_alloc(grown_w, pcls%w)
    call move_alloc(grown_species_id, pcls%species_id)
    call move_alloc(grown_source_element, pcls%source_element)
    call move_alloc(grown_alive, pcls%alive)
  end subroutine reserve_particles

  !> 確保済みの粒子枠数。未確保なら 0。
  pure integer(i32) function particle_capacity(pcls) result(capacity)
    type(particles_soa), intent(in) :: pcls

    capacity = 0_i32
    if (allocated(pcls%q)) capacity = int(size(pcls%q), i32)
  end function particle_capacity

  !> 先頭 `n` 個だけを残して配列を詰め直す。容量を先に確保して埋めた SoA を実数に合わせるときに使う。
  subroutine trim_particles(pcls, n)
    type(particles_soa), intent(inout) :: pcls
    integer(i32), intent(in) :: n
    real(dp), allocatable :: kept_x(:, :), kept_v(:, :), kept_q(:), kept_m(:), kept_w(:)
    integer(i32), allocatable :: kept_species_id(:), kept_source_element(:)
    logical, allocatable :: kept_alive(:)

    if (n < 0_i32) error stop 'particle trim count must be >= 0'
    if (.not. allocated(pcls%q)) error stop 'particle trim requires an allocated SoA'
    if (n > size(pcls%q)) error stop 'particle trim count exceeds capacity'
    pcls%n = n
    if (n == size(pcls%q)) return
    allocate (kept_x(3, n), kept_v(3, n), kept_q(n), kept_m(n), kept_w(n), kept_species_id(n), &
              kept_source_element(n), kept_alive(n))
    kept_x = pcls%x(:, :n)
    kept_v = pcls%v(:, :n)
    kept_q = pcls%q(:n)
    kept_m = pcls%m(:n)
    kept_w = pcls%w(:n)
    kept_species_id = pcls%species_id(:n)
    kept_source_element = pcls%source_element(:n)
    kept_alive = pcls%alive(:n)
    call move_alloc(kept_x, pcls%x)
    call move_alloc(kept_v, pcls%v)
    call move_alloc(kept_q, pcls%q)
    call move_alloc(kept_m, pcls%m)
    call move_alloc(kept_w, pcls%w)
    call move_alloc(kept_species_id, pcls%species_id)
    call move_alloc(kept_source_element, pcls%source_element)
    call move_alloc(kept_alive, pcls%alive)
  end subroutine trim_particles

end module bem_particles
