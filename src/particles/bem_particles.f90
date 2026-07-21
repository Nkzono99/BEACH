!> 粒子SoAデータ構造の初期化を提供するモジュール。
module bem_particles
  use bem_kinds, only: dp, i32
  use bem_types, only: particles_soa
  implicit none
contains

  !> 位置・速度・電荷・質量(と任意重み)配列から `particles_soa` を検証付きで構築する。
  !! @param[out] pcls 検証済み配列を内部に保持した `particles_soa` 構造体。
  !! @param[in] x 粒子位置配列 `x(3,n)` [m]。
  !! @param[in] v 粒子速度配列 `v(3,n)` [m/s]。
  !! @param[in] q 粒子電荷配列 `q(n)` [C]。
  !! @param[in] m 粒子質量配列 `m(n)` [kg]。
  !! @param[in] w マクロ粒子重み配列 `w(n)`（省略時は1）。
  !! @param[in] species_id 粒子種ID `species_id(n)`（省略時は0）。
  subroutine init_particles(pcls, x, v, q, m, w, species_id)
    type(particles_soa), intent(out) :: pcls
    real(dp), intent(in) :: x(:, :), v(:, :), q(:), m(:)
    real(dp), intent(in), optional :: w(:)
    integer(i32), intent(in), optional :: species_id(:)
    integer(i32) :: n

    n = size(q)
    if (size(x, 1) /= 3 .or. size(v, 1) /= 3) then
      error stop "particle input first dimension must be 3"
    end if
    if (size(x, 2) /= n .or. size(v, 2) /= n .or. size(m) /= n) then
      error stop "particle input size mismatch"
    end if

    pcls%n = n
    allocate (pcls%x(3, n), pcls%v(3, n), pcls%q(n), pcls%m(n), pcls%w(n), pcls%species_id(n), pcls%alive(n))
    pcls%x = x
    pcls%v = v
    pcls%q = q
    pcls%m = m
    if (present(w)) then
      if (size(w) /= n) error stop "w size mismatch"
      pcls%w = w
    else
      pcls%w = 1.0d0
    end if
    if (present(species_id)) then
      if (size(species_id) /= n) error stop "species_id size mismatch"
      pcls%species_id = species_id
    else
      pcls%species_id = 0_i32
    end if
    pcls%alive = .true.
  end subroutine init_particles

  !> 既存 SoA の末尾へ粒子群を追加し、既存の alive 状態を保持する。
  subroutine append_particles(pcls, x, v, q, m, w, species_id)
    type(particles_soa), intent(inout) :: pcls
    real(dp), intent(in) :: x(:, :), v(:, :), q(:), m(:), w(:)
    integer(i32), intent(in) :: species_id(:)
    real(dp), allocatable :: grown_x(:, :), grown_v(:, :), grown_q(:), grown_m(:), grown_w(:)
    integer(i32), allocatable :: grown_species_id(:)
    logical, allocatable :: grown_alive(:)
    integer(i32) :: old_count, added_count, new_count

    old_count = pcls%n
    added_count = int(size(q), i32)
    if (size(x, 1) /= 3 .or. size(v, 1) /= 3) error stop 'particle append first dimension must be 3'
    if (size(x, 2) /= added_count .or. size(v, 2) /= added_count .or. &
        size(m) /= added_count .or. size(w) /= added_count .or. size(species_id) /= added_count) then
      error stop 'particle append size mismatch'
    end if
    if (added_count == 0_i32) return
    if (old_count < 0_i32 .or. old_count > huge(0_i32) - added_count) then
      error stop 'particle append count overflow'
    end if
    if (old_count > 0_i32) then
      if (.not. allocated(pcls%x) .or. .not. allocated(pcls%v) .or. .not. allocated(pcls%q) .or. &
          .not. allocated(pcls%m) .or. .not. allocated(pcls%w) .or. .not. allocated(pcls%species_id) .or. &
          .not. allocated(pcls%alive)) error stop 'particle append requires a complete existing SoA'
    end if

    new_count = old_count + added_count
    allocate (grown_x(3, new_count), grown_v(3, new_count), grown_q(new_count), grown_m(new_count), &
              grown_w(new_count), grown_species_id(new_count), grown_alive(new_count))
    if (old_count > 0_i32) then
      grown_x(:, :old_count) = pcls%x
      grown_v(:, :old_count) = pcls%v
      grown_q(:old_count) = pcls%q
      grown_m(:old_count) = pcls%m
      grown_w(:old_count) = pcls%w
      grown_species_id(:old_count) = pcls%species_id
      grown_alive(:old_count) = pcls%alive
    end if
    grown_x(:, old_count + 1:new_count) = x
    grown_v(:, old_count + 1:new_count) = v
    grown_q(old_count + 1:new_count) = q
    grown_m(old_count + 1:new_count) = m
    grown_w(old_count + 1:new_count) = w
    grown_species_id(old_count + 1:new_count) = species_id
    grown_alive(old_count + 1:new_count) = .true.

    call move_alloc(grown_x, pcls%x)
    call move_alloc(grown_v, pcls%v)
    call move_alloc(grown_q, pcls%q)
    call move_alloc(grown_m, pcls%m)
    call move_alloc(grown_w, pcls%w)
    call move_alloc(grown_species_id, pcls%species_id)
    call move_alloc(grown_alive, pcls%alive)
    pcls%n = new_count
  end subroutine append_particles

end module bem_particles
