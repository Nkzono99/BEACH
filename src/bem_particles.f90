module bem_particles
  use bem_kinds, only: dp, i32
  use bem_types, only: particles_soa
  implicit none
contains

  subroutine init_particles(pcls, x, v, q, m, w)
    type(particles_soa), intent(out) :: pcls
    real(dp), allocatable, intent(in) :: x(:, :), v(:, :), q(:), m(:)
    real(dp), allocatable, intent(in), optional :: w(:)
    integer(i32) :: n

    n = size(q)
    if (size(x, 1) /= 3 .or. size(v, 1) /= 3) then
      error stop "particle input first dimension must be 3"
    end if
    if (size(x, 2) /= n .or. size(v, 2) /= n .or. size(m) /= n) then
      error stop "particle input size mismatch"
    end if

    pcls%n = n
    allocate(pcls%x(3, n), pcls%v(3, n), pcls%q(n), pcls%m(n), pcls%w(n), pcls%alive(n))
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
    pcls%alive = .true.
  end subroutine init_particles

end module bem_particles
