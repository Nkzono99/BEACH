!> 既存の境界イベントと表面放出で使う counter-based 乱数。実行順序や thread ID を使わない。
module bem_counter_rng
  use bem_kinds, only: dp, i32, i64
  implicit none
  private
  public :: counter_hash_mix, counter_uniform
contains
  pure real(dp) function counter_uniform(seed, counter, draw) result(uniform)
    integer(i32), intent(in) :: seed
    integer(i64), intent(in) :: counter(:), draw
    integer(i64) :: state
    integer :: j
    state = counter_hash_mix(104729_i64, int(seed, i64), 32452843_i64)
    do j = 1, size(counter)
      state = counter_hash_mix(state, counter(j), 7919_i64*int(j, i64))
    end do
    state = counter_hash_mix(state, draw, 433494437_i64)
    state = counter_hash_mix(state, draw, 13007_i64)
    uniform = (real(state, dp) - 0.5_dp)/2147483647.0_dp
  end function counter_uniform

  !> 31-bit素数体上の非線形写像でtuple成分をcounter hashへ混合する。
  pure integer(i64) function counter_hash_mix(state, value, salt) result(mixed)
    integer(i64), intent(in) :: state, value, salt
    integer(i64), parameter :: hash_modulus = 2147483647_i64
    integer(i64), parameter :: hash_multiplier = 48271_i64

    mixed = modulo(state*state, hash_modulus)
    mixed = modulo(mixed + modulo(hash_multiplier*state, hash_modulus), hash_modulus)
    mixed = modulo(mixed + modulo(value, hash_modulus), hash_modulus)
    mixed = modulo(mixed + modulo(salt, hash_modulus) + 1_i64, hash_modulus)
    if (mixed == 0_i64) mixed = 1_i64
  end function counter_hash_mix

end module bem_counter_rng
