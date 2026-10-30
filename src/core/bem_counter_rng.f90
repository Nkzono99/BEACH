!> PCG64 による counter-based 乱数。実行順序・rank・thread に依存せず、(seed, counter, draw) だけで値が決まる。
!>
!> counter（粒子 slot、境界 event など任意長の整数組）と seed から Murmur3 で 128 bit の鍵を作り、
!> 上位 64 bit を PCG64 の初期状態、下位 64 bit を stream 番号に使う。draw 番目の値は stream を draw 回進めて取り出す。
!> 鍵の計算は 32 bit の Murmur3 を 4 系統の seed で回すもので、整数演算は 16 bit 分割の乗算で
!> 符号付き 64 bit の範囲に収めている（Fortran の整数 overflow を起こさない）。
module bem_counter_rng
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_int64_t
  use bem_kinds, only: dp, i32, i64
  use m_pcg64_bind, only: pcg64_create, pcg64_destroy, pcg64_random_double53, pcg64_advance
  implicit none
  private

  public :: counter_uniform
  public :: counter_uniforms
  public :: counter_stream_type
  public :: counter_key
  public :: pcg64_open_uniform
  public :: open_uniform_from_double53
  public :: counter_child_id

  integer(i64), parameter :: mask32 = 4294967295_i64
  integer(i64), parameter :: murmur_c1 = int(z'cc9e2d51', i64)
  integer(i64), parameter :: murmur_c2 = int(z'1b873593', i64)
  integer(i64), parameter :: murmur_n = int(z'e6546b64', i64)
  integer(i64), parameter :: fmix_c1 = int(z'85ebca6b', i64)
  integer(i64), parameter :: fmix_c2 = int(z'c2b2ae35', i64)
  !> 4 系統の Murmur3 seed（鍵の 4 つの 32 bit 語）。
  integer(i64), parameter :: lane_seed(4) = [int(z'9e3779b9', i64), int(z'7f4a7c15', i64), &
                                             int(z'f39cc060', i64), int(z'5ced1c2b', i64)]

  !> 1 つの (seed, counter) の draw を連番で引く列。鍵の計算と PCG64 生成器の作成を 1 回で済ませ、
  !! `next_uniform` の n 回目は `counter_uniform(seed, counter, first_draw + n - 1)` と同じ値を返す。
  !! 生成器は C++ 側で確保するので、`open` したら `close` する。
  type :: counter_stream_type
    private
    type(c_ptr) :: handle = c_null_ptr
  contains
    procedure :: open => open_counter_stream
    procedure :: next_uniform => counter_stream_next_uniform
    procedure :: close => close_counter_stream
  end type counter_stream_type

contains

  !> (seed, counter) の draw `first_draw .. first_draw + size(uniforms) - 1` をまとめて返す。
  !! 値は同じ draw の `counter_uniform` と一致する。
  subroutine counter_uniforms(seed, counter, first_draw, uniforms)
    integer(i32), intent(in) :: seed
    integer(i64), intent(in) :: counter(:), first_draw
    real(dp), intent(out) :: uniforms(:)
    type(counter_stream_type) :: stream
    integer :: i

    call stream%open(seed, counter, first_draw)
    do i = 1, size(uniforms)
      uniforms(i) = stream%next_uniform()
    end do
    call stream%close()
  end subroutine counter_uniforms

  !> (seed, counter) の PCG64 を開き、`first_draw` 回進めた位置に置く。
  subroutine open_counter_stream(self, seed, counter, first_draw)
    class(counter_stream_type), intent(inout) :: self
    integer(i32), intent(in) :: seed
    integer(i64), intent(in) :: counter(:), first_draw
    integer(i64) :: state_seed, stream

    if (first_draw < 0_i64) error stop 'counter_uniform draw must be >= 0.'
    call self%close()
    call counter_key(seed, counter, state_seed, stream)
    self%handle = pcg64_create(int(state_seed, c_int64_t), int(stream, c_int64_t))
    if (.not. c_associated(self%handle)) error stop 'failed to create a PCG64 generator.'
    if (first_draw > 0_i64) call pcg64_advance(self%handle, int(first_draw, c_int64_t))
  end subroutine open_counter_stream

  !> 次の draw を (0, 1) で返す。
  real(dp) function counter_stream_next_uniform(self) result(uniform)
    class(counter_stream_type), intent(inout) :: self

    uniform = open_uniform_from_double53(pcg64_random_double53(self%handle))
  end function counter_stream_next_uniform

  !> 生成器を解放する（開いていなければ何もしない）。
  subroutine close_counter_stream(self)
    class(counter_stream_type), intent(inout) :: self

    if (c_associated(self%handle)) call pcg64_destroy(self%handle)
    self%handle = c_null_ptr
  end subroutine close_counter_stream

  !> (seed, counter, draw) で決まる (0, 1) の一様乱数。
  !! @param[in] seed 実行の乱数 seed（`run.rng_seed`）。
  !! @param[in] counter 乱数を使う対象を識別する整数組。
  !! @param[in] draw 同じ対象の中で何番目の値か（0 以上）。
  real(dp) function counter_uniform(seed, counter, draw) result(uniform)
    integer(i32), intent(in) :: seed
    integer(i64), intent(in) :: counter(:), draw
    integer(i64) :: state_seed, stream

    call counter_key(seed, counter, state_seed, stream)
    uniform = pcg64_open_uniform(state_seed, stream, draw)
  end function counter_uniform

  !> seed と counter から PCG64 の (初期状態 seed, stream) を作る。
  pure subroutine counter_key(seed, counter, state_seed, stream)
    integer(i32), intent(in) :: seed
    integer(i64), intent(in) :: counter(:)
    integer(i64), intent(out) :: state_seed, stream
    integer(i64) :: lanes(4)
    integer :: lane

    do lane = 1, 4
      lanes(lane) = murmur3_tuple(lane_seed(lane), int(seed, i64), counter)
    end do
    state_seed = ior(ishft(lanes(1), 32), lanes(2))
    stream = ior(ishft(lanes(3), 32), lanes(4))
  end subroutine counter_key

  !> 親の乱数 id と子の番号から、衝突しにくい 63 bit の子 id を作る（二次電子などの派生粒子用）。
  pure integer(i64) function counter_child_id(parent_id, ordinal) result(child_id)
    integer(i64), intent(in) :: parent_id, ordinal
    integer(i64) :: hi, lo

    hi = murmur3_tuple(lane_seed(1), parent_id, [ordinal, 4099_i64])
    lo = murmur3_tuple(lane_seed(2), parent_id, [ordinal, 4099_i64])
    child_id = iand(ior(ishft(hi, 32), lo), huge(0_i64))
  end function counter_child_id

  !> PCG64(state_seed, stream) を draw 回進めた次の値を (0, 1) で返す（変換は `open_uniform_from_double53`）。
  real(dp) function pcg64_open_uniform(state_seed, stream, draw) result(uniform)
    integer(i64), intent(in) :: state_seed, stream, draw
    type(c_ptr) :: handle

    if (draw < 0_i64) error stop 'counter_uniform draw must be >= 0.'
    handle = pcg64_create(int(state_seed, c_int64_t), int(stream, c_int64_t))
    if (.not. c_associated(handle)) error stop 'failed to create a PCG64 generator.'
    if (draw > 0_i64) call pcg64_advance(handle, int(draw, c_int64_t))
    uniform = open_uniform_from_double53(pcg64_random_double53(handle))
    call pcg64_destroy(handle)
  end function pcg64_open_uniform

  !> 53 bit 格子値 raw = k/2^53（k = 0..2^53-1）を (0, 1) へ写す。
  !! raw に 2^-54 を足し、binary64 の最近接丸めで 1.0 になる最上端 k = 2^53-1 だけを 1 の直前値 1-2^-53 に留める。
  !! それ以外の k では `raw + 2^-54` と同じ binary64 値を返す。
  !! @param[in] raw `pcg64_random_double53` の戻り値（[0, 1) の 53 bit 格子値）。
  !! @return uniform 0 < uniform < 1 の値。
  pure real(dp) function open_uniform_from_double53(raw) result(uniform)
    real(dp), intent(in) :: raw
    real(dp), parameter :: half_ulp = 0.5_dp/9007199254740992.0_dp
    real(dp), parameter :: largest_below_one = nearest(1.0_dp, -1.0_dp)

    uniform = min(raw + half_ulp, largest_below_one)
  end function open_uniform_from_double53

  !> seed と整数組を 32 bit 語の列として Murmur3_x86_32 にかける。
  pure integer(i64) function murmur3_tuple(lane, seed, counter) result(h)
    integer(i64), intent(in) :: lane, seed, counter(:)
    integer :: j
    integer(i64) :: nwords

    h = lane
    call murmur3_word(h, iand(seed, mask32))
    call murmur3_word(h, iand(ishft(seed, -32), mask32))
    do j = 1, size(counter)
      call murmur3_word(h, iand(counter(j), mask32))
      call murmur3_word(h, iand(ishft(counter(j), -32), mask32))
    end do
    nwords = 2_i64*int(size(counter) + 1, i64)
    h = ieor(h, iand(4_i64*nwords, mask32))
    h = fmix32(h)
  end function murmur3_tuple

  pure subroutine murmur3_word(h, word)
    integer(i64), intent(inout) :: h
    integer(i64), intent(in) :: word
    integer(i64) :: k

    k = mul32(word, murmur_c1)
    k = rotl32(k, 15)
    k = mul32(k, murmur_c2)
    h = ieor(h, k)
    h = rotl32(h, 13)
    h = iand(mul32(h, 5_i64) + murmur_n, mask32)
  end subroutine murmur3_word

  pure integer(i64) function fmix32(h0) result(h)
    integer(i64), intent(in) :: h0

    h = h0
    h = ieor(h, ishft(h, -16))
    h = mul32(h, fmix_c1)
    h = ieor(h, ishft(h, -13))
    h = mul32(h, fmix_c2)
    h = ieor(h, ishft(h, -16))
  end function fmix32

  !> 32 bit 語どうしの積 mod 2^32。b を 16 bit に分けて部分積を 2^48 未満に保つ。
  pure integer(i64) function mul32(a, b) result(r)
    integer(i64), intent(in) :: a, b
    integer(i64) :: b_lo, b_hi, t

    b_lo = iand(b, 65535_i64)
    b_hi = iand(ishft(b, -16), 65535_i64)
    t = iand(a*b_hi, 65535_i64)
    r = iand(ishft(t, 16) + a*b_lo, mask32)
  end function mul32

  pure integer(i64) function rotl32(x, s) result(r)
    integer(i64), intent(in) :: x
    integer, intent(in) :: s

    r = iand(ior(ishft(x, s), ishft(x, s - 32)), mask32)
  end function rotl32

end module bem_counter_rng
