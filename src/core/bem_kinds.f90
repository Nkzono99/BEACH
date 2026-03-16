!> 倍精度実数と32bit整数のkind定義を集約する基盤モジュール。
module bem_kinds
  use, intrinsic :: iso_fortran_env, only: real64, int32, int64
  implicit none
  integer, parameter :: dp = real64
  integer, parameter :: i32 = int32
  integer, parameter :: i64 = int64
end module bem_kinds
