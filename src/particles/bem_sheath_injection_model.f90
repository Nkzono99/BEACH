!> 互換性維持のためのシース注入ラッパモジュール。
module bem_sheath_injection_model
  use bem_sheath_runtime, only: sheath_injection_context, resolve_sheath_injection_context
  implicit none

  public :: sheath_injection_context
  public :: resolve_sheath_injection_context
end module bem_sheath_injection_model
