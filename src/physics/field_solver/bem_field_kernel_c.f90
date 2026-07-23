!> C ABI wrapper for the simulator-independent Coulomb FMM field kernel.
module bem_field_kernel_c
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, c_f_pointer, c_int, c_loc, c_null_char, &
                                                                            c_ptr, c_null_ptr, c_signed_char
  use bem_constants, only: k_coulomb
  use bem_coulomb_fmm_core, only: build_plan, build_panel_plan, destroy_plan, destroy_state, eval_points, &
                                  eval_direct_points, eval_potential_points, eval_direct_potential_points, &
                                  fmm_options_type, fmm_plan_type, fmm_state_type, update_state
  use bem_kinds, only: dp, i32, i64
  use bem_version, only: beach_build_info
  implicit none
  private

  integer(c_int), parameter, public :: beach_kernel_ok = 0_c_int
  integer(c_int), parameter, public :: beach_kernel_invalid_handle = 1_c_int
  integer(c_int), parameter, public :: beach_kernel_invalid_argument = 2_c_int
  integer(c_int), parameter, public :: beach_kernel_not_ready = 3_c_int
  integer(c_int), parameter, public :: beach_kernel_abi_major = 1_c_int
  integer(c_int), parameter, public :: beach_kernel_abi_minor = 0_c_int
  integer, parameter :: periodic_cache_path_max_bytes = 256
  character(len=*), parameter :: default_periodic_cache_dir = '.beach_cache/periodic2'
  real(dp), parameter :: default_periodic_generation_tolerance = 1.0d-8

  type :: field_kernel_handle
    type(fmm_plan_type) :: plan
    type(fmm_state_type) :: state
    logical :: built = .false.
    logical :: charged = .false.
    character(len=periodic_cache_path_max_bytes) :: periodic_cache_dir = default_periodic_cache_dir
    real(dp) :: periodic_generation_tolerance = default_periodic_generation_tolerance
  end type field_kernel_handle

  public :: beach_kernel_create
  public :: beach_kernel_get_abi_version
  public :: beach_kernel_get_build_info
  public :: beach_kernel_destroy
  public :: beach_kernel_build
  public :: beach_kernel_build_panel
  public :: beach_kernel_set_periodic_cache
  public :: beach_kernel_get_periodic_cache_info
  public :: beach_kernel_update_charges
  public :: beach_kernel_eval_e
  public :: beach_kernel_eval_e_direct
  public :: beach_kernel_eval_phi
  public :: beach_kernel_eval_phi_direct
  public :: beach_kernel_force_on_charges

contains

  integer(c_int) function beach_kernel_get_abi_version(major_ptr, minor_ptr) &
    bind(C, name='beach_kernel_get_abi_version') result(status)
    type(c_ptr), value :: major_ptr, minor_ptr
    integer(c_int), pointer :: major, minor

    if (.not. c_associated(major_ptr) .or. .not. c_associated(minor_ptr)) then
      status = beach_kernel_invalid_argument
      return
    end if
    call c_f_pointer(major_ptr, major)
    call c_f_pointer(minor_ptr, minor)
    major = beach_kernel_abi_major
    minor = beach_kernel_abi_minor
    status = beach_kernel_ok
  end function beach_kernel_get_abi_version

  integer(c_int) function beach_kernel_get_build_info(buffer_ptr, buffer_capacity, length_ptr) &
    bind(C, name='beach_kernel_get_build_info') result(status)
    type(c_ptr), value :: buffer_ptr, length_ptr
    integer(c_int), value :: buffer_capacity
    integer(c_int), pointer :: text_length
    integer :: required_length

    if (.not. c_associated(length_ptr)) then
      status = beach_kernel_invalid_argument
      return
    end if
    call c_f_pointer(length_ptr, text_length)
    required_length = len(beach_build_info)
    text_length = int(required_length, c_int)
    if (.not. c_associated(buffer_ptr) .or. buffer_capacity <= int(required_length, c_int)) then
      status = beach_kernel_invalid_argument
      return
    end if

    call copy_text_to_c_buffer(beach_build_info, required_length, buffer_ptr, buffer_capacity)
    status = beach_kernel_ok
  end function beach_kernel_get_build_info

  integer(c_int) function beach_kernel_create(handle) bind(C, name='beach_kernel_create') result(status)
    type(c_ptr), intent(out) :: handle
    type(field_kernel_handle), pointer :: kernel

    allocate (kernel)
    kernel%built = .false.
    kernel%charged = .false.
    kernel%periodic_cache_dir = default_periodic_cache_dir
    kernel%periodic_generation_tolerance = default_periodic_generation_tolerance
    handle = c_loc(kernel)
    status = beach_kernel_ok
  end function beach_kernel_create

  integer(c_int) function beach_kernel_destroy(handle) bind(C, name='beach_kernel_destroy') result(status)
    type(c_ptr), value :: handle
    type(field_kernel_handle), pointer :: kernel

    if (.not. c_associated(handle)) then
      status = beach_kernel_invalid_handle
      return
    end if

    call c_f_pointer(handle, kernel)
    call destroy_state(kernel%state)
    call destroy_plan(kernel%plan)
    deallocate (kernel)
    status = beach_kernel_ok
  end function beach_kernel_destroy

  integer(c_int) function beach_kernel_set_periodic_cache(handle, path_ptr, path_len, tolerance) &
    bind(C, name='beach_kernel_set_periodic_cache') result(status)
    type(c_ptr), value :: handle
    type(c_ptr), value :: path_ptr
    integer(c_int), value :: path_len
    real(c_double), value :: tolerance
    type(field_kernel_handle), pointer :: kernel
    character(kind=c_char), pointer :: path_bytes(:)
    character(len=periodic_cache_path_max_bytes) :: validated_path
    integer :: i

    status = get_kernel(handle, kernel)
    if (status /= beach_kernel_ok) return
    if (.not. c_associated(path_ptr) .or. path_len <= 0_c_int .or. &
        path_len > int(periodic_cache_path_max_bytes, c_int) .or. &
        .not. ieee_is_finite(tolerance) .or. tolerance <= 0.0_c_double) then
      status = beach_kernel_invalid_argument
      return
    end if

    call c_f_pointer(path_ptr, path_bytes, [int(path_len)])
    if (.not. valid_utf8_bytes(path_bytes) .or. c_byte_value(path_bytes(int(path_len))) == 32) then
      status = beach_kernel_invalid_argument
      return
    end if

    validated_path = ''
    do i = 1, int(path_len)
      validated_path(i:i) = transfer(path_bytes(i), ' ')
    end do
    kernel%periodic_cache_dir = validated_path
    kernel%periodic_generation_tolerance = real(tolerance, dp)
    status = beach_kernel_ok
  end function beach_kernel_set_periodic_cache

  integer(c_int) function beach_kernel_get_periodic_cache_info( &
    handle, hit_ptr, build_count_ptr, fingerprint_ptr, fingerprint_capacity, fingerprint_length_ptr, &
    path_ptr, path_capacity, path_length_ptr &
    ) bind(C, name='beach_kernel_get_periodic_cache_info') result(status)
    type(c_ptr), value :: handle
    type(c_ptr), value :: hit_ptr, build_count_ptr, fingerprint_ptr, fingerprint_length_ptr
    type(c_ptr), value :: path_ptr, path_length_ptr
    integer(c_int), value :: fingerprint_capacity, path_capacity
    type(field_kernel_handle), pointer :: kernel
    integer(c_int), pointer :: hit, build_count, fingerprint_length, path_length
    integer :: required_fingerprint_length, required_path_length
    logical :: cached_plan

    status = get_kernel(handle, kernel)
    if (status /= beach_kernel_ok) return
    if (.not. kernel%built) then
      status = beach_kernel_not_ready
      return
    end if
    if (.not. c_associated(hit_ptr) .or. .not. c_associated(build_count_ptr) .or. &
        .not. c_associated(fingerprint_length_ptr) .or. .not. c_associated(path_length_ptr)) then
      status = beach_kernel_invalid_argument
      return
    end if

    call c_f_pointer(hit_ptr, hit)
    call c_f_pointer(build_count_ptr, build_count)
    call c_f_pointer(fingerprint_length_ptr, fingerprint_length)
    call c_f_pointer(path_length_ptr, path_length)
    cached_plan = trim(kernel%plan%options%periodic_far_correction) == 'cached_kneq0'
    if (cached_plan) then
      required_fingerprint_length = len_trim(kernel%plan%periodic_cache_fingerprint)
      required_path_length = len_trim(kernel%plan%periodic_cache_path)
      hit = merge(1_c_int, 0_c_int, kernel%plan%periodic_cache_hit)
      build_count = int(kernel%plan%periodic_operator_build_count, c_int)
    else
      required_fingerprint_length = 0
      required_path_length = 0
      hit = 0_c_int
      build_count = 0_c_int
    end if
    fingerprint_length = int(required_fingerprint_length, c_int)
    path_length = int(required_path_length, c_int)

    if (.not. c_associated(fingerprint_ptr) .or. .not. c_associated(path_ptr) .or. &
        fingerprint_capacity <= int(required_fingerprint_length, c_int) .or. &
        path_capacity <= int(required_path_length, c_int)) then
      status = beach_kernel_invalid_argument
      return
    end if

    if (cached_plan) then
      call copy_text_to_c_buffer( &
        kernel%plan%periodic_cache_fingerprint, required_fingerprint_length, fingerprint_ptr, fingerprint_capacity &
        )
      call copy_text_to_c_buffer(kernel%plan%periodic_cache_path, required_path_length, path_ptr, path_capacity)
    else
      call copy_text_to_c_buffer('', 0, fingerprint_ptr, fingerprint_capacity)
      call copy_text_to_c_buffer('', 0, path_ptr, path_capacity)
    end if
    status = beach_kernel_ok
  end function beach_kernel_get_periodic_cache_info

  integer(c_int) function beach_kernel_build( &
    handle, nsrc, src_pos_ptr, theta, leaf_max, order, softening, use_periodic2, periodic_axes_ptr, &
    periodic_len_ptr, image_layers, far_correction, ewald_alpha, ewald_layers, box_min_ptr, box_max_ptr &
    ) bind(C, name='beach_kernel_build') result(status)
    type(c_ptr), value :: handle
    integer(c_int), value :: nsrc
    type(c_ptr), value :: src_pos_ptr
    real(c_double), value :: theta
    integer(c_int), value :: leaf_max
    integer(c_int), value :: order
    real(c_double), value :: softening
    integer(c_int), value :: use_periodic2
    type(c_ptr), value :: periodic_axes_ptr
    type(c_ptr), value :: periodic_len_ptr
    integer(c_int), value :: image_layers
    integer(c_int), value :: far_correction
    real(c_double), value :: ewald_alpha
    integer(c_int), value :: ewald_layers
    type(c_ptr), value :: box_min_ptr
    type(c_ptr), value :: box_max_ptr

    type(field_kernel_handle), pointer :: kernel
    real(c_double), pointer :: src_pos(:, :)
    integer(c_int), pointer :: periodic_axes(:)
    real(c_double), pointer :: periodic_len(:)
    real(c_double), pointer :: box_min(:)
    real(c_double), pointer :: box_max(:)
    type(fmm_options_type) :: options

    status = get_kernel(handle, kernel)
    if (status /= beach_kernel_ok) return
    if (nsrc <= 0_c_int .or. .not. c_associated(src_pos_ptr)) then
      status = beach_kernel_invalid_argument
      return
    end if
    if (theta <= 0.0_c_double .or. leaf_max <= 0_c_int .or. order < 0_c_int .or. softening < 0.0_c_double) then
      status = beach_kernel_invalid_argument
      return
    end if

    call c_f_pointer(src_pos_ptr, src_pos, [3, int(nsrc)])
    options%theta = real(theta, dp)
    options%leaf_max = int(leaf_max, i32)
    options%order = int(order, i32)
    options%softening = real(softening, dp)
    options%periodic_cache_dir = kernel%periodic_cache_dir
    options%periodic_generation_tolerance = kernel%periodic_generation_tolerance

    if (use_periodic2 /= 0_c_int) then
      if (.not. c_associated(periodic_axes_ptr) .or. .not. c_associated(periodic_len_ptr) .or. &
          .not. c_associated(box_min_ptr) .or. .not. c_associated(box_max_ptr)) then
        status = beach_kernel_invalid_argument
        return
      end if
      if (image_layers < 0_c_int .or. ewald_layers < 0_c_int .or. ewald_alpha < 0.0_c_double) then
        status = beach_kernel_invalid_argument
        return
      end if
      call c_f_pointer(periodic_axes_ptr, periodic_axes, [2])
      call c_f_pointer(periodic_len_ptr, periodic_len, [2])
      call c_f_pointer(box_min_ptr, box_min, [3])
      call c_f_pointer(box_max_ptr, box_max, [3])
      if (any(periodic_axes < 1_c_int) .or. any(periodic_axes > 3_c_int) .or. &
          periodic_axes(1) == periodic_axes(2) .or. any(periodic_len <= 0.0_c_double) .or. &
          any(box_max <= box_min)) then
        status = beach_kernel_invalid_argument
        return
      end if
      options%use_periodic2 = .true.
      options%periodic_axes = int(periodic_axes, i32)
      options%periodic_len = real(periodic_len, dp)
      options%periodic_image_layers = int(image_layers, i32)
      options%periodic_ewald_alpha = real(ewald_alpha, dp)
      options%periodic_ewald_layers = int(ewald_layers, i32)
      options%target_box_min = real(box_min, dp)
      options%target_box_max = real(box_max, dp)
      select case (far_correction)
      case (0_c_int)
        options%periodic_far_correction = 'auto'
      case (1_c_int)
        options%periodic_far_correction = 'none'
      case (2_c_int)
        status = beach_kernel_invalid_argument
        return
      case (3_c_int)
        if (order < 1_c_int .or. image_layers < 1_c_int .or. ewald_layers < 1_c_int) then
          status = beach_kernel_invalid_argument
          return
        end if
        options%periodic_far_correction = 'cached_kneq0'
      case default
        status = beach_kernel_invalid_argument
        return
      end select
    end if

    if (kernel%charged) call destroy_state(kernel%state)
    if (kernel%built) call destroy_plan(kernel%plan)
    call build_plan(kernel%plan, real(src_pos, dp), options)
    kernel%built = .true.
    kernel%charged = .false.
    status = beach_kernel_ok
  end function beach_kernel_build

  integer(c_int) function beach_kernel_build_panel( &
    handle, nsrc, v0_ptr, v1_ptr, v2_ptr, theta, leaf_max, order, softening, use_periodic2, periodic_axes_ptr, &
    periodic_len_ptr, image_layers, far_correction, ewald_alpha, ewald_layers, box_min_ptr, box_max_ptr &
    ) bind(C, name='beach_kernel_build_panel') result(status)
    type(c_ptr), value :: handle, v0_ptr, v1_ptr, v2_ptr
    integer(c_int), value :: nsrc, leaf_max, order, use_periodic2, image_layers, far_correction, ewald_layers
    real(c_double), value :: theta, softening, ewald_alpha
    type(c_ptr), value :: periodic_axes_ptr, periodic_len_ptr, box_min_ptr, box_max_ptr
    type(field_kernel_handle), pointer :: kernel
    real(c_double), pointer :: v0(:, :), v1(:, :), v2(:, :)
    integer(c_int), pointer :: periodic_axes(:)
    real(c_double), pointer :: periodic_len(:), box_min(:), box_max(:)
    type(fmm_options_type) :: options

    status = get_kernel(handle, kernel)
    if (status /= beach_kernel_ok) return
    if (nsrc <= 0_c_int .or. .not. c_associated(v0_ptr) .or. .not. c_associated(v1_ptr) .or. &
        .not. c_associated(v2_ptr) .or. theta <= 0.0_c_double .or. leaf_max <= 0_c_int .or. order < 0_c_int .or. &
        softening /= 0.0_c_double .or. far_correction == 2_c_int) then
      status = beach_kernel_invalid_argument
      return
    end if

    call c_f_pointer(v0_ptr, v0, [3, int(nsrc)])
    call c_f_pointer(v1_ptr, v1, [3, int(nsrc)])
    call c_f_pointer(v2_ptr, v2, [3, int(nsrc)])
    options%theta = real(theta, dp)
    options%leaf_max = int(leaf_max, i32)
    options%order = int(order, i32)
    options%softening = 0.0_dp
    options%periodic_cache_dir = kernel%periodic_cache_dir
    options%periodic_generation_tolerance = kernel%periodic_generation_tolerance

    if (use_periodic2 /= 0_c_int) then
      if (.not. c_associated(periodic_axes_ptr) .or. .not. c_associated(periodic_len_ptr) .or. &
          .not. c_associated(box_min_ptr) .or. .not. c_associated(box_max_ptr) .or. &
          image_layers < 0_c_int .or. ewald_layers < 0_c_int .or. ewald_alpha < 0.0_c_double) then
        status = beach_kernel_invalid_argument
        return
      end if
      call c_f_pointer(periodic_axes_ptr, periodic_axes, [2])
      call c_f_pointer(periodic_len_ptr, periodic_len, [2])
      call c_f_pointer(box_min_ptr, box_min, [3])
      call c_f_pointer(box_max_ptr, box_max, [3])
      if (any(periodic_axes < 1_c_int) .or. any(periodic_axes > 3_c_int) .or. &
          periodic_axes(1) == periodic_axes(2) .or. any(periodic_len <= 0.0_c_double) .or. &
          any(box_max <= box_min)) then
        status = beach_kernel_invalid_argument
        return
      end if
      options%use_periodic2 = .true.
      options%periodic_axes = int(periodic_axes, i32)
      options%periodic_len = real(periodic_len, dp)
      options%periodic_image_layers = int(image_layers, i32)
      options%periodic_ewald_alpha = real(ewald_alpha, dp)
      options%periodic_ewald_layers = int(ewald_layers, i32)
      options%target_box_min = real(box_min, dp)
      options%target_box_max = real(box_max, dp)
      select case (far_correction)
      case (0_c_int, 1_c_int)
        options%periodic_far_correction = 'none'
      case (3_c_int)
        if (order < 1_c_int .or. image_layers < 1_c_int .or. ewald_layers < 1_c_int) then
          status = beach_kernel_invalid_argument
          return
        end if
        options%periodic_far_correction = 'cached_kneq0'
      case default
        status = beach_kernel_invalid_argument
        return
      end select
    end if

    if (kernel%charged) call destroy_state(kernel%state)
    if (kernel%built) call destroy_plan(kernel%plan)
    call build_panel_plan(kernel%plan, real(v0, dp), real(v1, dp), real(v2, dp), options)
    kernel%built = .true.
    kernel%charged = .false.
    status = beach_kernel_ok
  end function beach_kernel_build_panel

  integer(c_int) function beach_kernel_update_charges(handle, nsrc, src_q_ptr) &
    bind(C, name='beach_kernel_update_charges') result(status)
    type(c_ptr), value :: handle
    integer(c_int), value :: nsrc
    type(c_ptr), value :: src_q_ptr
    type(field_kernel_handle), pointer :: kernel
    real(c_double), pointer :: src_q(:)

    status = get_kernel(handle, kernel)
    if (status /= beach_kernel_ok) return
    if (.not. kernel%built) then
      status = beach_kernel_not_ready
      return
    end if
    if (nsrc /= kernel%plan%nsrc .or. .not. c_associated(src_q_ptr)) then
      status = beach_kernel_invalid_argument
      return
    end if

    call c_f_pointer(src_q_ptr, src_q, [int(nsrc)])
    call update_state(kernel%plan, kernel%state, real(src_q, dp))
    kernel%charged = .true.
    status = beach_kernel_ok
  end function beach_kernel_update_charges

  integer(c_int) function beach_kernel_eval_e(handle, ntarget, target_pos_ptr, e_ptr) &
    bind(C, name='beach_kernel_eval_e') result(status)
    type(c_ptr), value :: handle
    integer(c_int), value :: ntarget
    type(c_ptr), value :: target_pos_ptr
    type(c_ptr), value :: e_ptr
    type(field_kernel_handle), pointer :: kernel
    real(c_double), pointer :: target_pos(:, :)
    real(c_double), pointer :: e(:, :)

    status = require_charged_kernel(handle, kernel)
    if (status /= beach_kernel_ok) return
    if (ntarget < 0_c_int .or. .not. c_associated(target_pos_ptr) .or. .not. c_associated(e_ptr)) then
      status = beach_kernel_invalid_argument
      return
    end if
    if (ntarget == 0_c_int) return

    call c_f_pointer(target_pos_ptr, target_pos, [3, int(ntarget)])
    call c_f_pointer(e_ptr, e, [3, int(ntarget)])
    call eval_points(kernel%plan, kernel%state, target_pos, e)
    e = k_coulomb*e
    status = beach_kernel_ok
  end function beach_kernel_eval_e

  integer(c_int) function beach_kernel_eval_phi(handle, ntarget, target_pos_ptr, phi_ptr) &
    bind(C, name='beach_kernel_eval_phi') result(status)
    type(c_ptr), value :: handle
    integer(c_int), value :: ntarget
    type(c_ptr), value :: target_pos_ptr
    type(c_ptr), value :: phi_ptr
    type(field_kernel_handle), pointer :: kernel
    real(c_double), pointer :: target_pos(:, :)
    real(c_double), pointer :: phi(:)

    status = require_charged_kernel(handle, kernel)
    if (status /= beach_kernel_ok) return
    if (ntarget < 0_c_int .or. .not. c_associated(target_pos_ptr) .or. .not. c_associated(phi_ptr)) then
      status = beach_kernel_invalid_argument
      return
    end if
    if (ntarget == 0_c_int) return

    call c_f_pointer(target_pos_ptr, target_pos, [3, int(ntarget)])
    call c_f_pointer(phi_ptr, phi, [int(ntarget)])
    call eval_potential_points(kernel%plan, kernel%state, target_pos, phi)
    phi = k_coulomb*phi
    status = beach_kernel_ok
  end function beach_kernel_eval_phi

  integer(c_int) function beach_kernel_eval_e_direct(handle, ntarget, target_pos_ptr, e_ptr) &
    bind(C, name='beach_kernel_eval_e_direct') result(status)
    type(c_ptr), value :: handle
    integer(c_int), value :: ntarget
    type(c_ptr), value :: target_pos_ptr
    type(c_ptr), value :: e_ptr
    type(field_kernel_handle), pointer :: kernel
    real(c_double), pointer :: target_pos(:, :)
    real(c_double), pointer :: e(:, :)

    status = require_direct_kernel(handle, kernel)
    if (status /= beach_kernel_ok) return
    if (.not. target_count_is_addressable(ntarget) .or. .not. c_associated(target_pos_ptr) .or. &
        .not. c_associated(e_ptr)) then
      status = beach_kernel_invalid_argument
      return
    end if
    if (ntarget == 0_c_int) return

    call c_f_pointer(target_pos_ptr, target_pos, [3, int(ntarget)])
    if (any(.not. ieee_is_finite(target_pos))) then
      status = beach_kernel_invalid_argument
      return
    end if
    call c_f_pointer(e_ptr, e, [3, int(ntarget)])
    call eval_direct_points(kernel%plan, kernel%state, target_pos, e)
    e = k_coulomb*e
    status = beach_kernel_ok
  end function beach_kernel_eval_e_direct

  integer(c_int) function beach_kernel_eval_phi_direct(handle, ntarget, target_pos_ptr, phi_ptr) &
    bind(C, name='beach_kernel_eval_phi_direct') result(status)
    type(c_ptr), value :: handle
    integer(c_int), value :: ntarget
    type(c_ptr), value :: target_pos_ptr
    type(c_ptr), value :: phi_ptr
    type(field_kernel_handle), pointer :: kernel
    real(c_double), pointer :: target_pos(:, :)
    real(c_double), pointer :: phi(:)

    status = require_direct_kernel(handle, kernel)
    if (status /= beach_kernel_ok) return
    if (.not. target_count_is_addressable(ntarget) .or. .not. c_associated(target_pos_ptr) .or. &
        .not. c_associated(phi_ptr)) then
      status = beach_kernel_invalid_argument
      return
    end if
    if (ntarget == 0_c_int) return

    call c_f_pointer(target_pos_ptr, target_pos, [3, int(ntarget)])
    if (any(.not. ieee_is_finite(target_pos))) then
      status = beach_kernel_invalid_argument
      return
    end if
    call c_f_pointer(phi_ptr, phi, [int(ntarget)])
    call eval_direct_potential_points(kernel%plan, kernel%state, target_pos, phi)
    phi = k_coulomb*phi
    status = beach_kernel_ok
  end function beach_kernel_eval_phi_direct

  integer(c_int) function beach_kernel_force_on_charges( &
    handle, ntarget, target_pos_ptr, target_q_ptr, origin_ptr, force_ptr, torque_ptr &
    ) bind(C, name='beach_kernel_force_on_charges') result(status)
    type(c_ptr), value :: handle
    integer(c_int), value :: ntarget
    type(c_ptr), value :: target_pos_ptr
    type(c_ptr), value :: target_q_ptr
    type(c_ptr), value :: origin_ptr
    type(c_ptr), value :: force_ptr
    type(c_ptr), value :: torque_ptr

    type(field_kernel_handle), pointer :: kernel
    real(c_double), pointer :: target_pos(:, :)
    real(c_double), pointer :: target_q(:)
    real(c_double), pointer :: origin(:)
    real(c_double), pointer :: force(:)
    real(c_double), pointer :: torque(:)
    real(dp), allocatable :: e(:, :)
    real(dp) :: f_i(3), r_rel(3)
    integer(i32) :: i

    status = require_charged_kernel(handle, kernel)
    if (status /= beach_kernel_ok) return
    if (ntarget < 0_c_int .or. .not. c_associated(target_pos_ptr) .or. .not. c_associated(target_q_ptr) .or. &
        .not. c_associated(origin_ptr) .or. .not. c_associated(force_ptr) .or. .not. c_associated(torque_ptr)) then
      status = beach_kernel_invalid_argument
      return
    end if

    call c_f_pointer(force_ptr, force, [3])
    call c_f_pointer(torque_ptr, torque, [3])
    force = 0.0_c_double
    torque = 0.0_c_double
    if (ntarget == 0_c_int) return

    call c_f_pointer(target_pos_ptr, target_pos, [3, int(ntarget)])
    call c_f_pointer(target_q_ptr, target_q, [int(ntarget)])
    call c_f_pointer(origin_ptr, origin, [3])
    allocate (e(3, int(ntarget)))
    call eval_points(kernel%plan, kernel%state, target_pos, e)
    e = k_coulomb*e
    do i = 1_i32, int(ntarget, i32)
      f_i = real(target_q(i), dp)*e(:, i)
      r_rel = target_pos(:, i) - real(origin, dp)
      force = force + real(f_i, c_double)
      torque(1) = torque(1) + real(r_rel(2)*f_i(3) - r_rel(3)*f_i(2), c_double)
      torque(2) = torque(2) + real(r_rel(3)*f_i(1) - r_rel(1)*f_i(3), c_double)
      torque(3) = torque(3) + real(r_rel(1)*f_i(2) - r_rel(2)*f_i(1), c_double)
    end do
    deallocate (e)
    status = beach_kernel_ok
  end function beach_kernel_force_on_charges

  integer(c_int) function get_kernel(handle, kernel) result(status)
    type(c_ptr), intent(in) :: handle
    type(field_kernel_handle), pointer, intent(out) :: kernel

    if (.not. c_associated(handle)) then
      nullify (kernel)
      status = beach_kernel_invalid_handle
      return
    end if

    call c_f_pointer(handle, kernel)
    status = beach_kernel_ok
  end function get_kernel

  integer(c_int) function require_charged_kernel(handle, kernel) result(status)
    type(c_ptr), intent(in) :: handle
    type(field_kernel_handle), pointer, intent(out) :: kernel

    status = get_kernel(handle, kernel)
    if (status /= beach_kernel_ok) return
    if (.not. kernel%built .or. .not. kernel%charged) status = beach_kernel_not_ready
  end function require_charged_kernel

  integer(c_int) function require_direct_kernel(handle, kernel) result(status)
    type(c_ptr), intent(in) :: handle
    type(field_kernel_handle), pointer, intent(out) :: kernel

    status = get_kernel(handle, kernel)
    if (status /= beach_kernel_ok) return
    if (.not. kernel%built) then
      status = beach_kernel_not_ready
    else if (kernel%plan%options%use_periodic2) then
      status = beach_kernel_invalid_argument
    else if (.not. kernel%charged) then
      status = beach_kernel_not_ready
    end if
  end function require_direct_kernel

  pure logical function target_count_is_addressable(count) result(addressable)
    integer(c_int), intent(in) :: count
    integer(i64) :: shape_limit

    shape_limit = min(int(huge(0), i64), int(huge(0_i32), i64))
    addressable = count >= 0_c_int .and. int(count, i64) <= shape_limit/3_i64
  end function target_count_is_addressable

  subroutine copy_text_to_c_buffer(text, text_length, buffer_ptr, buffer_capacity)
    character(len=*), intent(in) :: text
    integer, intent(in) :: text_length
    type(c_ptr), intent(in) :: buffer_ptr
    integer(c_int), intent(in) :: buffer_capacity
    character(kind=c_char), pointer :: buffer(:)
    integer :: i

    call c_f_pointer(buffer_ptr, buffer, [int(buffer_capacity)])
    do i = 1, text_length
      buffer(i) = transfer(text(i:i), c_null_char)
    end do
    buffer(text_length + 1) = c_null_char
  end subroutine copy_text_to_c_buffer

  pure logical function valid_utf8_bytes(bytes) result(valid)
    character(kind=c_char), intent(in) :: bytes(:)
    integer :: i, b1, b2

    valid = .false.
    i = 1
    do while (i <= size(bytes))
      b1 = c_byte_value(bytes(i))
      select case (b1)
      case (1:127)
        i = i + 1
      case (194:223)
        if (i + 1 > size(bytes)) return
        if (.not. utf8_continuation(bytes(i + 1))) return
        i = i + 2
      case (224)
        if (i + 2 > size(bytes)) return
        b2 = c_byte_value(bytes(i + 1))
        if (b2 < 160 .or. b2 > 191 .or. .not. utf8_continuation(bytes(i + 2))) return
        i = i + 3
      case (225:236, 238:239)
        if (i + 2 > size(bytes)) return
        if (.not. utf8_continuation(bytes(i + 1)) .or. .not. utf8_continuation(bytes(i + 2))) return
        i = i + 3
      case (237)
        if (i + 2 > size(bytes)) return
        b2 = c_byte_value(bytes(i + 1))
        if (b2 < 128 .or. b2 > 159 .or. .not. utf8_continuation(bytes(i + 2))) return
        i = i + 3
      case (240)
        if (i + 3 > size(bytes)) return
        b2 = c_byte_value(bytes(i + 1))
        if (b2 < 144 .or. b2 > 191 .or. .not. utf8_continuation(bytes(i + 2)) .or. &
            .not. utf8_continuation(bytes(i + 3))) return
        i = i + 4
      case (241:243)
        if (i + 3 > size(bytes)) return
        if (.not. utf8_continuation(bytes(i + 1)) .or. .not. utf8_continuation(bytes(i + 2)) .or. &
            .not. utf8_continuation(bytes(i + 3))) return
        i = i + 4
      case (244)
        if (i + 3 > size(bytes)) return
        b2 = c_byte_value(bytes(i + 1))
        if (b2 < 128 .or. b2 > 143 .or. .not. utf8_continuation(bytes(i + 2)) .or. &
            .not. utf8_continuation(bytes(i + 3))) return
        i = i + 4
      case default
        return
      end select
    end do
    valid = .true.
  end function valid_utf8_bytes

  pure logical function utf8_continuation(byte) result(valid)
    character(kind=c_char), intent(in) :: byte
    integer :: value

    value = c_byte_value(byte)
    valid = value >= 128 .and. value <= 191
  end function utf8_continuation

  pure integer function c_byte_value(byte) result(value)
    character(kind=c_char), intent(in) :: byte
    integer(c_signed_char) :: signed_value

    signed_value = transfer(byte, signed_value)
    value = int(signed_value)
    if (value < 0) value = value + 256
  end function c_byte_value

end module bem_field_kernel_c
