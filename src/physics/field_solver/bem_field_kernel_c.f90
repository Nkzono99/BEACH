!> C ABI wrapper for the simulator-independent Coulomb FMM field kernel.
module bem_field_kernel_c
  use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_f_pointer, c_int, c_loc, c_ptr, c_null_ptr
  use bem_constants, only: k_coulomb
  use bem_coulomb_fmm_core, only: build_plan, destroy_plan, destroy_state, eval_points, eval_potential_points, &
                                  fmm_options_type, fmm_plan_type, fmm_state_type, update_state
  use bem_kinds, only: dp, i32
  implicit none
  private

  integer(c_int), parameter, public :: beach_kernel_ok = 0_c_int
  integer(c_int), parameter, public :: beach_kernel_invalid_handle = 1_c_int
  integer(c_int), parameter, public :: beach_kernel_invalid_argument = 2_c_int
  integer(c_int), parameter, public :: beach_kernel_not_ready = 3_c_int

  type :: field_kernel_handle
    type(fmm_plan_type) :: plan
    type(fmm_state_type) :: state
    logical :: built = .false.
    logical :: charged = .false.
  end type field_kernel_handle

  public :: beach_kernel_create
  public :: beach_kernel_destroy
  public :: beach_kernel_build
  public :: beach_kernel_update_charges
  public :: beach_kernel_eval_e
  public :: beach_kernel_eval_phi
  public :: beach_kernel_force_on_charges

contains

  integer(c_int) function beach_kernel_create(handle) bind(C, name='beach_kernel_create') result(status)
    type(c_ptr), intent(out) :: handle
    type(field_kernel_handle), pointer :: kernel

    allocate (kernel)
    kernel%built = .false.
    kernel%charged = .false.
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
        options%periodic_far_correction = 'm2l_root_oracle'
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

end module bem_field_kernel_c
