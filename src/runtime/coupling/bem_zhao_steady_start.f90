!> Fresh run を Zhao 零電流定常枝の水平平均電荷へ初期化する。
module bem_zhao_steady_start
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_types, only: mesh_type, sim_config
  use bem_physics_config_types, only: periodic2_physics_config, coupling_config
  use bem_outer_plasma_kinetic, only: kinetic_outer_plasma_options_type, &
                                      solve_outer_plasma_zhao_stationary
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, outer_plasma_invalid
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  public :: initialize_zhao_floating_steady_start

contains

  subroutine initialize_zhao_floating_steady_start( &
    mesh, mesh_mode, sim, periodic2, coupling, kinetic_options, state, seed_charge_c, status, message &
    )
    type(mesh_type), intent(inout) :: mesh
    character(len=*), intent(in) :: mesh_mode
    type(sim_config), intent(in) :: sim
    type(periodic2_physics_config), intent(in) :: periodic2
    type(coupling_config), intent(in) :: coupling
    type(kinetic_outer_plasma_options_type), intent(in) :: kinetic_options
    type(outer_plasma_state_type), intent(out) :: state
    real(dp), intent(out) :: seed_charge_c
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    integer(i32) :: element, selected_count, last_selected
    real(dp) :: area_xy, selected_area, seed_factor, z_min, z_max, z_scale, z_tolerance
    real(dp) :: x_min, x_max, y_min, y_max, horizontal_scale, horizontal_tolerance
    real(dp) :: area_relative_error, charge_roundoff

    state = outer_plasma_state_type()
    seed_charge_c = 0.0_dp
    status = outer_plasma_invalid
    message = ''
    if (trim(lower_ascii(mesh_mode)) /= 'template') then
      message = 'Zhao floating steady start currently requires mesh.mode="template"'
      return
    end if
    if (mesh%nelem < 1_i32 .or. .not. allocated(mesh%q_elem) .or. &
        .not. allocated(mesh%elem_mesh_id) .or. .not. allocated(mesh%panel_area) .or. &
        .not. allocated(mesh%normals) .or. .not. allocated(mesh%v0) .or. &
        .not. allocated(mesh%v1) .or. .not. allocated(mesh%v2)) then
      message = 'Zhao floating steady start requires a fully initialized mesh'
      return
    end if
    if (any(mesh%q_elem /= 0.0_dp)) then
      message = 'Zhao floating steady start requires an initially uncharged mesh'
      return
    end if
    area_xy = (sim%box_max(1) - sim%box_min(1))*(sim%box_max(2) - sim%box_min(2))
    if (.not. ieee_is_finite(area_xy) .or. area_xy <= 0.0_dp) then
      message = 'Zhao floating steady start requires a positive finite horizontal box area'
      return
    end if

    selected_count = 0_i32
    last_selected = 0_i32
    selected_area = 0.0_dp
    z_min = huge(1.0_dp)
    z_max = -huge(1.0_dp)
    x_min = huge(1.0_dp)
    x_max = -huge(1.0_dp)
    y_min = huge(1.0_dp)
    y_max = -huge(1.0_dp)
    do element = 1_i32, mesh%nelem
      if (mesh%elem_mesh_id(element) /= coupling%steady_start_mesh_id) cycle
      selected_count = selected_count + 1_i32
      last_selected = element
      if (.not. ieee_is_finite(mesh%panel_area(element)) .or. mesh%panel_area(element) <= 0.0_dp) then
        message = 'Zhao floating steady-start mesh contains a nonpositive panel area'
        return
      end if
      if (.not. all(ieee_is_finite([ &
                                   mesh%v0(:, element), mesh%v1(:, element), mesh%v2(:, element), &
                                   mesh%normals(:, element) &
                                   ]))) then
        message = 'Zhao floating steady-start mesh contains non-finite geometry'
        return
      end if
      selected_area = selected_area + mesh%panel_area(element)
      x_min = min(x_min, mesh%v0(1, element), mesh%v1(1, element), mesh%v2(1, element))
      x_max = max(x_max, mesh%v0(1, element), mesh%v1(1, element), mesh%v2(1, element))
      y_min = min(y_min, mesh%v0(2, element), mesh%v1(2, element), mesh%v2(2, element))
      y_max = max(y_max, mesh%v0(2, element), mesh%v1(2, element), mesh%v2(2, element))
      z_min = min(z_min, mesh%v0(3, element), mesh%v1(3, element), mesh%v2(3, element))
      z_max = max(z_max, mesh%v0(3, element), mesh%v1(3, element), mesh%v2(3, element))
      if (abs(abs(mesh%normals(3, element)) - 1.0_dp) > 1.0e-10_dp) then
        message = 'Zhao floating steady-start mesh must be horizontal'
        return
      end if
    end do
    if (selected_count == 0_i32) then
      message = 'Zhao floating steady-start mesh_id was not found'
      return
    end if
    z_scale = max(1.0_dp, abs(z_min), abs(z_max), abs(sim%box_max(3) - sim%box_min(3)))
    z_tolerance = 256.0_dp*epsilon(1.0_dp)*z_scale
    if (z_max - z_min > z_tolerance) then
      message = 'Zhao floating steady-start mesh must be one coplanar horizontal plane'
      return
    end if
    if (z_max >= sim%box_max(3) - z_tolerance) then
      message = 'Zhao floating steady-start plane must lie below the outer interface'
      return
    end if
    area_relative_error = abs(selected_area - area_xy)/area_xy
    if (.not. ieee_is_finite(selected_area) .or. area_relative_error > 1.0e-8_dp) then
      message = 'Zhao floating steady-start plane area must equal the periodic cell area'
      return
    end if
    horizontal_scale = max( &
                       tiny(1.0_dp), abs(sim%box_min(1)), abs(sim%box_max(1)), &
                       abs(sim%box_min(2)), abs(sim%box_max(2)), &
                       abs(sim%box_max(1) - sim%box_min(1)), &
                       abs(sim%box_max(2) - sim%box_min(2)) &
                       )
    horizontal_tolerance = 1.0e-8_dp*horizontal_scale
    if (abs(x_min - sim%box_min(1)) > horizontal_tolerance .or. &
        abs(x_max - sim%box_max(1)) > horizontal_tolerance .or. &
        abs(y_min - sim%box_min(2)) > horizontal_tolerance .or. &
        abs(y_max - sim%box_max(2)) > horizontal_tolerance) then
      message = 'Zhao floating steady-start plane must span the horizontal periodic cell'
      return
    end if

    call solve_outer_plasma_zhao_stationary(kinetic_options, state, status, message)
    if (status /= outer_plasma_ok) return
    select case (trim(lower_ascii(periodic2%lower_boundary_model)))
    case ('e_bottom_zero')
      seed_factor = 1.0_dp
    case ('symmetric_vacuum')
      seed_factor = 2.0_dp
    case default
      status = outer_plasma_invalid
      state%applicability_status = status
      message = 'Zhao floating steady start received an unsupported lower boundary model'
      return
    end select
    seed_charge_c = seed_factor*eps0*area_xy*state%interface_field
    if (.not. ieee_is_finite(seed_charge_c)) then
      status = outer_plasma_invalid
      state%applicability_status = status
      message = 'Zhao floating steady start produced a non-finite seed charge'
      return
    end if
    do element = 1_i32, mesh%nelem
      if (mesh%elem_mesh_id(element) == coupling%steady_start_mesh_id) then
        mesh%q_elem(element) = seed_charge_c*mesh%panel_area(element)/selected_area
      end if
    end do
    charge_roundoff = seed_charge_c - sum(mesh%q_elem)
    mesh%q_elem(last_selected) = mesh%q_elem(last_selected) + charge_roundoff

    state%interface_z = sim%box_max(3)
    if (allocated(state%z)) state%z = state%z + state%interface_z
    status = outer_plasma_ok
    state%applicability_status = status
    message = 'Zhao floating steady start initialized from the zero-current stationary branch'
  end subroutine initialize_zhao_floating_steady_start

end module bem_zhao_steady_start
