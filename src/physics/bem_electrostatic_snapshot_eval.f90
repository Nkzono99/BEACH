!> `bem_electrostatic_snapshot` の局所電場・電位評価を実装する submodule。
submodule(bem_electrostatic_snapshot) bem_electrostatic_snapshot_eval
  implicit none
contains

  module procedure eval_snapshot_local_e
  real(dp) :: zero_potential, zero_field

  if (self%use_cached_kneq0) then
    call self%nonzero_solver%eval_e(mesh, position, electric_field)
  else if (self%use_panel_spectral_reference) then
    call eval_periodic_nonzero_panel_reference( &
      mesh, position, self%periodic_length(1), self%periodic_length(2), &
      self%reference_mode_layers, self%panel_quadrature_order, zero_potential, electric_field &
      )
  else
    call self%nonzero_solver%eval_e(mesh, position, electric_field)
  end if
  if (self%use_zero_mode) then
    call eval_periodic_zero_mode( &
      self%zero_plan, self%zero_state, position(3), zero_mode_trace_plus, zero_potential, zero_field &
      )
    electric_field(3) = electric_field(3) + zero_field
  end if
  electric_field = electric_field + self%prescribed_e
  end procedure eval_snapshot_local_e

  module procedure eval_snapshot_local_phi
  real(dp) :: zero_potential, zero_field, electric_field_dummy(3)

  if (self%use_cached_kneq0) then
    call self%nonzero_solver%eval_potential(mesh, sim, position, potential)
  else if (self%use_panel_spectral_reference) then
    call eval_periodic_nonzero_panel_reference( &
      mesh, position, self%periodic_length(1), self%periodic_length(2), &
      self%reference_mode_layers, self%panel_quadrature_order, potential, electric_field_dummy &
      )
  else
    call self%nonzero_solver%eval_potential(mesh, sim, position, potential)
  end if
  if (self%use_zero_mode) then
    call eval_periodic_zero_mode( &
      self%zero_plan, self%zero_state, position(3), zero_mode_trace_plus, zero_potential, zero_field &
      )
    potential = potential + zero_potential
  end if
  potential = potential - dot_product(self%prescribed_e, position - self%prescribed_phi_origin)
  end procedure eval_snapshot_local_phi

  module procedure eval_snapshot_local_phi_without_primary_self
  type(panel_geometry_type) :: geometry
  real(dp) :: self_potential, self_field(3)
  integer(i32) :: geometry_status

  if (element < 1_i32 .or. element > mesh%nelem) then
    error stop 'snapshot self exclusion element index is out of range.'
  end if

  call self%eval_local_phi(mesh, sim, mesh%centers(:, element), potential)
  call init_panel_geometry( &
    mesh%v0(:, element), mesh%v1(:, element), mesh%v2(:, element), geometry, geometry_status &
    )
  if (geometry_status /= panel_geometry_ok) then
    error stop 'snapshot self exclusion received invalid panel geometry.'
  end if
  call panel_potential_field( &
    geometry, mesh%q_elem(element), mesh%centers(:, element), panel_side_principal_value, &
    self_potential, self_field &
    )
  potential = potential - self_potential
  end procedure eval_snapshot_local_phi_without_primary_self

  module procedure compute_snapshot_mesh_potential
  integer(i32) :: element

  if (size(potential) /= mesh%nelem) error stop 'snapshot mesh-potential size mismatch.'
  if (self%use_panel_spectral_reference .or. self%use_cached_kneq0) then
    do element = 1, mesh%nelem
      call self%eval_local_phi(mesh, sim, mesh%centers(:, element), potential(element))
    end do
  else
    call self%nonzero_solver%compute_mesh_potential(mesh, sim, potential)
    do element = 1, mesh%nelem
      potential(element) = potential(element) - &
                           dot_product(self%prescribed_e, mesh%centers(:, element) - self%prescribed_phi_origin)
    end do
  end if
  end procedure compute_snapshot_mesh_potential

end submodule bem_electrostatic_snapshot_eval
