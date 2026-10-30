!> Deterministic geometry-only switch selection and real-charge monitoring.
submodule(bem_electrostatic_snapshot) bem_electrostatic_snapshot_upper_fourier
  use bem_sorting, only: sort_real_ascending, order_by_i32_key
  implicit none
contains

  module procedure init_snapshot_upper_fourier_switch
  type(mesh_type) :: proxy
  real(dp) :: heights(16), mismatch(16), phi_rms, ez_rms, lo, hi
  real(dp), allocatable :: patterns(:, :), sorted(:)
  real(dp) :: median
  integer(i32), allocatable :: order(:)
  integer(i32) :: layers(16), nplane, plane, pattern, count, selected, first, last

  if (self%periodic_options%upper_fourier_switch == 'auto') then
    nplane = 16_i32
    call upper_fourier_auto_height_range(sim, maxval(mesh%bb_max(3, :)), lo, hi)
    do plane = 1, nplane
      heights(plane) = lo + (hi - lo)*real(plane - 1, dp)/real(nplane - 1, dp)
    end do
  else
    nplane = 1_i32
    heights(1) = self%periodic_options%upper_fourier_switch_height
  end if
  do plane = 1, nplane
    layers(plane) = switch_mode_layers(self, mesh, heights(plane))
  end do
  call self%upper_fourier_plan%build(mesh, self%periodic_length(1), self%periodic_length(2), maxval(layers(:nplane)))
  allocate (patterns(mesh%nelem, 3))
  patterns(:, 1) = 1.0_dp
  patterns(:, 2) = merge(1.0_dp, -1.0_dp, modulo(mesh%elem_mesh_id, 2_i32) == 0_i32)
  ! Each mesh uses its centroid-height median; ties take +1 so planar meshes stay normalized.
  ! Ordering the elements by mesh id once visits every mesh as one contiguous run.
  allocate (order(mesh%nelem))
  call order_by_i32_key(mesh%elem_mesh_id, order)
  first = 1_i32
  do while (first <= mesh%nelem)
    last = first
    do while (last < mesh%nelem)
      if (mesh%elem_mesh_id(order(last + 1_i32)) /= mesh%elem_mesh_id(order(first))) exit
      last = last + 1_i32
    end do
    sorted = mesh%centers(3, order(first:last))
    count = last - first + 1_i32
    call sort_real_ascending(sorted)
    median = 0.5_dp*(sorted((count + 1)/2) + sorted((count + 2)/2))
    patterns(order(first:last), 3) = merge(1.0_dp, -1.0_dp, mesh%centers(3, order(first:last)) >= median)
    first = last + 1_i32
  end do
  proxy = mesh
  mismatch = 0.0_dp
  do pattern = 1, 3
    ! Normalize each proxy pattern to 1 pC of total |q| so the reported mismatch has a familiar scale.
    proxy%q_elem = 1.0e-12_dp*patterns(:, pattern)/sum(abs(patterns(:, pattern)))
    call self%nonzero_solver%refresh(proxy)
    ! The source sums depend only on the charges, so one refresh serves every candidate height.
    call self%upper_fourier_plan%refresh(proxy%q_elem)
    do plane = 1, nplane
      call plane_mismatch(self, proxy, sim, heights(plane), layers(plane), phi_rms, ez_rms)
      mismatch(plane) = max(mismatch(plane), phi_rms)
    end do
  end do
  call self%nonzero_solver%refresh(mesh)
  selected = minloc(mismatch(:nplane), dim=1)
  self%diagnostics%upper_fourier_switch_active = .true.
  self%diagnostics%upper_fourier_switch_height_m = heights(selected)
  self%diagnostics%upper_fourier_switch_mode_layers = layers(selected)
  self%diagnostics%upper_fourier_switch_selected_mismatch_v = mismatch(selected)
  if (self%upper_fourier_plan%mode_layers /= layers(selected)) then
    call self%upper_fourier_plan%build(mesh, self%periodic_length(1), self%periodic_length(2), layers(selected))
  end if
  call self%upper_fourier_plan%refresh(mesh%q_elem)
  end procedure init_snapshot_upper_fourier_switch

  integer(i32) function switch_mode_layers(self, mesh, height) result(layers)
    class(electrostatic_snapshot_type), intent(in) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: height
    real(dp) :: gap

    layers = self%periodic_options%upper_fourier_mode_layers
    if (layers /= 0_i32) return
    gap = height - maxval(mesh%bb_max(3, :))
    layers = self%reference_mode_layers
    do while (layers < 32_i32)
      if (exp(-2.0_dp*pi*real(layers + 1, dp)*gap/maxval(self%periodic_length)) <= 1.0e-6_dp) exit
      layers = layers + 1_i32
    end do
  end function switch_mode_layers

  module procedure measure_snapshot_switch_mismatch
  phi_rms = 0.0_dp
  ez_rms = 0.0_dp
  if (.not. self%diagnostics%upper_fourier_switch_active) return
  call self%upper_fourier_plan%refresh(mesh%q_elem)
  call plane_mismatch(self, mesh, sim, self%diagnostics%upper_fourier_switch_height_m, &
                      self%diagnostics%upper_fourier_switch_mode_layers, phi_rms, ez_rms)
  end procedure measure_snapshot_switch_mismatch

  subroutine plane_mismatch(self, mesh, sim, height, layers, phi_rms, ez_rms)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    real(dp), intent(in) :: height
    integer(i32), intent(in) :: layers
    real(dp), intent(out) :: phi_rms, ez_rms
    real(dp) :: target(3), phi_fmm, phi_fourier, e_fmm(3), e_fourier(3)
    integer(i32) :: ix, iy, status

    ! The caller refreshes upper_fourier_plan with mesh%q_elem; the Fourier side uses those cached source sums.
    phi_rms = 0.0_dp
    ez_rms = 0.0_dp
    target(3) = height
    do iy = 1, 12
      target(2) = self%periodic_origin(2) + (real(iy, dp) - 0.5_dp)*self%periodic_length(2)/12.0_dp
      do ix = 1, 12
        target(1) = self%periodic_origin(1) + (real(ix, dp) - 0.5_dp)*self%periodic_length(1)/12.0_dp
        call self%nonzero_solver%eval_potential(mesh, sim, target, phi_fmm)
        call self%nonzero_solver%eval_e(mesh, target, e_fmm)
        call self%upper_fourier_plan%eval(mesh%q_elem, target, phi_fourier, e_fourier, status, mode_layers=layers, &
                                          use_cache=.true.)
        phi_rms = phi_rms + (phi_fmm - phi_fourier)**2
        ez_rms = ez_rms + (e_fmm(3) - e_fourier(3))**2
      end do
    end do
    phi_rms = sqrt(phi_rms/144.0_dp)
    ez_rms = sqrt(ez_rms/144.0_dp)
  end subroutine plane_mismatch

end submodule bem_electrostatic_snapshot_upper_fourier
