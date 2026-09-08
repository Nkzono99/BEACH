!> メッシュ・電荷・電位・台帳のCSV出力。
submodule(bem_output_writer) bem_output_writer_files
  use bem_types, only: surface_model_insulator, surface_model_conductor, surface_model_dielectric
  use bem_string_utils, only: lower_ascii
  implicit none
contains

  !> species 別の signed charge flux と粒子数を `charge_ledger.csv` に保存する。
  module procedure write_charge_ledger_file
  character(len=1024) :: path
  integer :: u, ios, species_idx

  if (ledger%nspecies < 1_i32 .or. .not. allocated(ledger%injected_from_remote)) then
    error stop 'write_charge_ledger_file requires an initialized ledger.'
  end if
  path = trim(out_dir)//'/charge_ledger.csv'
  open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'Failed to open charge_ledger.csv.'
  write (u, '(a)') &
    'batch,species_idx,injected_from_remote_C,emitted_from_surface_C,absorbed_on_surface_C,'// &
    'escaped_to_infinity_C,discarded_unresolved_C,'// &
    'neutral_return_correction_C,neutral_return_weight_scale,neutral_return_unresolved_fraction,'// &
    'fixed_absorbed_target_charge_C,fixed_absorbed_weight_scale,'// &
    'fixed_emission_target_charge_C,fixed_emission_weight_scale,fixed_current_correction_C,'// &
    'fixed_absorbed_applied_charge_C,fixed_emission_applied_charge_C,'// &
    'fixed_escape_target_charge_C,fixed_escape_applied_charge_C,fixed_escape_correction_C,'// &
    'injected_count,emitted_count,absorbed_count,escaped_count,discarded_unresolved_count'
  ! Keep the three applied-charge columns as file-format aliases for existing
  ! readers; the ledger stores only the identical target values.
  do species_idx = 1, ledger%nspecies
    write (u, '(i0,a,i0,18(a,es24.16),5(a,i0))') &
      ledger%batch_count, ',', species_idx, &
      ',', ledger%injected_from_remote(species_idx), &
      ',', ledger%emitted_from_surface(species_idx), &
      ',', ledger%absorbed_on_surface(species_idx), &
      ',', ledger%escaped_to_infinity(species_idx), &
      ',', ledger%discarded_unresolved(species_idx), &
      ',', ledger%neutral_return_correction(species_idx), &
      ',', ledger%neutral_return_weight_scale(species_idx), &
      ',', ledger%neutral_return_unresolved_fraction(species_idx), &
      ',', ledger%fixed_absorbed_target_charge(species_idx), &
      ',', ledger%fixed_absorbed_weight_scale(species_idx), &
      ',', ledger%fixed_emission_target_charge(species_idx), &
      ',', ledger%fixed_emission_weight_scale(species_idx), &
      ',', ledger%fixed_current_correction(species_idx), &
      ',', ledger%fixed_absorbed_target_charge(species_idx), &
      ',', ledger%fixed_emission_target_charge(species_idx), &
      ',', ledger%fixed_escape_target_charge(species_idx), &
      ',', ledger%fixed_escape_target_charge(species_idx), &
      ',', ledger%fixed_escape_correction(species_idx), &
      ',', ledger%injected_count(species_idx), &
      ',', ledger%emitted_count(species_idx), &
      ',', ledger%absorbed_count(species_idx), &
      ',', ledger%escaped_count(species_idx), &
      ',', ledger%discarded_unresolved_count(species_idx)
  end do
  close (u)
  end procedure write_charge_ledger_file

  !> 要素電荷を `charges.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素電荷を含むメッシュ情報。
  module procedure write_charges_file
  character(len=1024) :: charges_path
  integer :: u, ios, i

  charges_path = trim(out_dir)//'/charges.csv'
  open (newunit=u, file=trim(charges_path), status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'Failed to open charges file.'
  write (u, '(a)') 'elem_idx,charge_C'
  do i = 1, mesh%nelem
    write (u, '(i0,a,es24.16)') i, ',', mesh%q_elem(i)
  end do
  close (u)
  end procedure write_charges_file

  !> 事前計算済み電位を `mesh_potential.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素情報（要素数の検証用）。
  !! @param[in] potential_v 各要素重心での電位 [V]。
  module procedure write_mesh_potential_file
  character(len=1024) :: potential_path
  integer :: u, ios, i

  if (size(potential_v) /= mesh%nelem) error stop 'precomputed mesh potential size mismatch.'

  potential_path = trim(out_dir)//'/mesh_potential.csv'
  open (newunit=u, file=trim(potential_path), status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'Failed to open mesh_potential.csv.'
  write (u, '(a)') 'elem_idx,potential_V'
  do i = 1, mesh%nelem
    write (u, '(i0,a,es24.16)') i, ',', potential_v(i)
  end do
  close (u)
  end procedure write_mesh_potential_file

  !> 三角形メッシュを `mesh_triangles.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 頂点座標と要素電荷を含むメッシュ情報。
  module procedure write_mesh_file
  character(len=1024) :: mesh_path
  integer :: u, ios, i

  mesh_path = trim(out_dir)//'/mesh_triangles.csv'
  open (newunit=u, file=trim(mesh_path), status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'Failed to open mesh file.'
  write (u, '(a)') 'elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id'
  do i = 1, mesh%nelem
    write (u, '(i0,10(a,es24.16),a,i0)') i, ',', mesh%v0(1, i), ',', mesh%v0(2, i), ',', mesh%v0(3, i), &
      ',', mesh%v1(1, i), ',', mesh%v1(2, i), ',', mesh%v1(3, i), &
      ',', mesh%v2(1, i), ',', mesh%v2(2, i), ',', mesh%v2(3, i), ',', &
      mesh%q_elem(i), ',', mesh%elem_mesh_id(i)
  end do
  close (u)
  end procedure write_mesh_file

  !> メッシュ識別情報を `mesh_sources.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素ごとの `mesh_id` を含むメッシュ情報。
  !! @param[in] cfg 元の入力設定。
  module procedure write_mesh_sources_file
  character(len=1024) :: path
  character(len=16) :: mode_key, source_kind, template_kind
  character(len=16) :: surface_model
  real(dp) :: epsilon_r
  logical :: has_obj
  integer :: u, ios, i, mesh_id
  integer(i32) :: elem_count

  path = trim(out_dir)//'/mesh_sources.csv'
  open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'Failed to open mesh_sources.csv.'
  write (u, '(a)') 'mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count'

  mode_key = trim(lower_ascii(cfg%mesh_mode))
  if (mode_key == 'obj') then
    source_kind = 'obj'
  else if (mode_key == 'template') then
    source_kind = 'template'
  else
    inquire (file=trim(cfg%obj_path), exist=has_obj)
    if (has_obj) then
      source_kind = 'obj'
    else
      source_kind = 'template'
    end if
  end if

  if (source_kind == 'obj') then
    surface_model = mesh_surface_model_name(mesh, 1)
    epsilon_r = mesh_epsilon_r(mesh, 1)
    write (u, '(i0,a,a,a,a,a,a,a,es16.8,a,i0)') &
      1, ',', 'obj', ',', 'obj', ',', trim(surface_model), ',', epsilon_r, ',', mesh%nelem
    close (u)
    return
  end if

  mesh_id = 0
  do i = 1, size(cfg%templates)
    if (.not. cfg%templates(i)%enabled) cycle
    mesh_id = mesh_id + 1
    template_kind = trim(lower_ascii(cfg%templates(i)%kind))
    surface_model = mesh_surface_model_name(mesh, mesh_id)
    epsilon_r = mesh_epsilon_r(mesh, mesh_id)
    elem_count = int(count(mesh%elem_mesh_id == mesh_id), kind=i32)
    write (u, '(i0,a,a,a,a,a,a,a,es16.8,a,i0)') &
      mesh_id, ',', 'template', ',', trim(template_kind), ',', trim(surface_model), ',', epsilon_r, ',', elem_count
  end do
  close (u)
  end procedure write_mesh_sources_file

  !> mesh_id に対応する表面モデル名を返す。
  !! 複数モデルが混在している場合は最初の要素のモデル名を代表値として返す。
  function mesh_surface_model_name(mesh, mesh_id) result(name)
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: mesh_id
    character(len=16) :: name
    integer :: i

    name = 'insulator'
    if (.not. allocated(mesh%elem_surface_model)) return
    do i = 1, mesh%nelem
      if (mesh%elem_mesh_id(i) /= mesh_id) cycle
      select case (mesh%elem_surface_model(i))
      case (surface_model_insulator)
        name = 'insulator'
      case (surface_model_conductor)
        name = 'conductor'
      case (surface_model_dielectric)
        name = 'dielectric'
      case default
        name = 'unknown'
      end select
      return
    end do
  end function mesh_surface_model_name

  !> mesh_id に対応する相対誘電率を返す。
  function mesh_epsilon_r(mesh, mesh_id) result(epsilon_r)
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: mesh_id
    real(dp) :: epsilon_r
    integer :: i

    epsilon_r = 1.0d0
    if (.not. allocated(mesh%elem_epsilon_r)) return
    do i = 1, mesh%nelem
      if (mesh%elem_mesh_id(i) /= mesh_id) cycle
      epsilon_r = mesh%elem_epsilon_r(i)
      return
    end do
  end function mesh_epsilon_r

end submodule bem_output_writer_files
