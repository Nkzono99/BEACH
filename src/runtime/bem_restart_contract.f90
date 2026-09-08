!> 再開時のスキーマ・メッシュ・モデル整合性の検証。
submodule(bem_restart) bem_restart_contract
  use bem_model_fingerprint, only: model_fingerprint, mesh_fingerprint, species_fingerprint
  use bem_physics_config_types, only: &
    field_physics_config, &
    panel_kernel_config, &
    derive_field_panel_config, &
    validate_active_physics_config, &
    physics_config_ok
  use bem_checkpoint_contract, only: checkpoint_schema_version_current
  implicit none
contains

  !> schema v2 fingerprint を照合し、state mappingを壊す不一致と条件変更を区別する。
  module procedure validate_restart_contract
  integer :: u, ios, pos
  integer(i32) :: schema_version
  character(len=512) :: line
  character(len=64) :: key
  character(len=256) :: value
  character(len=16) :: saved_model, saved_mesh, saved_species
  character(len=16) :: current_model, current_mesh, current_species
  logical :: found_schema, found_model, found_mesh, found_species
  integer(i32) :: physics_status
  character(len=256) :: physics_message
  type(field_physics_config) :: field_config
  type(panel_kernel_config) :: panel_config

  status = restart_contract_ok
  message = ''
  schema_version = -1_i32
  saved_model = ''
  saved_mesh = ''
  saved_species = ''
  found_schema = .false.
  found_model = .false.
  found_mesh = .false.
  found_species = .false.

  open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    status = restart_contract_malformed
    message = 'cannot open summary.txt'
    return
  end if
  do
    read (u, '(A)', iostat=ios) line
    if (ios /= 0) exit
    pos = index(line, '=')
    if (pos <= 0) cycle
    key = trim(adjustl(line(:pos - 1)))
    value = trim(adjustl(line(pos + 1:)))
    select case (trim(key))
    case ('checkpoint_schema_version')
      read (value, *, iostat=ios) schema_version
      if (ios /= 0) then
        close (u)
        status = restart_contract_malformed
        message = 'invalid checkpoint_schema_version'
        return
      end if
      found_schema = .true.
    case ('model_fingerprint')
      saved_model = trim(value)
      found_model = .true.
    case ('mesh_fingerprint')
      saved_mesh = trim(value)
      found_mesh = .true.
    case ('species_fingerprint')
      saved_species = trim(value)
      found_species = .true.
    end select
  end do
  close (u)

  ! Legacy checkpoints predate fingerprints and are accepted only for implemented Phase 0 point-source modes.
  if (.not. found_schema) then
    call derive_field_panel_config(app%sim, field_config, panel_config)
    call validate_active_physics_config( &
      app%sim, field_config, app%periodic2, panel_config, physics_status, physics_message &
      )
    if (physics_status /= physics_config_ok) then
      status = restart_contract_mismatch
      message = 'legacy checkpoint is incompatible with this physics model'
    end if
    return
  end if
  if (schema_version < 2_i32 .or. schema_version > checkpoint_schema_version_current) then
    status = restart_contract_unsupported_schema
    message = 'unsupported checkpoint schema version'
    return
  end if
  if (.not. (found_model .and. found_mesh .and. found_species)) then
    status = restart_contract_malformed
    message = 'schema v2 summary is missing fingerprints'
    return
  end if
  current_model = model_fingerprint(app)
  current_mesh = mesh_fingerprint(mesh)
  current_species = species_fingerprint(app)
  if (saved_mesh /= current_mesh) then
    status = restart_contract_mismatch
    message = 'mesh fingerprint differs; saved element charges cannot be mapped safely'
  else if (saved_model /= current_model .and. saved_species /= current_species) then
    status = restart_contract_configuration_changed
    message = 'model and species fingerprints differ'
  else if (saved_model /= current_model) then
    status = restart_contract_configuration_changed
    message = 'model fingerprint differs'
  else if (saved_species /= current_species) then
    status = restart_contract_configuration_changed
    message = 'species fingerprint differs'
  end if
  end procedure validate_restart_contract

end submodule bem_restart_contract
