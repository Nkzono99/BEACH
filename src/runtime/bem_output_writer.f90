!> 実行サマリ・最終CSV・履歴CSVの出力を担当するモジュール。
module bem_output_writer
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_stats, surface_model_insulator, surface_model_conductor, surface_model_dielectric
  use bem_app_config_types, only: app_config
  use bem_charge_ledger, only: charge_ledger_type
  use bem_checkpoint_contract, only: checkpoint_schema_version_current
  use bem_electrostatic_snapshot, only: electrostatic_diagnostics_type
  use bem_external_boundary_contract, only: external_boundary_contract_type, external_boundary_ok, &
                                            external_inflow_none, external_inflow_scalar_barrier, &
                                            external_inflow_legacy_sheath, external_inflow_linear_profile, &
                                            external_inflow_kinetic_profile, external_open_escape, &
                                            external_open_potential_barrier, external_transport_none, &
                                            external_transport_linear_1d, external_transport_kinetic_1d, &
                                            external_transport_unified_3d, resolve_external_boundary_contract
  use bem_outer_plasma_photoelectron, only: photoelectron_histogram_state_type
  use bem_model_fingerprint, only: model_fingerprint, mesh_fingerprint, species_fingerprint
  use bem_version, only: beach_build_id, beach_source_commit, beach_version, beach_version_mode
  use bem_filesystem, only: create_directories, filesystem_empty_path, filesystem_not_directory, filesystem_os_error, &
                            filesystem_success
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  public :: open_history_writer
  public :: open_potential_history_writer
  public :: print_run_summary
  public :: write_result_files
  public :: ensure_output_dir

contains

  !> 履歴 CSV のオープンとヘッダ初期化を行う。
  !! @param[in] app 出力設定を含むアプリ設定。
  !! @param[in] resumed 再開実行かどうか。
  !! @param[out] history_opened 履歴ファイルを開けた場合に `.true.`。
  !! @param[out] history_unit 履歴CSVの出力ユニット番号（未使用時は `-1`）。
  subroutine open_history_writer(app, resumed, history_opened, history_unit)
    type(app_config), intent(in) :: app
    logical, intent(in) :: resumed
    logical, intent(out) :: history_opened
    integer, intent(out) :: history_unit
    character(len=1024) :: history_path
    integer :: ios
    logical :: history_exists

    history_opened = .false.
    history_unit = -1
    if (.not. app%write_output) return
    if (app%history_stride <= 0) return

    call ensure_output_dir(app%output_dir)

    history_path = trim(app%output_dir)//'/charge_history.csv'
    inquire (file=trim(history_path), exist=history_exists)
    if (resumed) then
      open (newunit=history_unit, file=trim(history_path), status='unknown', position='append', action='write', iostat=ios)
    else
      open (newunit=history_unit, file=trim(history_path), status='replace', action='write', iostat=ios)
    end if
    if (ios /= 0) error stop 'Failed to open charge history file.'

    if (.not. resumed .or. .not. history_exists) then
      ! 再開時は既存ファイルがある場合だけヘッダ追記を避ける。
      write (history_unit, '(a)') 'batch,processed_particles,rel_change,elem_idx,charge_C'
    end if
    history_opened = .true.
  end subroutine open_history_writer

  !> 電位履歴 CSV のオープンとヘッダ初期化を行う。
  !! @param[in] app 出力設定を含むアプリ設定。
  !! @param[in] resumed 再開実行かどうか。
  !! @param[out] potential_history_opened ファイルを開けた場合に `.true.`。
  !! @param[out] potential_history_unit 電位履歴CSVの出力ユニット番号（未使用時は `-1`）。
  subroutine open_potential_history_writer(app, resumed, potential_history_opened, potential_history_unit)
    type(app_config), intent(in) :: app
    logical, intent(in) :: resumed
    logical, intent(out) :: potential_history_opened
    integer, intent(out) :: potential_history_unit
    character(len=1024) :: path
    integer :: ios
    logical :: file_exists

    potential_history_opened = .false.
    potential_history_unit = -1
    if (.not. app%write_output) return
    if (.not. app%write_potential_history) return
    if (app%history_stride <= 0) return

    call ensure_output_dir(app%output_dir)

    path = trim(app%output_dir)//'/potential_history.csv'
    inquire (file=trim(path), exist=file_exists)
    if (resumed) then
      open (newunit=potential_history_unit, file=trim(path), status='unknown', position='append', action='write', iostat=ios)
    else
      open (newunit=potential_history_unit, file=trim(path), status='replace', action='write', iostat=ios)
    end if
    if (ios /= 0) error stop 'Failed to open potential history file.'

    if (.not. resumed .or. .not. file_exists) then
      write (potential_history_unit, '(a)') 'batch,elem_idx,potential_V'
    end if
    potential_history_opened = .true.
  end subroutine open_potential_history_writer

  !> 実行結果の主要統計を標準出力へ表示する。
  !! @param[in] mesh 実行後のメッシュ情報。
  !! @param[in] stats 実行後の統計値。
  subroutine print_run_summary(mesh, stats)
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    integer(i32) :: dielectric_count

    print '(a,i0)', 'mesh nelem=', mesh%nelem
    print '(a,i0)', 'processed_particles=', stats%processed_particles
    print '(a,i0)', 'absorbed=', stats%absorbed
    print '(a,i0)', 'escaped=', stats%escaped
    print '(a,i0)', 'batches=', stats%batches
    print '(a,i0)', 'escaped_boundary=', stats%escaped_boundary
    print '(a,i0)', 'survived_max_step=', stats%survived_max_step
    print '(a,i0)', 'multiple_box_events_soft_discarded=', stats%multiple_box_events_soft_discarded
    print '(a,es12.4)', 'multiple_box_events_soft_discarded_abs_charge_C=', &
      stats%multiple_box_events_soft_discarded_abs_charge
    print '(a,es12.4)', 'last_rel_change=', stats%last_rel_change
    print '(a,*(es12.4,1x))', 'mesh charges=', mesh%q_elem
    dielectric_count = count_dielectric_surfaces(mesh)
    if (dielectric_count > 0_i32) then
      print '(a,i0)', 'surface_model_dielectric_elem_count=', dielectric_count
      print '(a)', 'surface_model_note=dielectric surface models are metadata-only in this version.'
    end if
  end subroutine print_run_summary

  !> 解析結果を `summary.txt` / `charges.csv` / `mesh_triangles.csv` などとして保存する。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 書き出し対象のメッシュ。
  !! @param[in] stats 書き出し対象の統計値。
  !! @param[in] cfg 出力設定を含むアプリ設定。
  subroutine write_result_files( &
    out_dir, mesh, stats, cfg, mpi_world_size, mesh_potential_v, charge_ledger, electrostatic_diagnostics, &
    photoelectron_state &
    )
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    type(app_config), intent(in) :: cfg
    integer(i32), intent(in), optional :: mpi_world_size
    real(dp), intent(in), optional :: mesh_potential_v(:)
    type(charge_ledger_type), intent(in), optional :: charge_ledger
    type(electrostatic_diagnostics_type), intent(in), optional :: electrostatic_diagnostics
    type(photoelectron_histogram_state_type), intent(in), optional :: photoelectron_state

    if (cfg%outer_plasma%photoelectron_histogram_enabled) then
      if (.not. present(photoelectron_state)) then
        error stop 'photoelectron_histogram_enabled=true requires histogram state for output.'
      end if
      if (.not. photoelectron_state%ready) then
        error stop 'photoelectron_histogram_enabled=true requires ready histogram state for output.'
      end if
    end if
    call ensure_output_dir(out_dir)
    call write_summary_file( &
      out_dir, mesh, stats, cfg, mpi_world_size=mpi_world_size, charge_ledger=charge_ledger, &
      electrostatic_diagnostics=electrostatic_diagnostics, photoelectron_state=photoelectron_state &
      )
    call write_charges_file(out_dir, mesh)
    if (cfg%write_mesh_potential) then
      if (.not. present(mesh_potential_v)) then
        error stop 'write_result_files: mesh_potential_v is required when write_mesh_potential is enabled.'
      end if
      call write_mesh_potential_file(out_dir, mesh, mesh_potential_v)
    end if
    call write_mesh_file(out_dir, mesh)
    call write_mesh_sources_file(out_dir, mesh, cfg)
    if (present(charge_ledger)) call write_charge_ledger_file(out_dir, charge_ledger)
    if (cfg%outer_plasma%photoelectron_histogram_enabled) then
      call write_photoelectron_histogram_file(out_dir, photoelectron_state)
    end if
    if (present(electrostatic_diagnostics)) then
      if (allocated(electrostatic_diagnostics%outer_profile_z)) then
        call write_outer_plasma_profile_file(out_dir, electrostatic_diagnostics)
      end if
    end if
  end subroutine write_result_files

  subroutine write_outer_plasma_profile_file(out_dir, diagnostics)
    character(len=*), intent(in) :: out_dir
    type(electrostatic_diagnostics_type), intent(in) :: diagnostics
    character(len=1024) :: path
    integer :: u, ios
    integer(i32) :: point

    if (.not. allocated(diagnostics%outer_profile_potential) .or. &
        .not. allocated(diagnostics%outer_profile_field) .or. &
        .not. allocated(diagnostics%outer_profile_charge_density) .or. &
        size(diagnostics%outer_profile_z) /= size(diagnostics%outer_profile_potential) .or. &
        size(diagnostics%outer_profile_z) /= size(diagnostics%outer_profile_field) .or. &
        size(diagnostics%outer_profile_z) /= size(diagnostics%outer_profile_charge_density)) then
      error stop 'Outer-plasma diagnostic profile is incomplete.'
    end if
    path = trim(out_dir)//'/outer_plasma_profile.csv'
    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open outer_plasma_profile.csv.'
    write (u, '(a)') 'point,z_m,potential_V,field_V_m,charge_density_C_m3'
    do point = 1_i32, size(diagnostics%outer_profile_z)
      write (u, '(i0,4(a,es24.16))') point, ',', diagnostics%outer_profile_z(point), ',', &
        diagnostics%outer_profile_potential(point), ',', diagnostics%outer_profile_field(point), ',', &
        diagnostics%outer_profile_charge_density(point)
    end do
    close (u)
  end subroutine write_outer_plasma_profile_file

  subroutine write_photoelectron_histogram_file(out_dir, state)
    character(len=*), intent(in) :: out_dir
    type(photoelectron_histogram_state_type), intent(in) :: state
    character(len=1024) :: path
    integer :: u, ios
    integer(i32) :: bin
    real(dp) :: energy_low, energy_high

    path = trim(out_dir)//'/photoelectron_histogram.csv'
    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open photoelectron_histogram.csv.'
    write (u, '(a)') &
      'bin,energy_low_J,energy_high_J,signed_charge_C,kinetic_energy_J,tangential_momentum_x_kg_m_s,'// &
      'tangential_momentum_y_kg_m_s,count,previous_signed_charge_C,previous_kinetic_energy_J,'// &
      'previous_tangential_momentum_x_kg_m_s,previous_tangential_momentum_y_kg_m_s,previous_count'
    do bin = 1_i32, state%cumulative%nbins
      energy_low = real(bin - 1_i32, dp)*state%cumulative%energy_max/real(state%cumulative%nbins, dp)
      energy_high = real(bin, dp)*state%cumulative%energy_max/real(state%cumulative%nbins, dp)
      write (u, '(i0,6(a,es24.16),a,i0,4(a,es24.16),a,i0)') &
        bin, ',', energy_low, ',', energy_high, ',', state%cumulative%signed_charge(bin), ',', &
        state%cumulative%kinetic_energy(bin), ',', state%cumulative%tangential_momentum_x(bin), ',', &
        state%cumulative%tangential_momentum_y(bin), ',', state%cumulative%count(bin), ',', &
        state%previous_batch%signed_charge(bin), ',', state%previous_batch%kinetic_energy(bin), ',', &
        state%previous_batch%tangential_momentum_x(bin), ',', state%previous_batch%tangential_momentum_y(bin), ',', &
        state%previous_batch%count(bin)
    end do
    close (u)
  end subroutine write_photoelectron_histogram_file

  !> 出力ディレクトリを作成する。
  !! @param[in] out_dir 作成対象ディレクトリのパス。
  subroutine ensure_output_dir(out_dir)
    character(len=*), intent(in) :: out_dir
    integer :: status

    call create_directories(out_dir, status)
    select case (status)
    case (filesystem_success)
      return
    case (filesystem_empty_path)
      error stop 'Failed to create output directory: output path is empty.'
    case (filesystem_not_directory)
      error stop 'Failed to create output directory: a path component is not an accessible directory.'
    case (filesystem_os_error)
      error stop 'Failed to create output directory: operating-system directory creation failed.'
    case default
      error stop 'Failed to create output directory: unexpected filesystem status.'
    end select
  end subroutine ensure_output_dir

  !> 実行統計を `summary.txt` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh メッシュ情報（要素数を書き出す）。
  !! @param[in] stats 実行統計。
  subroutine write_summary_file( &
    out_dir, mesh, stats, cfg, mpi_world_size, charge_ledger, electrostatic_diagnostics, photoelectron_state &
    )
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    type(app_config), intent(in) :: cfg
    integer(i32), intent(in), optional :: mpi_world_size
    type(charge_ledger_type), intent(in), optional :: charge_ledger
    type(electrostatic_diagnostics_type), intent(in), optional :: electrostatic_diagnostics
    type(photoelectron_histogram_state_type), intent(in), optional :: photoelectron_state
    type(external_boundary_contract_type) :: resolved_boundary
    character(len=1024) :: summary_path
    character(len=256) :: boundary_message
    integer :: u, ios
    integer(i32) :: world_size, boundary_status

    if (cfg%outer_plasma%photoelectron_histogram_enabled) then
      if (.not. present(photoelectron_state)) then
        error stop 'photoelectron_histogram_enabled=true requires histogram state for summary output.'
      end if
      if (.not. photoelectron_state%ready) then
        error stop 'photoelectron_histogram_enabled=true requires ready histogram state for summary output.'
      end if
    end if
    call resolve_external_boundary_contract( &
      cfg%sim%reservoir_potential_model, cfg%sim%sheath_injection_model, cfg%sim%open_boundary_model, &
      cfg%outer_plasma%model, cfg%outer_plasma%kinetic_closure, cfg%outer_plasma%return_model, &
      cfg%coupling%particle_transfer_mode, cfg%coupling%outer_queue_enabled, resolved_boundary, &
      boundary_status, boundary_message &
      )
    if (boundary_status /= external_boundary_ok) then
      error stop 'write_summary_file: invalid external boundary contract: '//trim(boundary_message)
    end if
    summary_path = trim(out_dir)//'/summary.txt'
    open (newunit=u, file=trim(summary_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open summary file.'
    world_size = 1_i32
    if (present(mpi_world_size)) world_size = max(1_i32, mpi_world_size)
    write (u, '(a,i0)') 'checkpoint_schema_version=', checkpoint_schema_version_current
    write (u, '(a)') 'build_info_schema_version=1'
    write (u, '(a,a)') 'build_version=', beach_version
    write (u, '(a,a)') 'build_version_mode=', beach_version_mode
    write (u, '(a,a)') 'build_source_commit=', beach_source_commit
    write (u, '(a,a)') 'build_id=', beach_build_id
    write (u, '(a,a)') 'model_fingerprint=', model_fingerprint(cfg)
    write (u, '(a,a)') 'mesh_fingerprint=', mesh_fingerprint(mesh)
    write (u, '(a,a)') 'species_fingerprint=', species_fingerprint(cfg)
    write (u, '(a,i0)') 'mesh_nelem=', mesh%nelem
    write (u, '(a,i0)') 'mesh_count=', max(1_i32, maxval(mesh%elem_mesh_id))
    write (u, '(a,i0)') 'mpi_world_size=', world_size
    write (u, '(a,i0)') 'processed_particles=', stats%processed_particles
    write (u, '(a,i0)') 'absorbed=', stats%absorbed
    write (u, '(a,i0)') 'escaped=', stats%escaped
    write (u, '(a,i0)') 'batches=', stats%batches
    write (u, '(a,i0)') 'escaped_boundary=', stats%escaped_boundary
    write (u, '(a,i0)') 'survived_max_step=', stats%survived_max_step
    write (u, '(a,i0)') 'multiple_box_events_soft_discarded=', stats%multiple_box_events_soft_discarded
    write (u, '(a,es24.16)') 'multiple_box_events_soft_discarded_abs_charge_C=', &
      stats%multiple_box_events_soft_discarded_abs_charge
    write (u, '(a,es24.16)') 'last_rel_change=', stats%last_rel_change
    write (u, '(a)') 'particle_time_centering=same_time_midpoint_boris'
    write (u, '(a,a)') 'field_backend=', trim(cfg%field%backend)
    write (u, '(a,a)') 'field_normalization=', trim(cfg%field%normalization)
    write (u, '(a,a)') 'field_source_model=', trim(cfg%panel%source_model)
    write (u, '(a,a)') 'field_kernel_id=', trim(cfg%panel%kernel_id)
    write (u, '(a,a)') 'periodic2_nonzero_mode_backend=', trim(cfg%periodic2%nonzero_mode_backend)
    write (u, '(a,a)') 'periodic2_zero_mode_policy=', trim(cfg%periodic2%zero_mode_policy)
    write (u, '(a,a)') 'periodic2_lower_boundary_model=', trim(cfg%periodic2%lower_boundary_model)
    write (u, '(a,a)') 'periodic2_cache_dir=', trim(cfg%sim%field_periodic_cache_dir)
    write (u, '(a,es24.16)') 'periodic2_generation_tolerance=', &
      cfg%sim%field_periodic_generation_tolerance
    write (u, '(a,a)') 'outer_plasma_model=', trim(cfg%outer_plasma%model)
    write (u, '(a,a)') 'outer_plasma_kinetic_closure=', trim(cfg%outer_plasma%kinetic_closure)
    write (u, '(a,a)') 'outer_plasma_zhao_branch_requested=', trim(cfg%outer_plasma%zhao_branch)
    write (u, '(a,es24.16)') 'outer_plasma_photoelectron_source_scale=', &
      cfg%outer_plasma%photoelectron_source_scale
    write (u, '(a,a)') 'coupling_update_mode=', trim(cfg%coupling%update_mode)
    write (u, '(a,a)') 'coupling_particle_transfer_mode=', trim(cfg%coupling%particle_transfer_mode)
    write (u, '(a,a)') 'coupling_steady_start_mode=', trim(cfg%coupling%steady_start_mode)
    write (u, '(a,i0)') 'coupling_steady_start_mesh_id=', cfg%coupling%steady_start_mesh_id
    write (u, '(a,l1)') 'coupling_outer_queue_enabled=', cfg%coupling%outer_queue_enabled
    write (u, '(a,a)') 'external_inflow_map=', trim(external_inflow_map_name(resolved_boundary%inflow_map))
    write (u, '(a,a)') 'external_ordinary_open_model=', &
      trim(external_open_model_name(resolved_boundary%ordinary_open_model))
    write (u, '(a,a)') 'external_interface_transport=', &
      trim(external_transport_name(resolved_boundary%interface_transport))
    write (u, '(a,a)') 'outer_particle_mode_resolved=', trim(outer_particle_mode_name(resolved_boundary))
    if (present(electrostatic_diagnostics)) then
      write (u, '(a,l1)') 'electrostatic_split_periodic_active=', electrostatic_diagnostics%split_periodic_active
      write (u, '(a,l1)') 'electrostatic_applicable=', electrostatic_diagnostics%applicable
      write (u, '(a,a)') 'electrostatic_status=', trim(electrostatic_diagnostics%status)
      write (u, '(a,i0)') 'interface_sample_n=', electrostatic_diagnostics%interface_sample_n
      write (u, '(a,es24.16)') 'interface_potential_V=', electrostatic_diagnostics%interface_potential
      write (u, '(a,es24.16)') 'interface_normal_field_V_m=', electrostatic_diagnostics%interface_field
      write (u, '(a,es24.16)') 'interface_eta_phi_kneq0=', electrostatic_diagnostics%eta_phi_kneq0
      write (u, '(a,es24.16)') 'interface_eta_field_kneq0=', electrostatic_diagnostics%eta_field_kneq0
      write (u, '(a,es24.16)') 'interface_eta_gap=', electrostatic_diagnostics%eta_gap
      write (u, '(a,es24.16)') 'interface_eta_local_charge=', electrostatic_diagnostics%eta_local_charge
      write (u, '(a,es24.16)') 'gauss_residual_C=', electrostatic_diagnostics%gauss_residual
      write (u, '(a,es24.16)') 'outer_integrated_charge_C=', electrostatic_diagnostics%outer_integrated_charge
      write (u, '(a,i0)') 'outer_nonlinear_iterations=', electrostatic_diagnostics%outer_nonlinear_iterations
      write (u, '(a,es24.16)') 'outer_nonlinear_residual=', electrostatic_diagnostics%outer_nonlinear_residual
      write (u, '(a,i0)') 'outer_applicability_status=', electrostatic_diagnostics%outer_applicability_status
      write (u, '(a,es24.16)') 'outer_infinity_potential_V=', &
        electrostatic_diagnostics%outer_infinity_potential
      write (u, '(a,es24.16)') 'outer_debye_length_m=', electrostatic_diagnostics%outer_debye_length
      write (u, '(a,es24.16)') 'outer_linearity_ratio=', electrostatic_diagnostics%outer_linearity_ratio
      write (u, '(a,es24.16)') 'outer_max_linearity_ratio=', &
        electrostatic_diagnostics%outer_max_linearity_ratio
      write (u, '(a,es24.16)') 'outer_integrated_charge_per_area_C_m2=', &
        electrostatic_diagnostics%outer_integrated_charge_per_area
      write (u, '(a,es24.16)') 'outer_electron_current_density_A_m2=', &
        electrostatic_diagnostics%outer_electron_current_density
      write (u, '(a,es24.16)') 'outer_ion_current_density_A_m2=', &
        electrostatic_diagnostics%outer_ion_current_density
      write (u, '(a,es24.16)') 'outer_photoelectron_current_density_A_m2=', &
        electrostatic_diagnostics%outer_photoelectron_current_density
      write (u, '(a,es24.16)') 'outer_total_current_density_A_m2=', &
        electrostatic_diagnostics%outer_total_current_density
      if (trim(electrostatic_diagnostics%outer_kinetic_closure) == 'zhao_charge_driven') then
        write (u, '(a,a)') 'outer_plasma_zhao_branch_resolved=', electrostatic_diagnostics%outer_zhao_branch
        write (u, '(a,es24.16)') 'outer_plasma_zhao_phi0_V=', electrostatic_diagnostics%outer_zhao_phi0
        write (u, '(a,es24.16)') 'outer_plasma_zhao_phi_minimum_V=', &
          electrostatic_diagnostics%outer_zhao_phi_minimum
        write (u, '(a,es24.16)') 'outer_plasma_zhao_electron_density_infinity_m3=', &
          electrostatic_diagnostics%outer_zhao_electron_density_infinity
        if (cfg%coupling%outer_queue_enabled) then
          write (u, '(a,es24.16)') 'outer_photoelectron_population_fraction=', &
            electrostatic_diagnostics%outer_photoelectron_population_fraction
          write (u, '(a,es24.16)') 'outer_photoelectron_column_per_area_m2=', &
            electrostatic_diagnostics%outer_photoelectron_column_per_area
          write (u, '(a,es24.16)') 'outer_photoelectron_column_target_per_area_m2=', &
            electrostatic_diagnostics%outer_photoelectron_column_target_per_area
          write (u, '(a,es24.16)') 'outer_photoelectron_column_residual_per_area_m2=', &
            electrostatic_diagnostics%outer_photoelectron_column_residual_per_area
          write (u, '(a,i0)') 'outer_queue_event_count=', electrostatic_diagnostics%outer_queue_event_count
          write (u, '(a,es24.16)') 'outer_queue_signed_charge_C=', &
            electrostatic_diagnostics%outer_queue_signed_charge
          write (u, '(a,a)') 'outer_queue_fingerprint=', &
            trim(electrostatic_diagnostics%outer_queue_fingerprint)
        end if
      end if
      write (u, '(a,es24.16)') 'outer_accessible_fraction_min=', &
        electrostatic_diagnostics%accessible_fraction_min
      write (u, '(a,es24.16)') 'outer_accessible_fraction_max=', &
        electrostatic_diagnostics%accessible_fraction_max
      write (u, '(a,es24.16)') 'outer_accessible_fraction_refinement_error=', &
        electrostatic_diagnostics%accessible_fraction_refinement_error
      write (u, '(a,es24.16)') 'outer_nonzero_tail_linearity=', &
        electrostatic_diagnostics%nonzero_tail_linearity
      write (u, '(a,es24.16)') 'outer_response_start_z_m=', &
        electrostatic_diagnostics%response_start_z
      write (u, '(a,i0)') 'last_outer_update_batch=', electrostatic_diagnostics%last_outer_update_batch
      write (u, '(a,es24.16)') 'max_outer_flight_time_s=', electrostatic_diagnostics%max_outer_flight_time
      write (u, '(a,es24.16)') 'max_outer_frozen_field_ratio=', electrostatic_diagnostics%max_frozen_field_ratio
      write (u, '(a,es24.16)') 'max_outer_energy_relative_error=', &
        electrostatic_diagnostics%max_outer_energy_relative_error
      write (u, '(a,l1)') 'periodic2_cache_hit=', electrostatic_diagnostics%periodic_cache_hit
      write (u, '(a,i0)') 'periodic2_operator_build_count=', &
        electrostatic_diagnostics%periodic_operator_build_count
      write (u, '(a,a)') 'periodic2_cache_fingerprint=', &
        trim(electrostatic_diagnostics%periodic_cache_fingerprint)
      write (u, '(a,a)') 'periodic2_cache_path=', trim(electrostatic_diagnostics%periodic_cache_path)
    end if
    if (cfg%outer_plasma%photoelectron_histogram_enabled) then
      write (u, '(a,i0)') 'photoelectron_histogram_bins=', photoelectron_state%cumulative%nbins
      write (u, '(a,es24.16)') 'photoelectron_histogram_energy_max_J=', &
        photoelectron_state%cumulative%energy_max
      write (u, '(a,i0)') 'photoelectron_last_completed_batch=', photoelectron_state%last_completed_batch
      write (u, '(a,es24.16)') 'photoelectron_cumulative_signed_charge_C=', &
        photoelectron_state%cumulative%total_signed_charge()
      write (u, '(a,es24.16)') 'photoelectron_cumulative_kinetic_energy_J=', &
        photoelectron_state%cumulative%total_kinetic_energy()
      write (u, '(a,i0)') 'photoelectron_cumulative_count=', photoelectron_state%cumulative%total_count()
      if (cfg%sim%batch_duration > 0.0_dp) then
        write (u, '(a,es24.16)') 'photoelectron_previous_signed_current_A=', &
          photoelectron_state%previous_batch%total_signed_charge()/cfg%sim%batch_duration
      end if
      if (cfg%outer_plasma%photoelectron_ambient_charge_scale > 0.0_dp) then
        write (u, '(a,es24.16)') 'photoelectron_previous_charge_ratio=', &
          abs(photoelectron_state%previous_batch%total_signed_charge())/ &
          cfg%outer_plasma%photoelectron_ambient_charge_scale
        write (u, '(a,es24.16)') 'photoelectron_max_charge_ratio=', &
          cfg%outer_plasma%max_photoelectron_charge_ratio
        write (u, '(a)') 'photoelectron_linear_applicability_status=applicable'
      end if
    end if
    if (present(charge_ledger)) then
      write (u, '(a,i0)') 'charge_ledger_nspecies=', charge_ledger%nspecies
      write (u, '(a,i0)') 'charge_ledger_batch_count=', charge_ledger%batch_count
      write (u, '(a,es24.16)') 'charge_ledger_surface_charge_before_C=', charge_ledger%surface_charge_before
      write (u, '(a,es24.16)') 'charge_ledger_surface_charge_after_C=', charge_ledger%surface_charge_after
      write (u, '(a,es24.16)') 'charge_ledger_local_flight_charge_before_C=', &
        charge_ledger%local_flight_charge_before
      write (u, '(a,es24.16)') 'charge_ledger_local_flight_charge_after_C=', &
        charge_ledger%local_flight_charge_after
      write (u, '(a,es24.16)') 'charge_ledger_outer_flight_charge_before_C=', &
        charge_ledger%outer_flight_charge_before
      write (u, '(a,es24.16)') 'charge_ledger_outer_flight_charge_after_C=', &
        charge_ledger%outer_flight_charge_after
      write (u, '(a,es24.16)') 'charge_ledger_unresolved_stock_before_C=', &
        charge_ledger%unresolved_stock_before
      write (u, '(a,es24.16)') 'charge_ledger_unresolved_stock_after_C=', &
        charge_ledger%unresolved_stock_after
      write (u, '(a,es24.16)') 'charge_ledger_residual_C=', charge_ledger%residual()
      write (u, '(a,es24.16)') 'charge_ledger_discarded_unresolved_abs_C=', &
        charge_ledger%discarded_unresolved_abs()
    end if
    if (count_dielectric_surfaces(mesh) > 0_i32) then
      write (u, '(a,i0)') 'surface_model_dielectric_elem_count=', count_dielectric_surfaces(mesh)
      write (u, '(a)') 'surface_model_note=metadata_only_dielectric_present'
    end if
    close (u)
  end subroutine write_summary_file

  !> 解決済みの流入モデルを summary 用の安定した語彙へ変換する。
  function external_inflow_map_name(inflow_map) result(name)
    integer(i32), intent(in) :: inflow_map
    character(len=32) :: name

    select case (inflow_map)
    case (external_inflow_none)
      name = 'source_vdf'
    case (external_inflow_scalar_barrier)
      name = 'infinity_barrier'
    case (external_inflow_legacy_sheath)
      name = 'legacy_sheath'
    case (external_inflow_linear_profile)
      name = 'linear_profile'
    case (external_inflow_kinetic_profile)
      name = 'kinetic_profile'
    case default
      error stop 'write_summary_file: unknown resolved external inflow map.'
    end select
  end function external_inflow_map_name

  !> 解決済みの通常 open 境界モデルを summary 用の安定した語彙へ変換する。
  function external_open_model_name(open_model) result(name)
    integer(i32), intent(in) :: open_model
    character(len=32) :: name

    select case (open_model)
    case (external_open_escape)
      name = 'escape'
    case (external_open_potential_barrier)
      name = 'potential_barrier'
    case default
      error stop 'write_summary_file: unknown resolved ordinary open model.'
    end select
  end function external_open_model_name

  !> 解決済みの interface 輸送モデルを summary 用の安定した語彙へ変換する。
  function external_transport_name(interface_transport) result(name)
    integer(i32), intent(in) :: interface_transport
    character(len=32) :: name

    select case (interface_transport)
    case (external_transport_none)
      name = 'none'
    case (external_transport_linear_1d)
      name = 'linear_1d'
    case (external_transport_kinetic_1d)
      name = 'kinetic_1d'
    case (external_transport_unified_3d)
      name = 'unified_3d'
    case default
      error stop 'write_summary_file: unknown resolved external interface transport.'
    end select
  end function external_transport_name

  !> 粒子外部境界の処理時機を、queue ownership を含めて一意な語彙へ変換する。
  function outer_particle_mode_name(contract) result(name)
    type(external_boundary_contract_type), intent(in) :: contract
    character(len=16) :: name

    if (contract%queue_enabled) then
      name = 'zhao_queue'
    else if (contract%interface_transport /= external_transport_none) then
      name = 'same_batch'
    else
      name = 'local_source'
    end if
  end function outer_particle_mode_name

  !> species 別の signed charge flux と粒子数を `charge_ledger.csv` に保存する。
  subroutine write_charge_ledger_file(out_dir, ledger)
    character(len=*), intent(in) :: out_dir
    type(charge_ledger_type), intent(in) :: ledger
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
      'escaped_to_infinity_C,discarded_unresolved_C,interface_outward_gross_C,interface_returned_gross_C,'// &
      'injected_count,emitted_count,absorbed_count,escaped_count,discarded_unresolved_count'
    do species_idx = 1, ledger%nspecies
      write (u, '(i0,a,i0,7(a,es24.16),5(a,i0))') &
        ledger%batch_count, ',', species_idx, &
        ',', ledger%injected_from_remote(species_idx), &
        ',', ledger%emitted_from_surface(species_idx), &
        ',', ledger%absorbed_on_surface(species_idx), &
        ',', ledger%escaped_to_infinity(species_idx), &
        ',', ledger%discarded_unresolved(species_idx), &
        ',', ledger%interface_outward_gross(species_idx), &
        ',', ledger%interface_returned_gross(species_idx), &
        ',', ledger%injected_count(species_idx), &
        ',', ledger%emitted_count(species_idx), &
        ',', ledger%absorbed_count(species_idx), &
        ',', ledger%escaped_count(species_idx), &
        ',', ledger%discarded_unresolved_count(species_idx)
    end do
    close (u)
  end subroutine write_charge_ledger_file

  !> 要素電荷を `charges.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素電荷を含むメッシュ情報。
  subroutine write_charges_file(out_dir, mesh)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
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
  end subroutine write_charges_file

  !> 事前計算済み電位を `mesh_potential.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素情報（要素数の検証用）。
  !! @param[in] potential_v 各要素重心での電位 [V]。
  subroutine write_mesh_potential_file(out_dir, mesh, potential_v)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: potential_v(:)
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
  end subroutine write_mesh_potential_file

  !> 三角形メッシュを `mesh_triangles.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 頂点座標と要素電荷を含むメッシュ情報。
  subroutine write_mesh_file(out_dir, mesh)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
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
  end subroutine write_mesh_file

  !> メッシュ識別情報を `mesh_sources.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素ごとの `mesh_id` を含むメッシュ情報。
  !! @param[in] cfg 元の入力設定。
  subroutine write_mesh_sources_file(out_dir, mesh, cfg)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(app_config), intent(in) :: cfg
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
  end subroutine write_mesh_sources_file

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

  !> 現行物理モデルで未分岐の dielectric 要素数を数える。
  integer(i32) function count_dielectric_surfaces(mesh) result(n)
    type(mesh_type), intent(in) :: mesh

    n = 0_i32
    if (.not. allocated(mesh%elem_surface_model)) return
    n = int(count(mesh%elem_surface_model == surface_model_dielectric), kind=i32)
  end function count_dielectric_surfaces

end module bem_output_writer
