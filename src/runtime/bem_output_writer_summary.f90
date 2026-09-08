!> 実行サマリと再構築用メタデータの出力。
submodule(bem_output_writer) bem_output_writer_summary
  use bem_types, only: surface_model_dielectric
  use bem_charge_ledger, only: finite_charge_sum
  use bem_checkpoint_contract, only: checkpoint_schema_version_current
  use bem_field_solver, only: &
    field_solver_fmm_expansion_order, &
    resolve_field_solver_mode, &
    resolve_field_solver_tree_params
  use bem_external_boundary_contract, only: &
    external_boundary_contract_type, &
    external_boundary_ok, &
    external_inflow_none, &
    external_inflow_scalar_barrier, &
    external_open_escape, &
    external_open_potential_barrier, &
    resolve_external_boundary_contract
  use bem_model_fingerprint, only: model_fingerprint, mesh_fingerprint, species_fingerprint
  use bem_matching_plane_response, only: get_matching_plane_response_content_fingerprint, matching_plane_response_ok
  use bem_physics_config_types, only: field_physics_config, panel_kernel_config, derive_field_panel_config
  use bem_surface_current_model, only: surface_current_model_result_type, evaluate_surface_current_model
  use bem_version, only: beach_build_id, beach_source_commit, beach_version, beach_version_mode
  use bem_string_utils, only: lower_ascii
  implicit none
contains

  !> 実行統計を `summary.txt` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh メッシュ情報（要素数を書き出す）。
  !! @param[in] stats 実行統計。
  module procedure write_summary_file
  type(external_boundary_contract_type) :: resolved_boundary
  type(surface_current_model_result_type) :: current_model
  type(field_physics_config) :: field_config
  type(panel_kernel_config) :: panel_config
  character(len=1024) :: summary_path
  character(len=512) :: matching_response_message
  character(len=256) :: boundary_message
  character(len=16) :: resolved_field_solver, matching_response_fingerprint
  integer :: u, ios
  integer(i32) :: world_size, boundary_status, resolved_tree_leaf_max, matching_response_status
  real(dp) :: resolved_tree_theta

  call resolve_external_boundary_contract( &
    cfg%sim%reservoir_potential_model, cfg%sim%open_boundary_model, &
    resolved_boundary, boundary_status, boundary_message &
    )
  if (boundary_status /= external_boundary_ok) then
    error stop 'write_summary_file: invalid local boundary contract: '//trim(boundary_message)
  end if
  call evaluate_surface_current_model(cfg, current_model)
  call derive_field_panel_config(cfg%sim, field_config, panel_config)
  call resolve_field_solver_tree_params( &
    mesh%nelem, cfg%sim, resolved_tree_theta, resolved_tree_leaf_max &
    )
  resolved_field_solver = resolve_field_solver_mode(mesh%nelem, cfg%sim)
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
  write (u, '(a,i0)') 'checkpoint_stride=', cfg%checkpoint_stride
  write (u, '(a,i0)') 'escaped_boundary=', stats%escaped_boundary
  write (u, '(a,i0)') 'survived_max_step=', stats%survived_max_step
  write (u, '(a,i0)') 'multiple_box_events_retry_attempted=', stats%multiple_box_events_retry_attempted
  write (u, '(a,i0)') 'multiple_box_events_retry_resolved=', stats%multiple_box_events_retry_resolved
  write (u, '(a,i0)') 'multiple_box_events_soft_discarded=', stats%multiple_box_events_soft_discarded
  write (u, '(a,es24.16)') 'multiple_box_events_soft_discard_fraction=', soft_discard_fraction(stats)
  write (u, '(a,es24.16)') 'multiple_box_events_soft_discarded_abs_charge_C=', &
    stats%multiple_box_events_soft_discarded_abs_charge
  write (u, '(a,es24.16)') 'last_rel_change=', stats%last_rel_change
  write (u, '(a,es24.16)') 'simulated_time_s=', stats%simulated_time
  write (u, '(a,i0)') 'adaptive_nonzero_mode_rejected_trials=', &
    stats%adaptive_nonzero_mode_rejected_trials
  write (u, '(a,es24.16)') 'adaptive_nonzero_mode_last_batch_duration_s=', &
    stats%adaptive_nonzero_mode_last_batch_duration
  write (u, '(a,es24.16)') 'adaptive_nonzero_mode_last_potential_step_V=', &
    stats%adaptive_nonzero_mode_last_potential_step
  write (u, '(a,i0)') 'adaptive_nonzero_mode_omp_threads=', &
    stats%adaptive_nonzero_mode_omp_threads
  write (u, '(a,l1)') 'matching_plane_state_valid=', stats%matching_plane_state_valid
  write (u, '(a,es24.16)') 'matching_plane_displacement_C_m2=', stats%matching_plane_displacement_c_m2
  write (u, '(a,es24.16)') 'matching_plane_phi_V=', stats%matching_plane_phi_v
  write (u, '(a,es24.16)') 'matching_plane_electron_inward_flux_m2_s=', stats%matching_plane_response(2)
  write (u, '(a,es24.16)') 'matching_plane_ion_inward_flux_m2_s=', stats%matching_plane_response(3)
  write (u, '(a,es24.16)') 'matching_plane_electron_access_potential_V=', stats%matching_plane_response(4)
  write (u, '(a,es24.16)') 'matching_plane_ion_access_potential_V=', stats%matching_plane_response(5)
  write (u, '(a,es24.16)') 'matching_plane_photoelectron_barrier_potential_V=', &
    stats%matching_plane_response(6)
  write (u, '(a,es24.16)') 'matching_plane_photoelectron_outward_flux_m2_s=', &
    stats%matching_plane_feedback(1)
  write (u, '(a,es24.16)') 'matching_plane_photoelectron_mean_normal_energy_eV=', &
    stats%matching_plane_feedback(2)
  write (u, '(a,es24.16)') 'matching_plane_electron_outward_flux_m2_s=', &
    stats%matching_plane_feedback(3)
  write (u, '(a,es24.16)') 'matching_plane_ion_outward_flux_m2_s=', &
    stats%matching_plane_feedback(4)
  write (u, '(a,es24.16)') 'matching_plane_photoelectron_return_flux_m2_s=', &
    stats%matching_plane_photoelectron_return_flux_m2_s
  write (u, '(a,es24.16)') 'matching_plane_photoelectron_escape_flux_m2_s=', &
    stats%matching_plane_photoelectron_escape_flux_m2_s
  write (u, '(a,i0)') 'matching_plane_iterations=', stats%matching_plane_iterations
  write (u, '(a,es24.16)') 'matching_plane_residual=', stats%matching_plane_residual
  write (u, '(a)') 'particle_time_centering=same_time_midpoint_boris'
  write (u, '(a,a)') 'field_backend=', trim(field_config%backend)
  write (u, '(a,a)') 'field_normalization=', trim(field_config%normalization)
  write (u, '(a)') 'field_source_model=triangle_p0'
  write (u, '(a,a)') 'field_kernel_id=', trim(panel_config%kernel_id)
  write (u, '(a)') 'field_reconstruction_schema_version=2'
  write (u, '(a,a)') 'field_reconstruction_resolved_field_solver=', trim(resolved_field_solver)
  write (u, '(a,i0)') 'field_reconstruction_fmm_expansion_order=', field_solver_fmm_expansion_order
  write (u, '(a,a)') 'field_reconstruction_field_bc_mode=', trim(cfg%sim%field_bc_mode)
  write (u, '(a,es24.16)') 'field_reconstruction_tree_theta=', resolved_tree_theta
  write (u, '(a,i0)') 'field_reconstruction_tree_leaf_max=', resolved_tree_leaf_max
  write (u, '(a,3(1x,es24.16))') 'field_reconstruction_e0_V_m=', cfg%sim%e0
  write (u, '(a,l1)') 'field_reconstruction_use_box=', cfg%sim%use_box
  write (u, '(a,3(1x,es24.16))') 'field_reconstruction_box_min_m=', cfg%sim%box_min
  write (u, '(a,3(1x,es24.16))') 'field_reconstruction_box_max_m=', cfg%sim%box_max
  write (u, '(a,3(1x,i0))') 'field_reconstruction_boundary_low=', cfg%sim%bc_low
  write (u, '(a,3(1x,i0))') 'field_reconstruction_boundary_high=', cfg%sim%bc_high
  write (u, '(a,i0)') 'field_reconstruction_periodic_image_layers=', cfg%sim%field_periodic_image_layers
  write (u, '(a,a)') 'field_reconstruction_periodic_far_correction=', &
    trim(cfg%sim%field_periodic_far_correction)
  write (u, '(a,a)') 'field_reconstruction_periodic_nonzero_mode_backend=', &
    trim(cfg%periodic2%nonzero_mode_backend)
  write (u, '(a,a)') 'field_reconstruction_periodic_zero_mode_policy=', &
    trim(cfg%periodic2%zero_mode_policy)
  write (u, '(a,a)') 'field_reconstruction_periodic_lower_boundary_model=', &
    trim(cfg%periodic2%lower_boundary_model)
  write (u, '(a,i0)') 'field_reconstruction_periodic_reference_mode_layers=', &
    cfg%periodic2%reference_mode_layers
  write (u, '(a,i0)') 'field_reconstruction_periodic_panel_quadrature_order=', &
    cfg%periodic2%panel_quadrature_order
  write (u, '(a,es24.16)') 'field_reconstruction_periodic_ewald_alpha=', cfg%sim%field_periodic_ewald_alpha
  write (u, '(a,i0)') 'field_reconstruction_periodic_ewald_layers=', cfg%sim%field_periodic_ewald_layers
  write (u, '(a,a)') 'field_reconstruction_periodic_cache_dir=', trim(cfg%sim%field_periodic_cache_dir)
  write (u, '(a,es24.16)') 'field_reconstruction_periodic_generation_tolerance=', &
    cfg%sim%field_periodic_generation_tolerance
  write (u, '(a,a)') 'periodic2_nonzero_mode_backend=', trim(cfg%periodic2%nonzero_mode_backend)
  write (u, '(a,a)') 'periodic2_zero_mode_policy=', trim(cfg%periodic2%zero_mode_policy)
  write (u, '(a,a)') 'periodic2_lower_boundary_model=', trim(cfg%periodic2%lower_boundary_model)
  write (u, '(a,es24.16)') 'periodic2_max_nonzero_mode_potential_step_V=', &
    cfg%periodic2%max_nonzero_mode_potential_step
  write (u, '(a,a)') 'periodic2_cache_dir=', trim(cfg%sim%field_periodic_cache_dir)
  write (u, '(a,es24.16)') 'periodic2_generation_tolerance=', &
    cfg%sim%field_periodic_generation_tolerance
  write (u, '(a,a)') 'reservoir_inflow_map=', trim(external_inflow_map_name(resolved_boundary%inflow_map))
  write (u, '(a,a)') 'particle_ordinary_open_model=', &
    trim(external_open_model_name(resolved_boundary%ordinary_open_model))
  write (u, '(a,a)') 'surface_current_model=', trim(current_model%model)
  if (current_model%active) then
    write (u, '(a,a)') 'surface_current_model_kinetic_contract=', trim(current_model%kinetic_contract)
    if (trim(current_model%model) == 'matching_plane_quasistatic') then
      write (u, '(a,a)') 'surface_current_model_response_backend=', &
        trim(lower_ascii(cfg%surface_current%response_backend))
      write (u, '(a,l1)') 'surface_current_model_implicit_zero_mode=', &
        cfg%surface_current%implicit_zero_mode
      select case (trim(lower_ascii(cfg%surface_current%response_backend)))
      case ('table')
        call get_matching_plane_response_content_fingerprint( &
          trim(cfg%surface_current%response_table_path), matching_response_fingerprint, &
          matching_response_status, matching_response_message &
          )
        if (matching_response_status /= matching_plane_response_ok) then
          error stop 'write_summary_file: matching-plane response fingerprint failed: '// &
            trim(matching_response_message)
        end if
        write (u, '(a,a)') 'surface_current_model_response_table_path=', &
          trim(cfg%surface_current%response_table_path)
        write (u, '(a,a)') 'surface_current_model_response_content_fingerprint=', &
          matching_response_fingerprint
      case ('zhao_online')
        write (u, '(a)') 'surface_current_model_response_contract=matching_plane_zhao_online_v1'
        write (u, '(a,a)') 'surface_current_model_zhao_branch=', &
          trim(lower_ascii(cfg%surface_current%zhao_branch))
        write (u, '(a,a)') 'surface_current_model_zhao_root_selection=', &
          trim(lower_ascii(cfg%surface_current%zhao_root_selection))
        write (u, '(a)') &
          'surface_current_model_outer_solver=charge_driven_finite_h_sagdeev'
        write (u, '(a)') &
          'surface_current_model_photoelectron_closure=moment_matched_half_maxwellian'
        write (u, '(a)') &
          'surface_current_model_ambient_outward_feedback=transparent'
        if (trim(lower_ascii(cfg%surface_current%zhao_root_selection)) == 'continuation') then
          write (u, '(a)') &
            'surface_current_model_outer_solver_state=accepted_endpoint_continuation_v2'
        else
          write (u, '(a)') 'surface_current_model_outer_solver_state=stateless'
        end if
      case default
        error stop 'write_summary_file: unknown matching-plane response backend.'
      end select
      write (u, '(a,es24.16)') 'surface_current_model_matching_plane_z_m=', cfg%sim%box_max(3)
      write (u, '(a,a)') 'surface_current_model_electron_species=', &
        trim(cfg%surface_current%electron_species)
      write (u, '(a,a)') 'surface_current_model_ion_species=', trim(cfg%surface_current%ion_species)
      write (u, '(a,a)') 'surface_current_model_photoelectron_species=', &
        trim(cfg%surface_current%photoelectron_species)
      write (u, '(a,es24.16)') 'surface_current_model_coupling_rtol=', cfg%surface_current%coupling_rtol
      write (u, '(a,4(1x,es24.16))') 'surface_current_model_coupling_atol=', &
        cfg%surface_current%coupling_atol
      write (u, '(a,i0)') 'surface_current_model_coupling_max_iterations=', &
        cfg%surface_current%coupling_max_iterations
      write (u, '(a,es24.16)') 'surface_current_model_coupling_relaxation=', &
        cfg%surface_current%coupling_relaxation
      write (u, '(a)') 'surface_current_model_dynamic_state_source=accepted_batch_fixed_point'
    else
      write (u, '(a,a)') 'surface_current_model_zhao_branch=', current_model%zhao_branch
      write (u, '(a,l1)') 'surface_current_model_photoelectron_active=', current_model%photoelectron_active
      write (u, '(a,es24.16)') 'surface_current_model_reference_area_m2=', current_model%reference_area_m2
      write (u, '(a,es24.16)') 'surface_current_model_phi0_V=', current_model%phi0_v
      write (u, '(a,es24.16)') 'surface_current_model_phi_m_V=', current_model%phi_m_v
      write (u, '(a,es24.16)') 'surface_current_model_ambient_electron_density_m3=', &
        current_model%ambient_electron_density_m3
      write (u, '(a,es24.16)') 'surface_current_model_electron_current_density_A_m2=', &
        current_model%electron_current_density_a_m2
      write (u, '(a,es24.16)') 'surface_current_model_ion_current_density_A_m2=', &
        current_model%ion_current_density_a_m2
      write (u, '(a,es24.16)') 'surface_current_model_pe_emission_current_density_A_m2=', &
        current_model%photoelectron_emission_current_density_a_m2
      write (u, '(a,es24.16)') 'surface_current_model_pe_escape_current_density_A_m2=', &
        current_model%photoelectron_escape_current_density_a_m2
      write (u, '(a,es24.16)') 'surface_current_model_pe_return_current_density_A_m2=', &
        current_model%photoelectron_return_current_density_a_m2
      write (u, '(a,es24.16)') 'surface_current_model_net_current_density_A_m2=', &
        current_model%net_current_density_a_m2
      write (u, '(a)') 'surface_current_model_current_budget_contract=surface_targets_plus_external_escape'
      write (u, '(a,es24.16)') 'surface_current_model_pe_budget_residual_current_density_A_m2=', &
        current_model%photoelectron_budget_residual_current_density_a_m2
      write (u, '(a,es24.16)') 'surface_current_model_surface_budget_residual_current_density_A_m2=', &
        current_model%surface_budget_residual_current_density_a_m2
      if (current_model%photoelectron_active) then
        write (u, '(a,es24.16)') 'surface_current_model_pe_escape_particle_current_A=', &
          current_model%escaped_particle_current_a(current_model%photoelectron_species_idx)
      else
        write (u, '(a,es24.16)') 'surface_current_model_pe_escape_particle_current_A=', 0.0_dp
      end if
      write (u, '(a,es24.16)') 'surface_current_model_electron_inflow_reservoir_potential_V=', &
        current_model%inflow_reservoir_potential_v(current_model%electron_species_idx)
      write (u, '(a,es24.16)') 'surface_current_model_electron_inflow_access_potential_V=', &
        current_model%inflow_access_potential_v(current_model%electron_species_idx)
      write (u, '(a,i0)') 'surface_current_model_electron_inflow_face=', &
        current_model%inflow_kinetic_face(current_model%electron_species_idx)
      if (current_model%photoelectron_active) then
        write (u, '(a,es24.16)') 'surface_current_model_pe_outflow_barrier_potential_V=', &
          current_model%outflow_barrier_potential_v(current_model%photoelectron_species_idx)
        write (u, '(a,i0)') 'surface_current_model_pe_outflow_barrier_face=', &
          current_model%outflow_barrier_face(current_model%photoelectron_species_idx)
      else
        write (u, '(a,es24.16)') 'surface_current_model_pe_outflow_barrier_potential_V=', 0.0_dp
        write (u, '(a,i0)') 'surface_current_model_pe_outflow_barrier_face=', 0_i32
      end if
    end if
  end if
  if (present(electrostatic_diagnostics)) then
    write (u, '(a,l1)') 'top_reference_available=', electrostatic_diagnostics%top_reference_available
    if (electrostatic_diagnostics%top_reference_available) then
      write (u, '(a)') 'top_reference_definition=box_z_high_plane_mean'
      write (u, '(a,i0)') 'top_reference_last_batch=', electrostatic_diagnostics%top_reference_last_batch
      write (u, '(a,es24.16)') 'top_reference_simulated_time_s=', electrostatic_diagnostics%top_reference_simulated_time
      write (u, '(a,es24.16)') 'top_reference_z_high_m=', electrostatic_diagnostics%top_reference_z_high
      write (u, '(a,i0)') 'top_reference_sample_n=', electrostatic_diagnostics%top_reference_sample_n
      write (u, '(a,es24.16)') 'top_reference_potential_mean_V=', electrostatic_diagnostics%top_reference_potential_mean
      write (u, '(a,es24.16)') 'top_reference_potential_std_V=', electrostatic_diagnostics%top_reference_potential_std
      write (u, '(a,es24.16)') 'top_reference_potential_min_V=', electrostatic_diagnostics%top_reference_potential_min
      write (u, '(a,es24.16)') 'top_reference_potential_max_V=', electrostatic_diagnostics%top_reference_potential_max
    end if
    write (u, '(a,l1)') 'electrostatic_split_periodic_active=', electrostatic_diagnostics%split_periodic_active
    write (u, '(a,a)') 'electrostatic_status=', trim(electrostatic_diagnostics%status)
    write (u, '(a,es24.16)') 'gauss_residual_C=', electrostatic_diagnostics%gauss_residual
    write (u, '(a,l1)') 'periodic2_cache_hit=', electrostatic_diagnostics%periodic_cache_hit
    write (u, '(a,i0)') 'periodic2_operator_build_count=', electrostatic_diagnostics%periodic_operator_build_count
    write (u, '(a,a)') 'periodic2_cache_fingerprint=', trim(electrostatic_diagnostics%periodic_cache_fingerprint)
    write (u, '(a,a)') 'periodic2_cache_path=', trim(electrostatic_diagnostics%periodic_cache_path)
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
    write (u, '(a,es24.16)') 'charge_ledger_unresolved_stock_before_C=', &
      charge_ledger%unresolved_stock_before
    write (u, '(a,es24.16)') 'charge_ledger_unresolved_stock_after_C=', &
      charge_ledger%unresolved_stock_after
    write (u, '(a,es24.16)') 'charge_ledger_residual_C=', charge_ledger%residual()
    write (u, '(a,es24.16)') 'charge_ledger_discarded_unresolved_abs_C=', &
      charge_ledger%discarded_unresolved_abs()
    write (u, '(a,es24.16)') 'charge_ledger_neutral_return_correction_C=', &
      finite_charge_sum(charge_ledger%neutral_return_correction, 'summary neutral-return correction')
    write (u, '(a,es24.16)') 'charge_ledger_fixed_current_correction_C=', &
      finite_charge_sum(charge_ledger%fixed_current_correction, 'summary fixed-current correction')
    write (u, '(a,es24.16)') 'charge_ledger_fixed_absorbed_applied_charge_C=', &
      finite_charge_sum(charge_ledger%fixed_absorbed_target_charge, 'summary fixed absorbed applied charge')
    write (u, '(a,es24.16)') 'charge_ledger_fixed_emission_applied_charge_C=', &
      finite_charge_sum(charge_ledger%fixed_emission_target_charge, 'summary fixed emission applied charge')
    write (u, '(a,es24.16)') 'charge_ledger_raw_escape_charge_C=', &
      finite_charge_sum(charge_ledger%escaped_to_infinity, 'summary raw escape charge')
    write (u, '(a,es24.16)') 'charge_ledger_fixed_escape_target_charge_C=', &
      finite_charge_sum(charge_ledger%fixed_escape_target_charge, 'summary fixed escape target charge')
    write (u, '(a,es24.16)') 'charge_ledger_fixed_escape_applied_charge_C=', &
      finite_charge_sum(charge_ledger%fixed_escape_target_charge, 'summary fixed escape applied charge')
    write (u, '(a,es24.16)') 'charge_ledger_fixed_escape_correction_C=', &
      finite_charge_sum(charge_ledger%fixed_escape_correction, 'summary fixed escape correction')
    write (u, '(a,es24.16)') 'charge_ledger_fixed_applied_surface_net_charge_C=', &
      finite_charge_sum( &
      [charge_ledger%fixed_absorbed_target_charge, charge_ledger%fixed_emission_target_charge], &
      'summary fixed applied surface net charge' &
      )
    if (current_model%active .and. current_model%photoelectron_active) then
      write (u, '(a,es24.16)') 'charge_ledger_fixed_pe_continuity_residual_C=', &
        finite_charge_sum( &
        [ &
        charge_ledger%fixed_absorbed_target_charge(current_model%photoelectron_species_idx), &
        charge_ledger%fixed_emission_target_charge(current_model%photoelectron_species_idx), &
        charge_ledger%fixed_escape_target_charge(current_model%photoelectron_species_idx) &
        ], &
        'summary fixed photoelectron continuity residual' &
        )
    end if
  end if
  if (count_dielectric_surfaces(mesh) > 0_i32) then
    write (u, '(a,i0)') 'surface_model_dielectric_elem_count=', count_dielectric_surfaces(mesh)
    write (u, '(a)') 'surface_model_note=metadata_only_dielectric_present'
  end if
  close (u)
  end procedure write_summary_file

  !> 解決済みの流入モデルを summary 用の安定した語彙へ変換する。
  function external_inflow_map_name(inflow_map) result(name)
    integer(i32), intent(in) :: inflow_map
    character(len=32) :: name

    select case (inflow_map)
    case (external_inflow_none)
      name = 'source_vdf'
    case (external_inflow_scalar_barrier)
      name = 'infinity_barrier'
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

  !> 処理済みmacro particleに対する累積soft-discard率を返す。
  module procedure soft_discard_fraction

  fraction = 0.0_dp
  if (stats%processed_particles <= 0) return
  fraction = real(stats%multiple_box_events_soft_discarded, dp)/real(stats%processed_particles, dp)
  end procedure soft_discard_fraction

  !> 現行物理モデルで未分岐の dielectric 要素数を数える。
  module procedure count_dielectric_surfaces

  n = 0_i32
  if (.not. allocated(mesh%elem_surface_model)) return
  n = int(count(mesh%elem_surface_model == surface_model_dielectric), kind=i32)
  end procedure count_dielectric_surfaces

end submodule bem_output_writer_summary
