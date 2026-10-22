!> app_config のメッシュ・粒子・電位評価を束ねる後方互換ファサード。
module bem_app_config_runtime
  use bem_app_config_mesh_runtime, only: &
    build_mesh_from_config, apply_panel_surface_config, apply_obj_transform, &
    build_parts_mesh, build_one_template, &
    append_triangles, append_mesh_ids, append_real_values, surface_model_id_from_name
  use bem_app_config_particle_runtime, only: &
    particle_source_plan_type, seed_particles_from_config, &
    batch_injection_layout_type, build_batch_injection_layout, generate_batch_slot, &
    slot_kind_none, slot_kind_volume, slot_kind_photo_ray, slot_kind_reservoir_face, &
    build_particle_source_plan, init_particle_batch_from_config, sample_species_state, &
    has_boundary_inflow, boundary_inflow_face_enabled, &
    make_boundary_inflow_spec, boundary_face_name, &
    sample_photo_species_state, compute_macro_particles_for_species, reservoir_face_velocity_correction, &
    species_number_density_m3, resolve_parallel_rank_size, species_temperature_k
  use bem_app_config_potential_runtime, only: &
    compute_face_average_potential, compute_wall_box_potential_statistics, &
    warn_face_average_potential_variation, resolve_face_sampling_geometry
  implicit none
end module bem_app_config_runtime
