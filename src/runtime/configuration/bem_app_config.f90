!> 設定型・TOMLパーサ・実行時変換ロジックを束ねる後方互換ファサード。
module bem_app_config
  use bem_app_config_types, only: &
    max_templates, max_particle_species, particle_bc_inherit, particle_inflow_none, particle_inflow_reservoir, &
    particle_species_spec, template_spec, app_config, &
    particles_per_batch_from_config, default_app_config, species_from_defaults
  use bem_physics_config_types, only: &
    field_physics_config, periodic2_physics_config, panel_kernel_config
  use bem_app_config_parser, only: load_app_config, load_toml_config
  use bem_string_utils, only: lower_ascii
  use bem_app_config_runtime, only: &
    build_mesh_from_config, seed_particles_from_config, &
    init_particle_batch_from_config, build_parts_mesh, build_one_template, &
    batch_injection_layout_type, build_batch_injection_layout, generate_batch_slot, &
    slot_kind_none, slot_kind_volume, slot_kind_photo_ray, slot_kind_reservoir_face, &
    append_triangles
  implicit none
end module bem_app_config
