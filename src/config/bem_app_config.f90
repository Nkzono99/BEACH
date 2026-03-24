!> 設定型・TOMLパーサ・実行時変換ロジックを束ねる後方互換ファサード。
module bem_app_config
  use bem_app_config_types, only: &
    max_templates, max_particle_species, particle_species_spec, template_spec, app_config, &
    particles_per_batch_from_config, total_particles_from_config, default_app_config, species_from_defaults
  use bem_app_config_parser, only: &
    load_app_config, load_toml_config, apply_sim_kv, apply_particles_kv, apply_particles_species_kv, &
    apply_mesh_kv, apply_template_kv, apply_output_kv, split_key_value, parse_real, parse_int, &
    parse_logical, parse_string, parse_real3, parse_boundary_mode, strip_comment, ends_with
  use bem_string_utils, only: lower_ascii
  use bem_app_config_runtime, only: &
    build_mesh_from_config, init_particles_from_config, seed_particles_from_config, &
    init_particle_batch_from_config, sample_species_state, build_template_mesh, build_one_template, &
    append_triangles
  implicit none
end module bem_app_config
