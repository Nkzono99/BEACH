!> 領域別の TOML 読み取り。意味検証と派生値の確定は preflight が担当する。
submodule(bem_app_config_parser) bem_app_config_parser_read_surface
  implicit none
contains

  module procedure apply_surface_current_model_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('model')
      call get_toml_string(table, keys(ikey), cfg%surface_current%model, 'surface_current_model.model')
      cfg%surface_current%model = lower_ascii(trim(cfg%surface_current%model))
    case ('response_backend')
      call get_toml_string( &
        table, keys(ikey), cfg%surface_current%response_backend, 'surface_current_model.response_backend' &
        )
      cfg%surface_current%response_backend = lower_ascii(trim(cfg%surface_current%response_backend))
      cfg%surface_current%has_response_backend = .true.
    case ('zhao_branch')
      call get_toml_string(table, keys(ikey), cfg%surface_current%zhao_branch, 'surface_current_model.zhao_branch')
      cfg%surface_current%zhao_branch = lower_ascii(trim(cfg%surface_current%zhao_branch))
      cfg%surface_current%has_zhao_branch = .true.
    case ('zhao_root_selection')
      call get_toml_string( &
        table, keys(ikey), cfg%surface_current%zhao_root_selection, &
        'surface_current_model.zhao_root_selection' &
        )
      cfg%surface_current%zhao_root_selection = lower_ascii(trim(cfg%surface_current%zhao_root_selection))
      cfg%surface_current%has_zhao_root_selection = .true.
    case ('electron_species')
      call get_toml_string( &
        table, keys(ikey), cfg%surface_current%electron_species, 'surface_current_model.electron_species' &
        )
      cfg%surface_current%has_electron_species = .true.
    case ('ion_species')
      call get_toml_string(table, keys(ikey), cfg%surface_current%ion_species, 'surface_current_model.ion_species')
      cfg%surface_current%has_ion_species = .true.
    case ('photoelectron_species')
      call get_toml_string( &
        table, keys(ikey), cfg%surface_current%photoelectron_species, &
        'surface_current_model.photoelectron_species' &
        )
      cfg%surface_current%has_photoelectron_species = .true.
    case ('solar_elevation_deg')
      call get_toml_real( &
        table, keys(ikey), cfg%surface_current%solar_elevation_deg, 'surface_current_model.solar_elevation_deg' &
        )
      cfg%surface_current%has_solar_elevation_deg = .true.
    case ('photoelectron_ref_density_m3')
      call get_toml_real( &
        table, keys(ikey), cfg%surface_current%photoelectron_ref_density_m3, &
        'surface_current_model.photoelectron_ref_density_m3' &
        )
      cfg%surface_current%has_photoelectron_ref_density_m3 = .true.
    case ('photoelectron_source_scale')
      call get_toml_real( &
        table, keys(ikey), cfg%surface_current%photoelectron_source_scale, &
        'surface_current_model.photoelectron_source_scale' &
        )
      cfg%surface_current%has_photoelectron_source_scale = .true.
    case ('reference_area_m2')
      call get_toml_real( &
        table, keys(ikey), cfg%surface_current%reference_area_m2, 'surface_current_model.reference_area_m2' &
        )
      cfg%surface_current%has_reference_area_m2 = .true.
    case ('response_table_path')
      call get_toml_string( &
        table, keys(ikey), cfg%surface_current%response_table_path, &
        'surface_current_model.response_table_path' &
        )
      cfg%surface_current%has_response_table_path = .true.
    case ('implicit_zero_mode')
      call get_toml_logical( &
        table, keys(ikey), cfg%surface_current%implicit_zero_mode, &
        'surface_current_model.implicit_zero_mode' &
        )
      cfg%surface_current%has_implicit_zero_mode = .true.
    case ('coupling_rtol')
      call get_toml_real( &
        table, keys(ikey), cfg%surface_current%coupling_rtol, 'surface_current_model.coupling_rtol' &
        )
      cfg%surface_current%has_coupling_rtol = .true.
    case ('coupling_atol')
      call get_toml_real4( &
        table, keys(ikey), cfg%surface_current%coupling_atol, 'surface_current_model.coupling_atol' &
        )
      cfg%surface_current%has_coupling_atol = .true.
    case ('coupling_max_iterations')
      call get_toml_int( &
        table, keys(ikey), cfg%surface_current%coupling_max_iterations, &
        'surface_current_model.coupling_max_iterations' &
        )
      cfg%surface_current%has_coupling_max_iterations = .true.
    case ('coupling_relaxation')
      call get_toml_real( &
        table, keys(ikey), cfg%surface_current%coupling_relaxation, &
        'surface_current_model.coupling_relaxation' &
        )
      cfg%surface_current%has_coupling_relaxation = .true.
    case default
      error stop 'Unknown key in [surface_current_model]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_surface_current_model_toml_table

end submodule bem_app_config_parser_read_surface
