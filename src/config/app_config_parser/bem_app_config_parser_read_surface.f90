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
      call get_toml_string(table, keys(ikey), cfg%surface_current%model, 'sheath.closure')
      cfg%surface_current%model = lower_ascii(trim(cfg%surface_current%model))
    case ('zhao_branch')
      call get_toml_string(table, keys(ikey), cfg%surface_current%zhao_branch, 'sheath.zhao.branch')
      cfg%surface_current%zhao_branch = lower_ascii(trim(cfg%surface_current%zhao_branch))
    case ('electron_species')
      call get_toml_string( &
        table, keys(ikey), cfg%surface_current%electron_species, 'sheath.species.electron' &
        )
    case ('ion_species')
      call get_toml_string(table, keys(ikey), cfg%surface_current%ion_species, 'sheath.species.ion')
    case ('secondary_electron_species')
      call get_toml_string(table, keys(ikey), cfg%surface_current%secondary_electron_species, &
                           'sheath.species.secondary_electron')
    case ('photoelectron_species')
      call get_toml_string( &
        table, keys(ikey), cfg%surface_current%photoelectron_species, &
        'sheath.species.photoelectron' &
        )
    case ('solar_elevation_deg')
      call get_toml_real( &
        table, keys(ikey), cfg%surface_current%solar_elevation_deg, 'sheath.photoelectrons.solar_elevation_deg' &
        )
    case ('photoelectron_ref_density_m3')
      call get_toml_real( &
        table, keys(ikey), cfg%surface_current%photoelectron_ref_density_m3, &
        'sheath.photoelectrons.ref_density_m3' &
        )
    case ('photoelectron_source_scale')
      call get_toml_real( &
        table, keys(ikey), cfg%surface_current%photoelectron_source_scale, &
        'sheath.photoelectrons.source_scale' &
        )
    case ('reference_area_m2')
      call get_toml_real( &
        table, keys(ikey), cfg%surface_current%reference_area_m2, 'sheath.reference_area_m2' &
        )
      cfg%surface_current%has_reference_area_m2 = .true.
    case ('outflow_refresh_batches')
      call get_toml_int( &
        table, keys(ikey), cfg%surface_current%outflow_refresh_batches, &
        'sheath.coupling.outflow_refresh_batches' &
        )
    case ('include_reflected_electrons')
      call get_toml_logical( &
        table, keys(ikey), cfg%surface_current%include_reflected_electrons, &
        'sheath.coupling.include_reflected_electrons' &
        )
    case ('zhao_upstream_band_tolerance')
      call get_toml_real( &
        table, keys(ikey), cfg%surface_current%zhao_upstream_band_tolerance, &
        'sheath.zhao.upstream_band_tolerance' &
        )
    case ('zhao_upstream_field_floor_v_m')
      call get_toml_real( &
        table, keys(ikey), cfg%surface_current%zhao_upstream_field_floor_v_m, &
        'sheath.zhao.upstream_field_floor_v_m' &
        )
    case ('zhao_ion_pressure_factor')
      call get_toml_real( &
        table, keys(ikey), cfg%surface_current%zhao_ion_pressure_factor, &
        'sheath.zhao.ion_pressure_factor' &
        )
    case default
      error stop 'Unknown key in [sheath]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_surface_current_model_toml_table

end submodule bem_app_config_parser_read_surface
