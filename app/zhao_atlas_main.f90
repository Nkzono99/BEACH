!> Zhao A/B/C solvability atlasを生成するCLI。
program zhao_atlas_main
  use bem_kinds, only: i32
  use bem_app_config, only: app_config, default_app_config, load_app_config
  use bem_matching_plane_zhao_atlas, only: generate_matching_plane_zhao_atlas, &
                                           matching_plane_atlas_ok
  implicit none

  type(app_config) :: cfg
  character(len=1024) :: config_path, query_path, output_path
  character(len=512) :: message
  integer(i32) :: status

  if (command_argument_count() == 1) then
    call get_command_argument(1, config_path)
    if (trim(config_path) == '--help' .or. trim(config_path) == '-h') then
      call print_usage()
      stop
    end if
  end if
  if (command_argument_count() /= 3) then
    call print_usage()
    error stop 'beach-zhao-atlas requires exactly three positional arguments.'
  end if

  call get_command_argument(1, config_path)
  call get_command_argument(2, query_path)
  call get_command_argument(3, output_path)
  call default_app_config(cfg)
  call load_app_config(trim(config_path), cfg)
  call generate_matching_plane_zhao_atlas( &
    cfg, trim(query_path), trim(output_path), status, message &
    )
  if (status /= matching_plane_atlas_ok) error stop trim(message)
  print '(a,a)', 'matching-plane Zhao atlas written to ', trim(output_path)

contains

  subroutine print_usage()
    print '(a)', 'usage: beach-zhao-atlas <beach.toml> <query-grid.csv> <atlas.csv>'
  end subroutine print_usage

end program zhao_atlas_main
