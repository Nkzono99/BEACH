!> 7グループの入力を既存の内部TOMLへ変換する。物理検証と既定値は既存readerが所有する。
!! 公開済みのflat入力は1.xで維持し、2.0で読み取り互換を削除する。
module bem_config_layout
  use bem_kinds, only: i32
  use bem_string_utils, only: lower_ascii
  use bem_config_toml, only: stop_config_error, require_toml_success, get_toml_int
  use bem_config_layout_paths, only: legacy_layout_path, grouped_layout_table
  use tomlf, only: toml_table, toml_array, toml_keyval, toml_value, toml_key, get_value, set_value, &
                   new_table, add_table, toml_len => len
  use tomlf_type, only: new_array
  implicit none
  private
  public :: normalize_config_layout
contains

  !> 新形式なら正規化したdocumentを返す。旧形式なら元documentを維持する。
  subroutine normalize_config_layout(document)
    type(toml_table), allocatable, intent(inout) :: document
    type(toml_table), allocatable :: normalized

    if (.not. has_grouped_keys(document, '')) return
    allocate (normalized)
    call new_table(normalized)
    call translate_table(document, '', normalized, .false.)
    call document%destroy()
    call move_alloc(normalized, document)
  end subroutine normalize_config_layout

  recursive logical function has_grouped_keys(table, prefix) result(found)
    type(toml_table), intent(inout) :: table
    character(len=*), intent(in) :: prefix
    type(toml_key), allocatable :: keys(:)
    class(toml_value), pointer :: node, item
    character(len=:), allocatable :: name, path
    integer :: i, j

    found = .false.
    call table%get_keys(keys)
    do i = 1, size(keys)
      name = lower_ascii(keys(i)%key)
      path = joined_path(prefix, name)
      select case (path)
      case ('run', 'fields', 'sheath', 'mesh.obj', 'output.enabled', 'output.history', &
            'output.final', 'output.checkpoint', 'output.diagnostics', 'particles.boundary', &
            'particles.reservoir', 'particles.tracking', 'particles.raycast', &
            'particles.species.charge_c', 'particles.species.mass_kg', 'particles.species.distribution', &
            'particles.species.source', 'particles.species.sampling', 'particles.species.inflow', &
            'particles.species.charging')
        found = .true.
        return
      end select
      call table%get(keys(i)%key, node)
      select type (node)
      type is (toml_table)
        if (path == 'particles' .or. path == 'mesh' .or. path == 'output') then
          if (has_grouped_keys(node, path)) then
            found = .true.
            return
          end if
        end if
      type is (toml_array)
        if (path /= 'particles.species') cycle
        do j = 1, toml_len(node)
          call node%get(j, item)
          select type (item)
          type is (toml_table)
            if (has_grouped_keys(item, path)) then
              found = .true.
              return
            end if
          end select
        end do
      end select
    end do
  end function has_grouped_keys

  recursive subroutine translate_table(table, prefix, target, species_context)
    type(toml_table), intent(inout) :: table
    character(len=*), intent(in) :: prefix
    type(toml_table), intent(inout), target :: target
    logical, intent(in) :: species_context
    type(toml_key), allocatable :: keys(:)
    class(toml_value), pointer :: node
    character(len=:), allocatable :: path, legacy, text
    integer :: i, stat

    if (.not. grouped_layout_table(prefix)) call stop_config_error('Unknown or mixed-layout table: '//prefix)
    if (prefix == 'run.restart') call put_logical(target, 'output.resume', .true.)
    if (prefix == 'particles.species') call validate_source_presence(table)
    if (prefix == 'particles.species.source') then
      call find_node(table, 'mode', node)
      if (.not. associated(node)) call stop_config_error('particles.species.source requires mode')
    end if
    call table%get_keys(keys)
    do i = 1, size(keys)
      path = joined_path(prefix, lower_ascii(keys(i)%key))
      call table%get(keys(i)%key, node)
      legacy = legacy_layout_path(path)
      if (path == 'particles.species') then
        select type (node)
        type is (toml_array)
          call translate_species(node, target)
        class default
          call stop_config_error('particles.species must be an array of tables')
        end select
      else if (path == 'fields.periodic.backend') then
        call get_value(table, keys(i), text, stat=stat)
        call require_toml_success(stat, path)
        if (.not. allocated(text)) call stop_config_error(path//' must be a string')
        call translate_backend(lower_ascii(trim(text)), target)
      else if (path == 'sheath.closure') then
        call get_value(table, keys(i), text, stat=stat)
        call require_toml_success(stat, path)
        if (.not. allocated(text)) call stop_config_error(path//' must be a string')
        call translate_choice(path, lower_ascii(trim(text)), target)
      else if (len(legacy) > 0) then
        if (species_context) legacy = legacy(len('particles.species.') + 1:)
        call put_node(target, legacy, node)
      else
        select type (node)
        type is (toml_table)
          call translate_table(node, path, target, species_context)
        class default
          call stop_config_error('Unknown or mixed-layout key: '//path)
        end select
      end if
    end do
    ! No split settings may silently activate a backend different from the input.
    if (prefix == 'fields.periodic') then
      call find_node(table, 'backend', node)
      text = 'finite_images'
      if (associated(node)) then
        call get_value(table, node%key, text, stat=stat)
        call require_toml_success(stat, 'fields.periodic.backend')
        text = lower_ascii(trim(text))
      end if
      if (text == 'finite_images') then
        call find_node(table, 'lower_boundary_model', node)
        if (associated(node)) call stop_config_error('split periodic settings require an explicit split backend')
        call find_node(table, 'reference', node)
        if (associated(node)) call stop_config_error('split periodic settings require an explicit split backend')
      end if
    end if
    if (prefix == 'sheath') then
      call find_node(table, 'closure', node)
      text = 'none'
      if (associated(node)) then
        call get_value(table, node%key, text, stat=stat)
        call require_toml_success(stat, 'sheath.closure')
      end if
      if (lower_ascii(trim(text)) == 'zero_current') then
        call find_node(table, 'response', node)
        if (associated(node)) then
          call get_value(table, node%key, text, stat=stat)
          call require_toml_success(stat, 'sheath.response')
          if (lower_ascii(trim(text)) /= 'zhao') call stop_config_error('zero_current requires sheath.response=zhao')
        end if
      end if
    end if
  end subroutine translate_table

  subroutine translate_species(array, target)
    type(toml_array), intent(inout) :: array
    type(toml_table), intent(inout), target :: target
    type(toml_array) :: normalized
    class(toml_value), pointer :: item
    class(toml_value), allocatable :: flat
    integer :: i, stat

    call new_array(normalized)
    do i = 1, toml_len(array)
      call array%get(i, item)
      allocate (toml_table :: flat)
      select type (flat)
      type is (toml_table)
        call new_table(flat)
        select type (item)
        type is (toml_table)
          call translate_table(item, 'particles.species', flat, .true.)
        class default
          call stop_config_error('particles.species entries must be tables')
        end select
      end select
      call normalized%push_back(flat, stat)
      call require_toml_success(stat, 'particles.species')
    end do
    call put_node(target, 'particles.species', normalized)
  end subroutine translate_species

  subroutine translate_backend(backend, target)
    character(len=*), intent(in) :: backend
    type(toml_table), intent(inout), target :: target

    select case (backend)
    case ('finite_images')
      call put_string(target, 'sim.field_periodic_far_correction', 'none')
    case ('cached_kneq0', 'panel_spectral_reference')
      if (backend == 'cached_kneq0') then
        call put_string(target, 'sim.field_periodic_far_correction', 'cached_kneq0')
      else
        call put_string(target, 'sim.field_periodic_far_correction', 'none')
      end if
      call put_string(target, 'periodic2.nonzero_mode_backend', backend)
      call put_string(target, 'periodic2.zero_mode_policy', 'exclude_k0')
    case default
      call stop_config_error('invalid fields.periodic.backend')
    end select
  end subroutine translate_backend

  subroutine translate_choice(path, value, target)
    character(len=*), intent(in) :: path, value
    type(toml_table), intent(inout), target :: target

    select case (path)
    case ('sheath.closure')
      select case (value)
      case ('none')
        call put_string(target, 'surface_current_model.model', 'none')
      case ('zero_current')
        call put_string(target, 'surface_current_model.model', 'zhao_stationary')
      case default
        call stop_config_error('invalid sheath.closure')
      end select
    end select
  end subroutine translate_choice

  subroutine validate_source_presence(table)
    type(toml_table), intent(inout) :: table
    class(toml_value), pointer :: node, count
    type(toml_key) :: key
    integer(i32) :: volume_count

    call find_node(table, 'source', node)
    if (associated(node)) return
    call find_node(table, 'sampling', node)
    if (.not. associated(node)) return
    select type (node)
    type is (toml_table)
      call find_node(node, 'volume_macro_particles_per_batch', count)
      if (.not. associated(count)) return
      key%key = count%key
      call get_toml_int(node, key, volume_count, 'particles.species.sampling.volume_macro_particles_per_batch')
      if (volume_count /= 0_i32) call stop_config_error('volume samples require particles.species.source')
    end select
  end subroutine validate_source_presence

  subroutine find_node(table, name, node)
    type(toml_table), intent(inout) :: table
    character(len=*), intent(in) :: name
    class(toml_value), pointer, intent(out) :: node
    type(toml_key), allocatable :: keys(:)
    integer :: i

    call table%get_keys(keys)
    nullify (node)
    do i = 1, size(keys)
      if (lower_ascii(keys(i)%key) /= name) cycle
      call table%get(keys(i)%key, node)
      return
    end do
  end subroutine find_node

  function joined_path(prefix, key) result(path)
    character(len=*), intent(in) :: prefix, key
    character(len=:), allocatable :: path
    path = key
    if (len(prefix) > 0) path = prefix//'.'//key
  end function joined_path

  subroutine put_node(target, path, node)
    type(toml_table), intent(inout), target :: target
    character(len=*), intent(in) :: path
    class(toml_value), intent(inout) :: node
    class(toml_value), allocatable :: copied
    type(toml_table), pointer :: parent, child
    integer :: start, finish, stat
    character(len=:), allocatable :: key

    parent => target
    start = 1
    do
      finish = index(path(start:), '.')
      if (finish == 0) exit
      key = path(start:start + finish - 2)
      call get_value(parent, key, child, requested=.false., stat=stat)
      call require_toml_success(stat, path)
      if (.not. associated(child)) then
        call add_table(parent, key, child, stat=stat)
        call require_toml_success(stat, path)
      end if
      parent => child
      start = start + finish
    end do
    call clone_node(node, copied)
    copied%key = path(start:)
    call parent%push_back(copied, stat)
    call require_toml_success(stat, 'duplicate or invalid setting: '//path)
  end subroutine put_node

  !> 再帰的なpolymorphic containerを構成要素ごとに複製し、所有権を分離する。
  recursive subroutine clone_node(node, copied)
    class(toml_value), intent(inout) :: node
    class(toml_value), allocatable, intent(out) :: copied
    class(toml_value), pointer :: child
    class(toml_value), allocatable :: copied_child
    type(toml_key), allocatable :: keys(:)
    integer :: i, stat

    select type (node)
    type is (toml_table)
      allocate (toml_table :: copied)
      select type (copied)
      type is (toml_table)
        call new_table(copied)
        call node%get_keys(keys)
        do i = 1, size(keys)
          call node%get(keys(i)%key, child)
          call clone_node(child, copied_child)
          call copied%push_back(copied_child, stat)
          call require_toml_success(stat, 'copy TOML table')
        end do
      end select
    type is (toml_array)
      allocate (toml_array :: copied)
      select type (copied)
      type is (toml_array)
        call new_array(copied)
        copied%inline = node%inline
        do i = 1, toml_len(node)
          call node%get(i, child)
          call clone_node(child, copied_child)
          call copied%push_back(copied_child, stat)
          call require_toml_success(stat, 'copy TOML array')
        end do
      end select
    type is (toml_keyval)
      allocate (copied, source=node)
    class default
      call stop_config_error('Unsupported TOML value')
    end select
    if (allocated(node%key)) copied%key = node%key
    copied%origin = node%origin
  end subroutine clone_node

  subroutine put_string(target, path, value)
    type(toml_table), intent(inout), target :: target
    character(len=*), intent(in) :: path, value
    type(toml_table) :: scratch
    class(toml_value), pointer :: node
    integer :: stat

    call new_table(scratch)
    call set_value(scratch, 'value', value, stat=stat)
    call require_toml_success(stat, path)
    call scratch%get('value', node)
    call put_node(target, path, node)
    call scratch%destroy()
  end subroutine put_string

  subroutine put_logical(target, path, value)
    type(toml_table), intent(inout), target :: target
    character(len=*), intent(in) :: path
    logical, intent(in) :: value
    type(toml_table) :: scratch
    class(toml_value), pointer :: node
    integer :: stat

    call new_table(scratch)
    call set_value(scratch, 'value', value, stat=stat)
    call require_toml_success(stat, path)
    call scratch%get('value', node)
    call put_node(target, path, node)
    call scratch%destroy()
  end subroutine put_logical
end module bem_config_layout
