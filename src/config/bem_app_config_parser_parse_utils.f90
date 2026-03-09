submodule (bem_app_config_parser) bem_app_config_parser_parse_utils
  implicit none
contains

  module procedure split_key_value
    integer :: p

    p = index(line, '=')
    if (p <= 0) then
      key = ''
      value = ''
      return
    end if
    key = lower(trim(adjustl(line(:p - 1))))
    value = trim(adjustl(line(p + 1:)))
  end procedure split_key_value

  module procedure parse_real
    read (text, *) out
  end procedure parse_real

  module procedure parse_int
    read (text, *) out
  end procedure parse_int

  module procedure parse_logical
    character(len=32) :: t

    t = lower(trim(adjustl(text)))
    select case (t)
    case ('true', '.true.')
      out = .true.
    case ('false', '.false.')
      out = .false.
    case default
      error stop 'Logical value must be true/false (or .true./.false.).'
    end select
  end procedure parse_logical

  module procedure parse_string
    character(len=:), allocatable :: tmp

    tmp = trim(adjustl(text))
    if (len(tmp) >= 2 .and. tmp(1:1) == '"' .and. tmp(len(tmp):len(tmp)) == '"') then
      tmp = tmp(2:len(tmp) - 1)
    end if
    out = trim(tmp)
  end procedure parse_string

  module procedure parse_real3
    character(len=256) :: t

    t = trim(adjustl(text))
    if (t(1:1) == '[') t = t(2:)
    if (t(len_trim(t):len_trim(t)) == ']') t = t(:len_trim(t) - 1)
    read (t, *) out(1), out(2), out(3)
  end procedure parse_real3

  module procedure parse_boundary_mode
    character(len=64) :: mode

    call parse_string(text, mode)
    select case (trim(lower(mode)))
    case ('open', 'outflow', 'escape')
      out = bc_open
    case ('reflect', 'reflection')
      out = bc_reflect
    case ('periodic')
      out = bc_periodic
    case default
      error stop 'Unknown boundary condition mode in [sim].'
    end select
  end procedure parse_boundary_mode

  module procedure strip_comment
    integer :: p

    p = index(line, '#')
    if (p > 0) then
      out = trim(line(:p - 1))
    else
      out = trim(line)
    end if
  end procedure strip_comment

  module procedure lower
    integer :: i, c

    o = s
    do i = 1, len(s)
      c = iachar(s(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) o(i:i) = achar(c + 32)
    end do
  end procedure lower

  module procedure ends_with
    integer :: ls, lf

    ls = len_trim(s)
    lf = len_trim(suffix)
    if (lf > ls) then
      ends_it = .false.
    else
      ends_it = (s(ls - lf + 1:ls) == suffix(1:lf))
    end if
  end procedure ends_with

end submodule bem_app_config_parser_parse_utils
