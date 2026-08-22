!> ASCII 文字列操作ユーティリティ。
module bem_string_utils
  implicit none
  private

  public :: lower_ascii
  public :: is_decimal_real_token, is_decimal_integer_token, is_logical_token

contains

  !> ASCII 英字を小文字化した文字列を返す。
  pure function lower_ascii(s) result(out)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: out
    integer :: i, code

    out = s
    do i = 1, len(s)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        out(i:i) = achar(code + 32)
      end if
    end do
  end function lower_ascii

  !> List-directed control syntaxを除外した十進実数tokenだけを許可する。
  pure logical function is_decimal_real_token(token) result(valid)
    character(len=*), intent(in) :: token
    integer :: position, token_length
    logical :: has_mantissa_digit, has_exponent_digit

    valid = .false.
    token_length = len_trim(token)
    if (token_length <= 0) return
    position = 1
    if (token(position:position) == '+' .or. token(position:position) == '-') position = position + 1
    if (position > token_length) return

    has_mantissa_digit = .false.
    do while (position <= token_length)
      if (.not. is_decimal_digit(token(position:position))) exit
      has_mantissa_digit = .true.
      position = position + 1
    end do
    if (position <= token_length) then
      if (token(position:position) == '.') then
        position = position + 1
        do while (position <= token_length)
          if (.not. is_decimal_digit(token(position:position))) exit
          has_mantissa_digit = .true.
          position = position + 1
        end do
      end if
    end if
    if (.not. has_mantissa_digit) return

    if (position <= token_length) then
      if (index('eEdD', token(position:position)) > 0) then
        position = position + 1
        if (position <= token_length) then
          if (token(position:position) == '+' .or. token(position:position) == '-') position = position + 1
        end if
        has_exponent_digit = .false.
        do while (position <= token_length)
          if (.not. is_decimal_digit(token(position:position))) exit
          has_exponent_digit = .true.
          position = position + 1
        end do
        if (.not. has_exponent_digit) return
      end if
    end if
    valid = position > token_length
  end function is_decimal_real_token

  !> 符号付き十進整数tokenだけを許可する。
  pure logical function is_decimal_integer_token(token) result(valid)
    character(len=*), intent(in) :: token
    integer :: position, token_length

    valid = .false.
    token_length = len_trim(token)
    if (token_length <= 0) return
    position = 1
    if (token(position:position) == '+' .or. token(position:position) == '-') position = position + 1
    if (position > token_length) return
    do while (position <= token_length)
      if (.not. is_decimal_digit(token(position:position))) return
      position = position + 1
    end do
    valid = .true.
  end function is_decimal_integer_token

  !> Checkpoint writerが生成するlogical tokenと標準的な長形式だけを許可する。
  pure logical function is_logical_token(token) result(valid)
    character(len=*), intent(in) :: token
    character(len=len(token)) :: normalized

    normalized = lower_ascii(trim(token))
    valid = trim(normalized) == 't' .or. trim(normalized) == 'f' .or. &
            trim(normalized) == '.true.' .or. trim(normalized) == '.false.'
  end function is_logical_token

  pure logical function is_decimal_digit(character_value) result(is_digit)
    character(len=1), intent(in) :: character_value

    is_digit = character_value >= '0' .and. character_value <= '9'
  end function is_decimal_digit

end module bem_string_utils
