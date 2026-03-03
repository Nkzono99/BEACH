module bem_importers
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type
  use bem_mesh, only: init_mesh
  implicit none
contains

  subroutine load_obj_mesh(path, mesh)
    character(len=*), intent(in) :: path
    type(mesh_type), intent(out) :: mesh
    integer(i32) :: nvert, ntri
    real(dp), allocatable :: vertices(:, :)
    integer(i32), allocatable :: faces(:, :)

    call scan_obj(path, nvert, ntri)
    if (nvert == 0 .or. ntri == 0) error stop "OBJ has no vertices/faces"
    allocate(vertices(3, nvert), faces(3, ntri))
    call parse_obj(path, vertices, faces)
    call build_mesh_from_indexed(vertices, faces, mesh)
  end subroutine load_obj_mesh

  subroutine build_mesh_from_indexed(vertices, faces, mesh)
    real(dp), intent(in) :: vertices(:, :)
    integer(i32), intent(in) :: faces(:, :)
    type(mesh_type), intent(out) :: mesh
    real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)
    integer(i32) :: i, ntri

    if (size(vertices, 1) /= 3 .or. size(faces, 1) /= 3) then
      error stop "vertices/faces shape mismatch"
    end if

    ntri = size(faces, 2)
    allocate(v0(3, ntri), v1(3, ntri), v2(3, ntri))
    do i = 1, ntri
      v0(:, i) = vertices(:, faces(1, i))
      v1(:, i) = vertices(:, faces(2, i))
      v2(:, i) = vertices(:, faces(3, i))
    end do
    call init_mesh(mesh, v0, v1, v2)
  end subroutine build_mesh_from_indexed

  subroutine scan_obj(path, nvert, ntri)
    character(len=*), intent(in) :: path
    integer(i32), intent(out) :: nvert, ntri
    character(len=1024) :: line
    integer :: u, ios, ntok

    nvert = 0
    ntri = 0
    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop "failed to open OBJ"

    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (is_vertex_line(line)) nvert = nvert + 1
      if (is_face_line(line)) then
        ntok = count_face_tokens(line)
        if (ntok >= 3) ntri = ntri + (ntok - 2)
      end if
    end do
    close(u)
  end subroutine scan_obj

  subroutine parse_obj(path, vertices, faces)
    character(len=*), intent(in) :: path
    real(dp), intent(out) :: vertices(:, :)
    integer(i32), intent(out) :: faces(:, :)
    character(len=1024) :: line
    integer :: u, ios, i
    integer(i32) :: iv, itri, ntok, idx(512)

    iv = 0
    itri = 0
    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop "failed to open OBJ"

    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (is_vertex_line(line)) then
        iv = iv + 1
        call parse_vertex_line(line, vertices(:, iv))
      else if (is_face_line(line)) then
        call parse_face_line(line, iv, idx, ntok)
        if (ntok >= 3) then
          do i = 2, ntok - 1
            itri = itri + 1
            faces(:, itri) = [idx(1), idx(i), idx(i + 1)]
          end do
        end if
      end if
    end do
    close(u)
  end subroutine parse_obj

  logical pure function is_vertex_line(line)
    character(len=*), intent(in) :: line
    is_vertex_line = len_trim(adjustl(line)) > 1 .and. adjustl(line)(1:2) == 'v '
  end function is_vertex_line

  logical pure function is_face_line(line)
    character(len=*), intent(in) :: line
    is_face_line = len_trim(adjustl(line)) > 1 .and. adjustl(line)(1:2) == 'f '
  end function is_face_line

  integer(i32) pure function count_face_tokens(line)
    character(len=*), intent(in) :: line
    character(len=1024) :: s
    integer :: pos, n, i

    s = trim(adjustl(line))
    pos = 3
    n = len_trim(s)
    count_face_tokens = 0
    do while (pos <= n)
      do while (pos <= n .and. s(pos:pos) == ' ')
        pos = pos + 1
      end do
      if (pos > n) exit
      count_face_tokens = count_face_tokens + 1
      do i = pos, n
        if (s(i:i) == ' ') then
          pos = i + 1
          exit
        end if
        if (i == n) pos = n + 1
      end do
    end do
  end function count_face_tokens

  subroutine parse_vertex_line(line, p)
    character(len=*), intent(in) :: line
    real(dp), intent(out) :: p(3)
    character(len=1024) :: s
    s = trim(adjustl(line))
    read(s(3:), *) p(1), p(2), p(3)
  end subroutine parse_vertex_line

  subroutine parse_face_line(line, nvert, idx, ntok)
    character(len=*), intent(in) :: line
    integer(i32), intent(in) :: nvert
    integer(i32), intent(out) :: idx(:)
    integer(i32), intent(out) :: ntok
    character(len=1024) :: s, tok
    integer :: pos, n, slash, i
    integer(i32) :: vi

    s = trim(adjustl(line))
    pos = 3
    n = len_trim(s)
    ntok = 0

    do while (pos <= n)
      do while (pos <= n .and. s(pos:pos) == ' ')
        pos = pos + 1
      end do
      if (pos > n) exit
      tok = ''
      do i = pos, n
        if (s(i:i) == ' ') then
          tok = s(pos:i - 1)
          pos = i + 1
          exit
        end if
        if (i == n) then
          tok = s(pos:n)
          pos = n + 1
        end if
      end do

      slash = index(tok, '/')
      if (slash > 0) tok = tok(1:slash - 1)
      read(tok, *) vi
      if (vi < 0) vi = nvert + vi + 1
      if (vi <= 0 .or. vi > nvert) error stop "OBJ face index out of range"
      ntok = ntok + 1
      idx(ntok) = vi
    end do
  end subroutine parse_face_line

end module bem_importers
