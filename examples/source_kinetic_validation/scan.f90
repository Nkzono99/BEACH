!> Normalized offline survey driver. Compile with bem_kinds/constants/source_kinetic_sheath.
program source_kinetic_scan
  use bem_kinds, only: dp
  use bem_source_kinetic_sheath
  implicit none
  type(source_kinetic_result) :: result
  type(source_kinetic_root) :: r
  real(dp) :: mach, tau, emission, field, depth
  integer :: query, points, ios, i
  print '(a)', 'query_id,status,root_count,branch,phi_h,phi_min,amplitude,q,current,field_residual,outer_residual'
  do
    read (*, *, iostat=ios) query, mach, tau, emission, field, points, depth
    if (ios < 0) exit
    if (ios /= 0) error stop 'invalid normalized survey input'
    call solve_source_kinetic(mach, tau, emission, field, 1836.15267343_dp, result, points, depth)
    if (size(result%roots) == 0) then
      print '(i0,",",i0,",0,,,,,,,,")', query, result%status
    else
      do i = 1, size(result%roots)
        r = result%roots(i)
        print '(3(i0,","),a,7(",",es24.16))', query, result%status, size(result%roots), trim(r%branch), &
          r%phi_h, r%phi_min, r%amplitude, r%escaping_flux, r%current, r%field_squared_residual, r%outer_residual
      end do
    end if
  end do
end program
