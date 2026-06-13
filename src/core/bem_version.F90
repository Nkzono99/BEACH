!> BEACH build metadata stamped by build.sh.
module bem_version
  implicit none
  private

#ifndef __BEACH_VERSION__
#define __BEACH_VERSION__ 'unknown'
#endif

#ifndef __BEACH_VERSION_MODE__
#define __BEACH_VERSION_MODE__ 'unknown'
#endif

  character(len=*), parameter, public :: beach_version = __BEACH_VERSION__
  character(len=*), parameter, public :: beach_version_mode = __BEACH_VERSION_MODE__
end module bem_version
