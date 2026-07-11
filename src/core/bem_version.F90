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

#ifndef __BEACH_SOURCE_COMMIT__
#define __BEACH_SOURCE_COMMIT__ 'unknown'
#endif

#ifndef __BEACH_BUILD_ID__
#define __BEACH_BUILD_ID__ 'unknown'
#endif

  character(len=*), parameter, public :: beach_version = __BEACH_VERSION__
  character(len=*), parameter, public :: beach_version_mode = __BEACH_VERSION_MODE__
  character(len=*), parameter, public :: beach_source_commit = __BEACH_SOURCE_COMMIT__
  character(len=*), parameter, public :: beach_build_id = __BEACH_BUILD_ID__
  character(len=*), parameter, public :: beach_build_info = &
                                         'build_info_schema_version=1'//new_line('a')// &
                                         'build_version='//beach_version//new_line('a')// &
                                         'build_version_mode='//beach_version_mode//new_line('a')// &
                                         'build_source_commit='//beach_source_commit//new_line('a')// &
                                         'build_id='//beach_build_id
end module bem_version
