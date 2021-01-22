# Distributed under the OSI-approved BSD 3-Clause License.

#.rst:
# FindFFTW3
# ---------
#
# Try to find FFTW3
#
# Uses FFTW3_ROOT in the cache variables or in the environment as a hint
# where to search
#
# IMPORTED Targets
# ^^^^^^^^^^^^^^^^
#
# This module defines :prop_tgt:`IMPORTED` target ``FFTW3::FFTW3``, if
# FFTW3 has been found.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module defines the following variables:
#
# ``FFTW3_FOUND``
#   system has FFTW3
# ``FFTW3_INCLUDE_DIRS``
#   the FFTW3 include directories
# ``FFTW3_LIBRARIES``
#   Link these to use FFTW3

set(_hints "")
if(FFTW3_ROOT)
  set(_hints "${FFTW3_ROOT}")
else()
  if(NOT ("$ENV{FFTW3_ROOT}" STREQUAL ""))
    set(_hints "$ENV{FFTW3_ROOT}")
    # module cray-fftw sets FFTW_DIR, which includes trailing .../lib
  elseif(NOT ("$ENV{FFTW_DIR}" STREQUAL ""))
    # drop trailing .../lib
    get_filename_component(_hints "$ENV{FFTW_DIR}" DIRECTORY)
  endif()
endif()

find_path(FFTW3_INCLUDE_DIRS NAMES fftw3.h HINTS ${_hints} PATH_SUFFIXES include)
find_library(FFTW3_LIBRARIES NAMES fftw3 HINTS ${_hints} PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3
                                  REQUIRED_VARS FFTW3_LIBRARIES FFTW3_INCLUDE_DIRS)

if (FFTW3_FOUND)
    if(NOT TARGET FFTW3::FFTW3)
      add_library(FFTW3::FFTW3 INTERFACE IMPORTED)
      set_property(TARGET FFTW3::FFTW3 PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}")
      set_property(TARGET FFTW3::FFTW3 PROPERTY
        INTERFACE_LINK_LIBRARIES "${FFTW3_LIBRARIES}")
    endif()
endif()
