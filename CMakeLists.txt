cmake_minimum_required(VERSION 3.13)
project(wdm_coupling VERSION 0.0.1 LANGUAGES CXX)

list(INSERT CMAKE_MODULE_PATH 0 ${wdm_coupling_SOURCE_DIR}/CMake)
if(CMAKE_CXX_COMPILER_ID MATCHES "PGI")
include(find_dependencies_summit)
endif()
enable_language(CXX)
set(MPI_CXX_SKIP_MPICXX ON)
find_package(MPI REQUIRED)
#adios2 adds C and Fortran depending on how it was built
find_package(ADIOS2 2.5 REQUIRED)
if(CMAKE_CXX_COMPILER_ID MATCHES "PGI")
#find_package(Kokkos MODULE REQUIRED)
else()
find_package(Kokkos REQUIRED)
endif()

## use pkgconfig since the fftw autoconf install produces
## broken cmake config files
## https://github.com/FFTW/fftw3/issues/130
if(CMAKE_CXX_COMPILER_ID MATCHES "PGI")
find_package(FFTW3 REQUIRED fftw3)
else()
find_package(PkgConfig REQUIRED)
pkg_check_modules(fftw REQUIRED IMPORTED_TARGET fftw3>=3.3)
endif()

include(CTest)
if(BUILD_TESTING)
  #TODO update the following logic after merging Kai's find_package(mpi) support
  get_filename_component(COMPILER_DIR "${CMAKE_CXX_COMPILER}" PATH) 
  find_program(MPIRUN NAMES mpirun PATHS "${COMPILER_DIR}")
  message(STATUS "MPIRUN: ${MPIRUN}")
  set(MPIRUN_PROCFLAG "-np" CACHE STRING
      "the command line flag to give process count to MPIRUN")
  message(STATUS "MPIRUN_PROCFLAG: ${MPIRUN_PROCFLAG}")
  set(TEST_DATA_DIR "" CACHE PATH
    "Path to a local copy of the wdmapp_coupling_data repo.")
  if(NOT EXISTS ${TEST_DATA_DIR})
    message(FATAL_ERROR "TEST_DATA_DIR \"${TEST_DATA_DIR}\" is not accessible")
  endif()
endif()

add_subdirectory(src)
add_subdirectory(test)
