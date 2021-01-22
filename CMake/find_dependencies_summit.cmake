#set(Kokkos_ROOT "/ccs/home/scheinberg/Software/Jan20/kokkos/cmake_install/lib/CMake/Kokkos")
#set(Kokkos_ROOT "/gpfs/alpine/world-shared/phy122/lib/install/summit/kokkos/pgi19.10/install/lib64/cmake/Kokkos")
#SET(Kokkos_INCLUDE_DIRS_RET "/gpfs/alpine/world-shared/phy122/lib/install/summit/kokkos/pgi19.10/install/include")
#include_directories(${Kokkos_INCLUDE_DIRS_RET})

set(ADIOS1_ROOT "/ccs/home/shku/Software/install/summit/adios/devel/pgi19.4")
#set(ADIOS2_ROOT "/gpfs/alpine/world-shared/csc143/jyc/summit/sw/adios2/devel/pgi")
#set(ADIOS2_ROOT "/gpfs/alpine/csc143/world-shared/jyc/summit/sw/adios2/v2.6.0-367-g15b8642/pgi19.9")
set(ADIOS2_ROOT /gpfs/alpine/fus123/world-shared/shuangxiz/adios2-pgi-19.9)

#set(FFTW3_ROOT "/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/pgi-19.9/fftw-3.3.8-wri566ow7qlimailhyznou5cfm4h362z")

set(LAPACK_ROOT "$ENV{OLCF_PGI_ROOT}/linuxpower/19.9")
