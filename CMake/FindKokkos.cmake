find_package(Kokkos CONFIG)
#if(Kokkos_FOUND)
#  # Kokkos' imported target doesn't give us its include path,
#  # so we add it here.
#  get_target_property(Kokkos_LIBRARY Kokkos LOCATION)
#  get_filename_component(Kokkos_LIBRARY_DIR ${Kokkos_LIBRARY} DIRECTORY)
#  get_filename_component(Kokkos_INSTALL_DIR ${Kokkos_LIBRARY_DIR} DIRECTORY)
#  find_path(Kokkos_INCLUDE_DIR NAMES KokkosCore_config.h PATHS ${Kokkos_INSTALL_DIR}/include)
#  set_target_properties(kokkos PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Kokkos_INCLUDE_DIR}")
#endif()
