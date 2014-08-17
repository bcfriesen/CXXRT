include_directories ($ENV{EIGEN})
include_directories ($ENV{YAML}/include)
include_directories ($ENV{BOOST_DIR}/include)
include_directories ($ENV{GSL_INCLUDE_DIR})

add_subdirectory (EOS)
add_subdirectory (lines)
add_subdirectory (misc)

set (CXXRT_SOURCES
     ALI.cc
     Lambda_iteration.cc
     build_internal_model.cc
     calc_ALO.cc
     grid.cc
     initialize_rays.cc
     main.cc
     planck_function.cc
     ray.cc
     read_mesa_model.cc
     rmsd.cc
     wavelength_grid.cc
    )

ADD_EXECUTABLE(CXXRT ${CXXRT_SOURCES})

target_link_libraries(CXXRT $ENV{YAML}/build/libyaml-cpp.a EOS lines $ENV{GSL_LIB_DIR}/libgsl.a $ENV{GSL_LIB_DIR}/libgslcblas.a misc)
