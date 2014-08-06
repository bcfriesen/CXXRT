include_directories ($ENV{EIGEN})
include_directories ($ENV{YAML}/include)
include_directories ($ENV{BOOST_DIR}/include)
include_directories ($ENV{GSL_DIR}/include)

add_subdirectory (EOS)
add_subdirectory (lines)

set (CXXRT_SOURCES
     Lambda_iteration.cc
     calc_ALO.cc
     grid.cc
     main.cc
     planck_function.cc
     ray.cc
     rmsd.cc
     wavelength_grid.cc
    )

ADD_EXECUTABLE(CXXRT ${CXXRT_SOURCES})

target_link_libraries(CXXRT $ENV{YAML}/build/libyaml-cpp.a EOS lines $ENV{GSL_DIR}/lib/libgsl.a $ENV{GSL_DIR}/lib/libgslcblas.a)
