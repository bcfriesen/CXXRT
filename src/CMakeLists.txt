include_directories ($ENV{EIGEN})
include_directories ($ENV{YAML}/include)
include_directories ($ENV{BOOST_DIR}/include)

add_subdirectory (EOS)

set (CXXRT_SOURCES
     Lambda_iteration.cc
     calc_ALO.cc
     grid.cc
     main.cc
     planck_function.cc
     ray.cc
     rmsd.cc
    )

ADD_EXECUTABLE(CXXRT ${CXXRT_SOURCES})

target_link_libraries(CXXRT $ENV{YAML}/build/libyaml-cpp.a EOS)
