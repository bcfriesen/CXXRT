include_directories (/home/brian/src//eigen-eigen-6b38706d90a9)

add_subdirectory (EOS)

set (CXXRT_SOURCES
     Lambda_iteration.cc
     calc_ALO.cc
     grid.cc
     main.cc
     planck_function.cc
     ray.cc
    )

ADD_EXECUTABLE(CXXRT ${CXXRT_SOURCES})

target_link_libraries(CXXRT yaml-cpp EOS)
