include_directories (/home/brian/src//eigen-eigen-6b38706d90a9)

set (CXXRT_SOURCES
     calc_ALO.cc
     grid.cc
     main.cc
     planck_function.cc
     ray.cc
    )

ADD_EXECUTABLE(CXXRT ${CXXRT_SOURCES})

target_link_libraries(CXXRT yaml-cpp)
