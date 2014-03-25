set (CXXRT_SOURCES
     main.cc
     planck_function.cc
     ray.cc
    )

ADD_EXECUTABLE(CXXRT ${CXXRT_SOURCES})

target_link_libraries(CXXRT yaml-cpp)
