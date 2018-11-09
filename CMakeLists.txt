# DACE Examples
#
# This file shows how to use the DACE in a cmake project
# Simply import the DACE cmake package. It automatically provides the correct
# locations for headers and libraries.
# If the DACE is installed in a non-default location, tell cmake where to look:
#   set(dace_DIR "/my/custom/root/usr/local/lib/cmake/dace")
# It exposes two targets:
#   dace::dace   (the dynamic DACE library)
#   dace::dace_s (the static DACE library)
# Make sure to link your executable with one of those and include the correct
# header (<dace/dace.h> or <dace/dace_s.h>) in the code.

cmake_minimum_required (VERSION 2.8.4)

project(Examples2 CXX)

# DACE requires C++11
set(CMAKE_CXX_STANDARD 11)

set(dace_DIR "~/Desktop/DACE/build/root/usr/local/lib/cmake/dace")
find_package(dace 2.0.0 REQUIRED)

add_executable(Test test.cpp)
target_link_libraries(Test PUBLIC dace::dace)

add_executable(Euler1D Euler1D.cpp)
target_link_libraries(Euler1D PUBLIC dace::dace)

add_executable(EulerND EulerND.cpp)
target_link_libraries(EulerND PUBLIC dace::dace)

add_executable(EulerNDplot EulerNDplot.cpp)
target_link_libraries(EulerNDplot PUBLIC dace::dace)
