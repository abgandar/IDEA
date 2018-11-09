# IDEA
#

cmake_minimum_required (VERSION 3.0.0)

project(IDEA CXX)

# DACE requires at least C++11
set(CMAKE_CXX_STANDARD 11)

# Set this to your DACE installation directory if it is not the system directory
#set(dace_DIR ".../root/usr/local/lib/cmake/dace")
find_package(dace 2.0.0 REQUIRED)

# globally set the include directory for library
include_directories(src/include)

# now process the subdirectories
add_subdirectory(src)
