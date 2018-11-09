# IDEA library
FILE(GLOB Sources lib/*.cpp)
add_library(idea SHARED ${Sources})
target_include_directories(idea INTERFACE include)
target_link_libraries(idea PUBLIC dace::dace)

# IDEA test program
add_executable(IDEATest test/ideatest.cpp)
target_link_libraries(IDEATest PUBLIC idea)
