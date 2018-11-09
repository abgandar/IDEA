# IDEA library
FILE(GLOB Sources lib/*.cpp)
add_library(idea SHARED ${Sources})
target_include_directories(idea INTERFACE include)

# IDEA test program
add_executable(IDEATest test/ideatest.cpp)
target_link_libraries(IDEATest PUBLIC idea)
