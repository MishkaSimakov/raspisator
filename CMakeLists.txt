cmake_minimum_required(VERSION 3.31)
project(raspisator)

set(CMAKE_CXX_STANDARD 23)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(load_libraries)

add_subdirectory(src)
add_subdirectory(tests)

add_executable(playground playground/main.cpp)
target_link_libraries(playground raspisator)