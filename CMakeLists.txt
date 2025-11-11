cmake_minimum_required(VERSION 3.31)
project(raspisator)

set(CMAKE_CXX_STANDARD 23)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(load_libraries)

add_subdirectory(tests)

add_executable(raspisator main.cpp)