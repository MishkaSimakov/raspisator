cmake_minimum_required(VERSION 3.31)
project(raspisator)

add_compile_definitions(LOG_FILEPATH="${CMAKE_CURRENT_LIST_DIR}/log")
add_compile_definitions(RESOURCES_FILEPATH="${CMAKE_CURRENT_LIST_DIR}/resources")

set(CMAKE_CXX_STANDARD 23)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

include(load_libraries)

add_subdirectory(src)
add_subdirectory(benchmark)
add_subdirectory(tests)
add_subdirectory(playground)
