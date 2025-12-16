cmake_minimum_required(VERSION 3.31)
project(raspisator)

set(CMAKE_CXX_STANDARD 23)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(load_libraries)

add_subdirectory(src)
add_subdirectory(benchmark)
add_subdirectory(tests)

add_executable(playground playground/main.cpp)
target_link_libraries(playground raspisator)
set_target_properties(playground PROPERTIES CMAKE_CXX_FLAGS " -fsanitize=address -O0")

add_custom_target(copy_resources
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/resources
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/resources ${CMAKE_BINARY_DIR}/resources
        COMMENT "copying ${CMAKE_SOURCE_DIR}/resources to ${CMAKE_BINARY_DIR}/resources"
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_dependencies(raspisator copy_resources)