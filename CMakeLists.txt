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

add_custom_target(copy_resources
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/resources
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/resources ${CMAKE_CURRENT_BINARY_DIR}/resources
        COMMENT "copying ${CMAKE_CURRENT_SOURCE_DIR}/resources to ${CMAKE_CURRENT_BINARY_DIR}/resources"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_dependencies(raspisator copy_resources)