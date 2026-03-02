# set cover
add_executable(play_setcover setcover/setcover.cpp)
include_directories(play_setcover setcover)
target_link_libraries(play_setcover raspisator)
set_target_properties(play_setcover PROPERTIES CMAKE_CXX_FLAGS " -fsanitize=address -O0")

# test problems
add_executable(play_test_problems test_problems.cpp)
target_link_libraries(play_test_problems raspisator)