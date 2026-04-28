function(add_playground_target target_name sources)
    add_executable(raspisator_${target_name})

    target_sources(raspisator_${target_name} PUBLIC ${sources})
    target_link_libraries(raspisator_${target_name} raspisator)

    set_target_properties(raspisator_${target_name} PROPERTIES CMAKE_CXX_FLAGS " -fsanitize=address -O0")
endfunction()

# set cover
add_executable(play_setcover setcover/setcover.cpp)
include_directories(play_setcover setcover)
target_link_libraries(play_setcover raspisator)
set_target_properties(play_setcover PROPERTIES CMAKE_CXX_FLAGS " -fsanitize=address -O0")

# test problems
add_playground_target(test_problems test_problems.cpp)

# MPS to NPY
add_playground_target(mps_to_npy mps_to_numpy.cpp)

# LUP analysis
add_playground_target(lup_analysis lup_error_measurement.cpp)