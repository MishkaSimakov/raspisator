add_executable(benchmark_lu linear/LU.cpp)
add_executable(benchmark_csc_constructor linear/CSCMatrixConstructor.cpp)
add_executable(benchmark_linalg linear/LinalgBenchmark.cpp)

target_link_libraries(benchmark_lu benchmark::benchmark raspisator)
target_link_libraries(benchmark_csc_constructor benchmark::benchmark raspisator)
target_link_libraries(benchmark_linalg benchmark::benchmark raspisator)
set(CMAKE_CXX_FLAGS "-O3")