add_executable(benchmark_linalg linear/LinalgBenchmark.cpp)
target_link_libraries(benchmark_linalg benchmark::benchmark raspisator)
target_compile_options(benchmark_linalg PRIVATE "-O3")

add_executable(benchmark_simplex_primal simplex/primal.cpp)
target_link_libraries(benchmark_simplex_primal raspisator)
target_compile_options(benchmark_simplex_primal PRIVATE "-O3")