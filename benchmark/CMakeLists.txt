add_executable(benchmark_linalg linear/LinalgBenchmark.cpp)
target_link_libraries(benchmark_linalg benchmark::benchmark raspisator)
target_compile_options(benchmark_linalg PRIVATE "-O3")