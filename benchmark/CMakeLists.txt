add_executable(BenchmarkLU linear/LU.cpp)

target_link_libraries(BenchmarkLU benchmark::benchmark raspisator)