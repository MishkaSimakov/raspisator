add_executable(BenchmarkSparseLU linear/SparseLU.cpp)

target_link_libraries(BenchmarkSparseLU benchmark::benchmark raspisator)