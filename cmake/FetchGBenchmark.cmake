set(BENCHMARK_ENABLE_TESTING NO)

include(FetchContent)
FetchContent_Declare(
        googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG eddb0241389718a23a42db6af5f0164b6e0139af
)
FetchContent_MakeAvailable(googlebenchmark)