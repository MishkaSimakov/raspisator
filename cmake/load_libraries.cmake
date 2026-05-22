include(FetchContent)

# Google Test
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/52eb8108c5bdec04579160ae17225d66034bd723.zip
)
FetchContent_MakeAvailable(googletest)


# Google Benchmark
set(BENCHMARK_ENABLE_TESTING NO)

include(FetchContent)
FetchContent_Declare(
        googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG eddb0241389718a23a42db6af5f0164b6e0139af
)
FetchContent_MakeAvailable(googlebenchmark)


# HiGHS MILP-solver for testing
include(FetchContent)

FetchContent_Declare(
        highs
        GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
        GIT_TAG v1.7.2
)

FetchContent_MakeAvailable(highs)

