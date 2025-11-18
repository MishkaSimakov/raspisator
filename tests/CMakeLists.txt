enable_testing()

add_executable(tests)

target_sources(tests PUBLIC
        SimplexMethodTests.cpp
        MatrixTests.cpp
        DegenerateSimplexMethodTests.cpp
        BranchAndBoundTests.cpp
        RowBasisTests.cpp
)
target_link_libraries(tests GTest::gtest_main raspisator)

include(GoogleTest)
gtest_discover_tests(tests)