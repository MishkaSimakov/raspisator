enable_testing()

add_executable(tests)

target_sources(tests PUBLIC
        SimplexMethodTests.cpp
        MatrixTests.cpp
        DegenerateSimplexMethodTests.cpp
        BranchAndBoundTests.cpp
        RowBasisTests.cpp
        ProblemBuilderTests.cpp
        BigIntegerTests.cpp
)
target_link_libraries(tests GTest::gtest_main raspisator)
set_target_properties(tests PROPERTIES CMAKE_CXX_FLAGS "-fsanitize=address")

include(GoogleTest)
gtest_discover_tests(tests)