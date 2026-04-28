enable_testing()

add_executable(tests)

target_sources(tests PUBLIC
        MatrixTests.cpp
        DegenerateSimplexMethodTests.cpp
        BranchAndBoundTests.cpp
        RowBasisTests.cpp
        SolutionCheckingTests.cpp
        ProblemBuilderTests.cpp
        BigIntegerTests.cpp
        NumericalStabilityTests.cpp
        BoundsTests.cpp

        simplex/RandomTests.cpp
        simplex/SmallTests.cpp
        simplex/BigTests.cpp
        simplex/UnboundedSimplexTests.cpp

        LU/CommonTests.cpp
        LU/SparseLUTests.cpp
        LU/DenseLUTests.cpp
        LU/EtaFileTests.cpp

        accumulators/ArgMinTests.cpp
        accumulators/ArgMaxTests.cpp
        accumulators/ArithmeticMeanTests.cpp

        feature/FullProblemSolvingTests.cpp
)
target_link_libraries(tests GTest::gtest_main raspisator)
set_target_properties(tests PROPERTIES CMAKE_CXX_FLAGS "-fsanitize=address")

include(GoogleTest)
gtest_discover_tests(tests)