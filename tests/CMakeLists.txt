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

        simplex/FreeVariablesTests.cpp
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

        mps/DataRecordTokenizerTests.cpp
        mps/IndicatorRecordTokenizerTests.cpp
        mps/RowsParserTests.cpp
        mps/ColumnsParserTests.cpp
        mps/RHSParserTests.cpp
        mps/RangesParserTests.cpp
        mps/BoundsParserTests.cpp
        mps/MPSParserTests.cpp
        mps/IntegerVariablesExtensionTests.cpp
        mps/MPSTests.cpp

        feature/FullProblemSolvingTests.cpp

        linalg/MatrixInitializationTests.cpp
        linalg/MatrixResizeTests.cpp
)
target_link_libraries(tests GTest::gtest_main raspisator)
set_target_properties(tests PROPERTIES CMAKE_CXX_FLAGS "-fsanitize=address")

include(GoogleTest)
gtest_discover_tests(tests)