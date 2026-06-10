enable_testing()

add_executable(tests)

target_sources(tests PUBLIC
        DegenerateSimplexMethodTests.cpp
        BranchAndBoundTests.cpp
        SolutionCheckingTests.cpp
        ProblemBuilderTests.cpp
        BigIntegerTests.cpp
        NumericalStabilityTests.cpp
        BoundsTests.cpp

        simplex/FreeVariablesTests.cpp
        simplex/RandomTests.cpp
        simplex/SmallTests.cpp
        simplex/UnboundedSimplexTests.cpp

        #        LU/CommonTests.cpp
        #        LU/SparseLUTests.cpp
        #        LU/DenseLUTests.cpp
        #        LU/EtaFileTests.cpp
        #
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
        linalg/VectorTests.cpp
        linalg/CSCMatrixTests.cpp
        linalg/UseCasesTests.cpp
        linalg/lu/EtaFileTests.cpp
        linalg/lu/FullPivotingLUTests.cpp
        linalg/lu/SparseLUTests.cpp
        linalg/RowBasisTests.cpp

        linalg/expr/AllTests.cpp
        linalg/expr/SubColsExprTests.cpp
        linalg/expr/OwningViewTests.cpp
        linalg/expr/MulExprTests.cpp
        linalg/expr/SumExprTests.cpp
        linalg/expr/TransposedExprTests.cpp
)
target_link_libraries(tests GTest::gtest_main raspisator)
set_target_properties(tests PROPERTIES CMAKE_CXX_FLAGS "-fsanitize=address")
target_include_directories(tests PRIVATE .)

include(GoogleTest)
gtest_discover_tests(tests)