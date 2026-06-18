enable_testing()

add_executable(tests)

target_sources(tests PUBLIC
        utils/MergeSortedRangeTests.cpp

        field/BigIntegerTests.cpp

        schedule/SolutionCheckingTests.cpp
#        schedule/FullProblemSolvingTests.cpp

        problem/ProblemBuilderTests.cpp
        problem/BoundsTests.cpp

        simplex/Phase1Tests.cpp
        simplex/CatalogTests.cpp
        simplex/RandomTests.cpp

        accumulators/ArgMinTests.cpp
        accumulators/ArgMaxTests.cpp
        accumulators/ArithmeticMeanTests.cpp
        accumulators/GeometricMeanTests.cpp

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

        presolve/RemoveLinearlyDependentEqualitiesTests.cpp
        presolve/TransformToEqualitiesTests.cpp
        presolve/ScalingTests.cpp

        linalg/MatrixInitializationTests.cpp
        linalg/MatrixResizeTests.cpp
        linalg/VectorTests.cpp
        linalg/CSCMatrixTests.cpp
        linalg/PermutationTests.cpp
        linalg/UseCasesTests.cpp
        linalg/RowBasisTests.cpp
        linalg/RRQRTests.cpp
        linalg/SparseRRQRTests.cpp

        linalg/lu/EtaFileTests.cpp
        linalg/lu/FullPivotingLUTests.cpp
        linalg/lu/LUPATests.cpp

        linalg/expr/AllTests.cpp
        linalg/expr/SubColsExprTests.cpp
        linalg/expr/OwningViewTests.cpp
        linalg/expr/MulExprTests.cpp
        linalg/expr/SumExprTests.cpp
        linalg/expr/TransposedExprTests.cpp

        simplex2/SmallTests.cpp
        simplex2/RandomTests.cpp
        simplex2/TraversalTests.cpp
)
target_link_libraries(tests GTest::gtest_main highs raspisator)
#target_compile_options(tests PRIVATE -fsanitize=address)
#target_link_options(tests PRIVATE -fsanitize=address)
target_include_directories(tests PUBLIC .)

include(GoogleTest)
gtest_discover_tests(tests)