#include <gtest/gtest.h>

#include "linear/BigInteger.h"
#include "linear/BranchAndBound.h"
#include "linear/Matrix.h"
#include "linear/SimplexMethod.h"

TEST(BranchAndBoundTests, SimpleProblems) {
  // problems are from:
  // http://web.tecnico.ulisboa.pt/mcasquilho/compute/_linpro/TaylorB_module_c.pdf

  {
    Matrix<Rational> A = {{8000, 4000, 1, 0}, {15, 30, 0, 1}};

    Matrix<Rational> b = {{40'000}, {200}};
    Matrix<Rational> c = {{100, 150, 0, 0}};

    std::vector<size_t> integer = {0, 1};

    MILPProblem problem(A, b, c, integer);
    BranchAndBound<Rational, SimplexMethod<Rational>> solver(problem);
    auto solution = std::get<FiniteMILPSolution<Rational>>(solver.solve());

    ASSERT_EQ(solution.value, 1000);
    ASSERT_EQ((solution.point[0, 0]), 1);
    ASSERT_EQ((solution.point[1, 0]), 6);
  }

  {
    Matrix<Rational> A = {
        {35'000, 10'000, 25'000, 90'000, 1, 0, 0, 0, 0, 0, 0},
        {4, 2, 7, 3, 0, 1, 0, 0, 0, 0, 0},
        {1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
    };

    Matrix<Rational> b = {{120'000}, {12}, {1}, {1}, {1}, {1}, {1}};
    Matrix<Rational> c = {{300, 90, 400, 150, 0, 0, 0, 0, 0, 0, 0}};

    std::vector<size_t> integer = {0, 1, 2, 3};

    MILPProblem problem(A, b, c, integer);
    BranchAndBound<Rational, SimplexMethod<Rational>> solver(problem);
    auto solution = std::get<FiniteMILPSolution<Rational>>(solver.solve());

    ASSERT_EQ(solution.value, 700);

    ASSERT_EQ((solution.point[0, 0]), 1);
    ASSERT_EQ((solution.point[1, 0]), 0);
    ASSERT_EQ((solution.point[2, 0]), 1);
    ASSERT_EQ((solution.point[3, 0]), 0);
  }
}
