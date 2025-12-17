#include <gtest/gtest.h>

#include "linear/BigInteger.h"
#include "linear/bb/FullStrongBranching.h"
#include "linear/matrix/Matrix.h"
#include "linear/model/MILP.h"
#include "linear/problem/VariableType.h"
#include "linear/simplex/BoundedSimplexMethod.h"

TEST(BranchAndBoundTests, SimpleProblems) {
  // problems are from:
  // http://web.tecnico.ulisboa.pt/mcasquilho/compute/_linpro/TaylorB_module_c.pdf

  {
    Matrix<Rational> A = {{8000, 4000, 1, 0}, {15, 30, 0, 1}};

    Matrix<Rational> b = {{40'000}, {200}};
    Matrix<Rational> c = {{100, 150, 0, 0}};

    std::vector variables = {VariableType::INTEGER, VariableType::INTEGER,
                             VariableType::SLACK, VariableType::SLACK};
    std::vector<Rational> lower(4, 0);
    std::vector<Rational> upper(4, 100'000);

    FullStrongBranchingBranchAndBound solver(A, b, c, lower, upper, variables);
    auto solution =
        std::get<FiniteMILPSolution<Rational>>(solver.solve().solution);

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

    std::vector variables = {
        VariableType::INTEGER, VariableType::INTEGER, VariableType::INTEGER,
        VariableType::INTEGER, VariableType::SLACK,   VariableType::SLACK,
        VariableType::SLACK,   VariableType::SLACK,   VariableType::SLACK,
        VariableType::SLACK,   VariableType::SLACK};

    std::vector<Rational> lower(11, 0);
    std::vector<Rational> upper(11, 100'000);

    FullStrongBranchingBranchAndBound solver(A, b, c, lower, upper, variables);
    auto solution = std::get<FiniteMILPSolution<Rational>>(solver.solve().solution);

    ASSERT_EQ(solution.value, 700);

    ASSERT_EQ((solution.point[0, 0]), 1);
    ASSERT_EQ((solution.point[1, 0]), 0);
    ASSERT_EQ((solution.point[2, 0]), 1);
    ASSERT_EQ((solution.point[3, 0]), 0);
  }
}

TEST(BranchAndBoundTests, NoFeasibleElements) {
  // x_2 -> max, s.t. 1/3 <= x_1 <= 2/3
  Matrix<Rational> A = {{-1, 0, 1, 0}, {1, 0, 0, 1}};

  Matrix<Rational> b = {{-Rational{1} / 3}, {Rational{2} / 3}};
  Matrix<Rational> c = {{0, 1, 0, 0}};

  std::vector variables = {VariableType::INTEGER, VariableType::INTEGER,
                           VariableType::SLACK, VariableType::SLACK};

  std::vector<Rational> lower(4, 0);
  std::vector<Rational> upper(4, 100);

  FullStrongBranchingBranchAndBound solver(A, b, c, lower, upper, variables);
  auto run_result = solver.solve();

  ASSERT_TRUE(std::holds_alternative<NoFiniteSolution>(run_result.solution));
}
