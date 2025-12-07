#include <gtest/gtest.h>

#include "../src/linear/simplex/BoundedSimplexMethod.h"
#include "linear/BigInteger.h"
#include "linear/bb/BranchAndBound.h"
#include "linear/matrix/Matrix.h"

TEST(BranchAndBoundTests, SimpleProblems) {
  // problems are from:
  // http://web.tecnico.ulisboa.pt/mcasquilho/compute/_linpro/TaylorB_module_c.pdf

  {
    Matrix<Rational> A = {{8000, 4000, 1, 0}, {15, 30, 0, 1}};

    Matrix<Rational> b = {{40'000}, {200}};
    Matrix<Rational> c = {{100, 150, 0, 0}};

    std::vector variables = {VariableType::INTEGER, VariableType::INTEGER,
                             VariableType::SLACK, VariableType::SLACK};

    MILPProblem problem(A, b, c, variables);
    BranchAndBound<Rational, simplex::BoundedSimplexMethod<Rational>> solver(problem);
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

    std::vector variables = {VariableType::INTEGER, VariableType::INTEGER,
                             VariableType::INTEGER, VariableType::INTEGER,
                             VariableType::SLACK,   VariableType::SLACK,
                             VariableType::SLACK,   VariableType::SLACK};

    MILPProblem problem(A, b, c, variables);
    BranchAndBound<Rational, simplex::BoundedSimplexMethod<Rational>> solver(problem);
    auto solution = std::get<FiniteMILPSolution<Rational>>(solver.solve());

    ASSERT_EQ(solution.value, 700);

    ASSERT_EQ((solution.point[0, 0]), 1);
    ASSERT_EQ((solution.point[1, 0]), 0);
    ASSERT_EQ((solution.point[2, 0]), 1);
    ASSERT_EQ((solution.point[3, 0]), 0);
  }
}

// TODO: next two tests must be updated when BranchAndBound distinguishes
// between InfiniteSolution and NoFiniteSolution
TEST(BranchAndBoundTests, InfiniteSolution) {
  {
    Matrix<Rational> A = {{0, 1}};

    Matrix<Rational> b = {{1}};
    Matrix<Rational> c = {{1, 0}};

    std::vector variables = {VariableType::INTEGER, VariableType::INTEGER};

    MILPProblem problem(A, b, c, variables);
    BranchAndBound<Rational, simplex::BoundedSimplexMethod<Rational>> solver(problem);
    auto solution = solver.solve();

    ASSERT_TRUE(std::holds_alternative<NoFiniteSolution>(solution));
  }

  {
    Matrix<Rational> A = {{1, -1, 1, 0}, {-1, 1, 0, 1}};

    Matrix<Rational> b = {{1}, {1}};
    Matrix<Rational> c = {{1, 0, 0, 0}};

    std::vector variables = {VariableType::INTEGER, VariableType::INTEGER,
                             VariableType::INTEGER, VariableType::INTEGER};

    MILPProblem problem(A, b, c, variables);
    BranchAndBound<Rational, simplex::BoundedSimplexMethod<Rational>> solver(problem);
    auto solution = solver.solve();

    ASSERT_TRUE(std::holds_alternative<NoFiniteSolution>(solution));
  }
}

TEST(BranchAndBoundTests, NoFeasibleElements) {
  // x_2 -> max, s.t. 1/3 <= x_1 <= 2/3
  Matrix<Rational> A = {{-1, 0, 1, 0}, {1, 0, 0, 1}};

  Matrix<Rational> b = {{-Rational{1} / 3}, {Rational{2} / 3}};
  Matrix<Rational> c = {{0, 1, 0, 0}};

  std::vector variables = {VariableType::INTEGER, VariableType::INTEGER,
                           VariableType::SLACK, VariableType::SLACK};

  MILPProblem problem(A, b, c, variables);
  BranchAndBound<Rational, simplex::BoundedSimplexMethod<Rational>> solver(problem);
  auto solution = solver.solve();

  ASSERT_TRUE(std::holds_alternative<NoFiniteSolution>(solution));
}
