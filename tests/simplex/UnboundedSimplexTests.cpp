#include <gtest/gtest.h>

#include <random>

#include "Assertions.h"
#include "linear/BigInteger.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/Random.h"
#include "linear/simplex/Simplex.h"

TEST(UnboundedSimplexTests, SimpleTest1) {
  CSCMatrix<Rational> A = {
      {1, -1, 1, 0},
      {2, 1, 0, 1},
  };

  Matrix<Rational> b = {{1}, {3}};
  Matrix<Rational> c = {{2, 1, 1, -1}};

  Bounds<Rational> bounds(4);
  bounds[0] = {0, 10};
  bounds[1] = {std::nullopt, 10};
  bounds[2] = {0, 10};
  bounds[3] = {0, std::nullopt};

  simplex::Simplex solver(A, b, c);

  auto feasible = solver.try_get_primal_feasible(bounds);

  ASSERT_TRUE(feasible.has_value());

  auto solution = solver.primal(bounds, *feasible).solution;

  Matrix<Rational> expected = {{0}, {3}, {4}, {0}};

  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).point, expected);
  ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).value, 7);
}

TEST(UnboundedSimplexTests, SimpleTest2) {
  Matrix<Rational> A = {
      {2, 1, 1, 0, 0},
      {3, 5, 0, 1, 0},
      {1, 3, 0, 0, 1},
  };

  Matrix<Rational> b = {{3}, {9}, {5}};
  Matrix<Rational> c = {{1, 4, 0, 0, 0}};

  Bounds<Rational> bounds(5);
  for (size_t i = 0; i < 5; ++i) {
    bounds[i] = {0, std::nullopt};
  }

  simplex::Simplex solver(CSCMatrix(A), b, c);

  auto feasible = solver.try_get_primal_feasible(bounds);
  ASSERT_TRUE(feasible.has_value());

  auto solution = solver.primal(bounds, *feasible).solution;

  ASSERT_TRUE(std::holds_alternative<FiniteLPSolution<Rational>>(solution));

  auto finite = std::get<FiniteLPSolution<Rational>>(solution);

  ASSERT_EQ(finite.value, Rational{20} / Rational{3});
  ASSERT_EQ((finite.point[0, 0]), 0);
  ASSERT_EQ((finite.point[1, 0]), Rational{5} / Rational{3});
}

TEST(UnboundedSimplexTests, UnboundedTest) {
  Matrix<Rational> A = {
      {1, -1, 1, 0},
      {-2, 1, 0, 1},
  };

  Matrix<Rational> b = {{1}, {2}};
  Matrix<Rational> c = {{3, 4, 0, 0}};

  Bounds<Rational> bounds(4);
  bounds[0] = {0, std::nullopt};
  bounds[1] = {0, std::nullopt};
  bounds[2] = {0, std::nullopt};
  bounds[3] = {0, std::nullopt};

  simplex::Simplex solver(CSCMatrix(A), b, c);

  auto feasible = solver.try_get_primal_feasible(bounds);
  ASSERT_TRUE(feasible.has_value());

  auto solution = solver.primal(bounds, *feasible).solution;

  ASSERT_TRUE(std::holds_alternative<Unbounded>(solution));
}
