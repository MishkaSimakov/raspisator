#include <gtest/gtest.h>

#include <random>
#include <set>

#include "ConstructSparse.h"
#include "PassAssertions.h"
#include "presolve/passes/TransformToEqualities.h"
#include "support/Highs.h"
#include "support/ProblemConstructors.h"
#include "support/RandomProblem.h"

TEST(TransformToEqualities, AddsSlackForInequalities) {
  auto matrix = sparse<Rational>({
      {1, 2, 3, 4},
      {0, 4, 1, 2},
      {1, -2, 2, 2},
  });

  auto problem = feasible_from_matrix(matrix);

  problem.rhs_bounds = {
      Bound<Rational>{1, 1},
      Bound<Rational>{2, 3},
      Bound<Rational>{-6, -5},
  };

  auto new_problem = presolve::TransformToEqualities<Rational>().apply(problem);
  new_problem.validate();

  ASSERT_EQ(new_problem.matrix.rows(), 3);
  ASSERT_EQ(new_problem.matrix.cols(), 6);

  for (size_t i = 0; i < 3; ++i) {
    ASSERT_TRUE(new_problem.rhs_bounds[i].is_fixed());
  }
}

TEST(TransformToEqualities, GeneratesCorrectNames) {
  auto matrix = sparse<Rational>({
      {1, 2, 3, 4},
      {0, 4, 1, 2},
  });

  auto problem = feasible_from_matrix(matrix);

  problem.row_names = {"r1", "r2"};
  problem.rhs_bounds = {
      Bound<Rational>{1, 1},
      Bound<Rational>{2, 3},
  };

  auto new_problem = presolve::TransformToEqualities<Rational>().apply(problem);
  new_problem.validate();

  ASSERT_EQ(new_problem.matrix.rows(), 2);
  ASSERT_EQ(new_problem.matrix.cols(), 5);

  ASSERT_EQ(new_problem.var_names[4], "r2_range");
}

TEST(TransformToEqualities, DontTouchEqualities) {
  auto matrix = sparse<Rational>({
      {1, 2, 3, 4},
      {0, 4, 1, 2},
  });

  auto problem = feasible_from_matrix(matrix);

  problem.row_names = {"r1", "r2"};
  problem.rhs_bounds = {
      Bound<Rational>{1, 1},
      Bound<Rational>{2, 2},
  };

  auto new_problem = presolve::TransformToEqualities<Rational>().apply(problem);
  new_problem.validate();

  ASSERT_EQ(new_problem.matrix.rows(), 2);
  ASSERT_EQ(new_problem.matrix.cols(), 4);
}

TEST(TransformToEqualities, RandomTests) {
  std::default_random_engine random;

  for (size_t i = 0; i < 1'000; ++i) {
    auto problem = random_feasible_problem<double>(10, 10, random);
    auto pass = presolve::TransformToEqualities<double>();

    ASSERT_PASS_CORRECT(problem, pass);
  }
}
