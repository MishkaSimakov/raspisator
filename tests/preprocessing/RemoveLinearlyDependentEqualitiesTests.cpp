#include <gtest/gtest.h>

#include <random>
#include <set>

#include "obfuscators/ShuffleRows.h"
#include "presolve/passes/RemoveLinearlyDependentEqualities.h"
#include "support/Highs.h"
#include "support/ProblemConstructors.h"
#include "support/RandomProblem.h"

template <typename Field, typename Gen>
  requires std::uniform_random_bit_generator<Gen>
void add_linearly_dependent_constraints(problem::MILP<Field>& problem,
                                        Gen& random) {
  const auto [n, d] = problem.matrix.shape();

  std::uniform_int_distribution<int> elements_distribution(-5, 5);

  const auto multiplier = linalg::random<Field>(
      n, n, [&] { return elements_distribution(random); });

  // add linearly dependent constraints and their bounds
  auto new_matrix = linalg::to_dense(problem.matrix);
  new_matrix = linalg::vstack(new_matrix, multiplier * new_matrix);
  problem.matrix = CSCMatrix(new_matrix);

  problem.rhs_bounds.resize(2 * n, Bound<Field>{0, 0});

  for (size_t row = 0; row < n; ++row) {
    for (size_t col = 0; col < n; ++col) {
      problem.rhs_bounds[n + row] +=
          multiplier[row, col] * problem.rhs_bounds[col];
    }
  }

  // widen rhs bounds, problem should still remain feasible
  std::uniform_int_distribution<int> bound_widening(0, 5);
  std::uniform_int_distribution<int> coin(0, 1);

  for (size_t row = 0; row < 2 * n; ++row) {
    if (coin(random) == 1) {
      continue;
    }

    if (problem.rhs_bounds[row].lower) {
      *problem.rhs_bounds[row].lower -= bound_widening(random);
    }
    if (problem.rhs_bounds[row].upper) {
      *problem.rhs_bounds[row].upper += bound_widening(random);
    }
  }

  problem.row_names.resize(2 * n);
}

TEST(RemoveLinearlyDependentConstraintsTests,
     RemovesLinearlyDependentEqualities) {
  // matrix[2] = matrix[0] - matrix[1]
  CSCMatrix<Rational> matrix = {
      {1, 2, 3, 4},
      {0, 4, 1, 2},
      {1, -2, 2, 2},
  };

  auto problem = feasible_from_matrix(matrix);

  // ensure that all rows are equalities
  problem.rhs_bounds = {
      Bound<Rational>{0, 0},
      Bound<Rational>{0, 0},
      Bound<Rational>{0, 0},
  };

  auto new_problem =
      presolve::RemoveLinearlyDependentEqualities<Rational>().apply(problem);
  new_problem.validate();

  ASSERT_EQ(new_problem.matrix.rows(), 2);
  ASSERT_EQ(new_problem.matrix.cols(), 4);
}

TEST(RemoveLinearlyDependentConstraintsTests, PreservesNames) {
  // matrix[2] = matrix[0] - matrix[1]
  CSCMatrix<Rational> matrix = {
      {1, 2, 3, 4},
      {0, 4, 1, 2},
      {1, -2, 2, 2},
  };

  auto problem = feasible_from_matrix(matrix);

  problem.row_names = {"r0", "r1", "r2"};
  // ensure that all rows are equalities
  problem.rhs_bounds = {
      Bound<Rational>{0, 0},
      Bound<Rational>{0, 0},
      Bound<Rational>{0, 0},
  };

  auto new_problem =
      presolve::RemoveLinearlyDependentEqualities<Rational>().apply(problem);

  std::set remaining_names(new_problem.row_names.begin(),
                           new_problem.row_names.end());

  ASSERT_TRUE((remaining_names == std::set<std::string>{"r0", "r1"} ||
               remaining_names == std::set<std::string>{"r0", "r2"} ||
               remaining_names == std::set<std::string>{"r1", "r2"}));
}

TEST(RemoveLinearlyDependentConstraintsTests, InfeasibilityDetection) {
  CSCMatrix<Rational> matrix = {
      {1, 2, 3},
      {2, 4, 6},
  };

  auto problem = feasible_from_matrix(matrix);

  problem.rhs_bounds = {
      Bound<Rational>{1, 1},
      Bound<Rational>{-3, -3},
  };

  auto new_problem =
      presolve::RemoveLinearlyDependentEqualities<Rational>().apply(problem);

  ASSERT_TRUE(new_problem.proven_infeasible);
}

TEST(RemoveLinearlyDependentConstraintsTests, RandomTests) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 20;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine random;

  for (size_t i = 0; i < kIterations; ++i) {
    auto problem =
        random_feasible_problem<double>(kSize, kElementMagnitude, random);

    add_linearly_dependent_constraints(problem, random);
    shuffle_rows(problem);

    // solve without preprocessing
    auto solution = highs::solve(highs::from_milp(problem));

    // preprocess
    presolve::RemoveLinearlyDependentEqualities<double> pass;

    auto new_problem = pass.apply(problem);

    new_problem.validate();

    // solve new problem
    auto new_solution = highs::solve(highs::from_milp(new_problem));

    ASSERT_EQ(solution.status, new_solution.status);

    if (solution.status == HighsModelStatus::kOptimal) {
      // Objective must match
      ASSERT_NEAR(solution.objective, new_solution.objective, 1e-6);
    }
  }
}
