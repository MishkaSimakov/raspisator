#include <gtest/gtest.h>

#include <random>
#include <set>

#include "obfuscators/ShuffleRows.h"
#include "presolve/passes/RemoveLinearlyDependentConstraints.h"
#include "support/Highs.h"
#include "support/ProblemConstructors.h"
#include "support/RandomProblem.h"

template <typename Field>
void add_linearly_dependent_constraints(problem::MILP<Field>& problem) {
  const auto [n, d] = problem.matrix.shape();

  std::default_random_engine random;
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

  problem.row_names.resize(2 * n);
}

TEST(RemoveLinearlyDependentConstraintsTests, RemovesLinearlyDependent) {
  // matrix[2] = matrix[0] - matrix[1]
  CSCMatrix<Rational> matrix = {
      {1, 2, 3, 4},
      {0, 4, 1, 2},
      {1, -2, 2, 2},
  };

  auto problem = feasible_from_matrix(matrix);

  auto new_problem =
      presolve::RemoveLinearlyDependentConstraints<Rational>().apply(problem);
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

  auto new_problem =
      presolve::RemoveLinearlyDependentConstraints<Rational>().apply(problem);

  std::set remaining_names(new_problem.row_names.begin(),
                           new_problem.row_names.end());

  ASSERT_TRUE((remaining_names == std::set<std::string>{"r0", "r1"} ||
               remaining_names == std::set<std::string>{"r0", "r2"} ||
               remaining_names == std::set<std::string>{"r1", "r2"}));
}

TEST(RemoveLinearlyDependentConstraintsTests, RandomTests) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 10;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine random;

  for (size_t i = 0; i < kIterations; ++i) {
    auto problem =
        random_feasible_problem<double>(kSize, kElementMagnitude, random);

    add_linearly_dependent_constraints(problem);

    shuffle_rows(problem);

    // solve without preprocessing
    auto solution = highs::solve(highs::from_milp(problem));

    // preprocess
    presolve::RemoveLinearlyDependentConstraints<double> pass;

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
