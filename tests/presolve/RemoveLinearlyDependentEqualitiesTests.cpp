#include <gtest/gtest.h>

#include <random>
#include <set>

#include "ConstructSparse.h"
#include "PassAssertions.h"
#include "obfuscators/ShuffleRows.h"
#include "obfuscators/AddLinearlyDependentRows.h"
#include "presolve/passes/RemoveLinearlyDependentEqualities.h"
#include "support/Highs.h"
#include "support/ProblemConstructors.h"
#include "support/RandomProblem.h"

TEST(RemoveLinearlyDependentConstraintsTests,
     RemovesLinearlyDependentEqualities) {
  // matrix[2] = matrix[0] - matrix[1]
  auto matrix = sparse<double>({
      {1, 2, 3, 4},
      {0, 4, 1, 2},
      {1, -2, 2, 2},
  });

  auto problem = feasible_from_matrix(matrix);

  // ensure that all rows are equalities
  problem.rhs_bounds = {
      Bound<double>{0, 0},
      Bound<double>{0, 0},
      Bound<double>{0, 0},
  };

  auto new_problem =
      presolve::RemoveLinearlyDependentEqualities<double>().apply(problem);
  new_problem.validate();

  ASSERT_EQ(new_problem.matrix.rows(), 2);
  ASSERT_EQ(new_problem.matrix.cols(), 4);
}

TEST(RemoveLinearlyDependentConstraintsTests, SmallTest) {
  auto matrix = sparse<double>({
      {1, 2, 3},
      {2, 4, 6},
  });

  auto problem = feasible_from_matrix(matrix);

  problem.rhs_bounds = {
      Bound<double>{1, 1},
      Bound<double>{2, 2},
  };

  auto new_problem =
      presolve::RemoveLinearlyDependentEqualities<double>().apply(problem);

  ASSERT_EQ(new_problem.matrix.rows(), 1);
  ASSERT_EQ(new_problem.matrix.cols(), 3);

  ASSERT_FALSE(new_problem.proven_infeasible);
}

TEST(RemoveLinearlyDependentConstraintsTests, SmallTest2) {
  auto matrix = sparse<double>({
      {1, 2, 3},
      {1, 2, 3},
      {2, 4, 6},
  });

  auto problem = feasible_from_matrix(matrix);

  problem.rhs_bounds = {
      Bound<double>{1, 1},
      Bound<double>{1, 1},
      Bound<double>{2, 2},
  };

  auto new_problem =
      presolve::RemoveLinearlyDependentEqualities<double>().apply(problem);

  ASSERT_EQ(new_problem.matrix.rows(), 1);
  ASSERT_EQ(new_problem.matrix.cols(), 3);

  ASSERT_FALSE(new_problem.proven_infeasible);
}

TEST(RemoveLinearlyDependentConstraintsTests, PreservesNames) {
  // matrix[2] = matrix[0] - matrix[1]
  auto matrix = sparse<double>({
      {1, 2, 3, 4},
      {0, 4, 1, 2},
      {1, -2, 2, 2},
  });

  auto problem = feasible_from_matrix(matrix);

  problem.row_names = {"r0", "r1", "r2"};
  // ensure that all rows are equalities
  problem.rhs_bounds = {
      Bound<double>{0, 0},
      Bound<double>{0, 0},
      Bound<double>{0, 0},
  };

  auto new_problem =
      presolve::RemoveLinearlyDependentEqualities<double>().apply(problem);

  std::set remaining_names(new_problem.row_names.begin(),
                           new_problem.row_names.end());

  ASSERT_TRUE((remaining_names == std::set<std::string>{"r0", "r1"} ||
               remaining_names == std::set<std::string>{"r0", "r2"} ||
               remaining_names == std::set<std::string>{"r1", "r2"}));
}

TEST(RemoveLinearlyDependentConstraintsTests, InfeasibilityDetection1) {
  auto matrix = sparse<double>({
      {1, 2, 3},
      {2, 4, 6},
  });

  auto problem = feasible_from_matrix(matrix);

  problem.rhs_bounds = {
      Bound<double>{1, 1},
      Bound<double>{-3, -3},
  };

  auto new_problem =
      presolve::RemoveLinearlyDependentEqualities<double>().apply(problem);

  ASSERT_TRUE(new_problem.proven_infeasible);
}

TEST(RemoveLinearlyDependentConstraintsTests, RandomTests) {
  std::default_random_engine random;

  for (size_t i = 0; i < 1'000; ++i) {
    auto problem = random_feasible_problem<double>(10, 10, random);

    add_linearly_dependent_rows(problem, random);
    shuffle_rows(problem, random);

    presolve::RemoveLinearlyDependentEqualities<double> pass;

    ASSERT_PASS_CORRECT(problem, pass);
  }
}
