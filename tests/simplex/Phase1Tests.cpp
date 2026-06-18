#include <gtest/gtest.h>

#include "faker/Catalog.h"
#include "field/BigInteger.h"
#include "linalg/RRQR.h"
#include "linalg/Rank.h"
#include "problem/mutations/RemoveRows.h"
#include "simplex/Feasibility.h"
#include "simplex/init/primal/Phase1.h"
#include "support/Highs.h"

TEST(PrimalPhase1Tests, CatalogProblems) {
  const auto problems = faker::catalog<Rational>()
                            .problem_type(faker::ProblemType::StandardLP)
                            .solution_type(faker::SolutionType::FEASIBLE)
                            .linearly_dependent_rows(false)
                            .all();

  for (auto instance : problems) {
    problem::StandardLP problem(instance.problem);

    auto result = simplex::primal_phase1(problem);

    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(simplex::is_primal_feasible(problem, result->states));
    ASSERT_TRUE(result->redundant_rows.empty());
  }
}

TEST(PrimalPhase1Tests, CatalogInfeasibleProblems) {
  const auto problems = faker::catalog<Rational>()
                            .problem_type(faker::ProblemType::StandardLP)
                            .linearly_dependent_rows(false)
                            .solution_type(faker::SolutionType::INFEASIBLE)
                            .all();

  for (auto instance : problems) {
    problem::StandardLP problem(instance.problem);

    auto states = simplex::primal_phase1(problem);

    ASSERT_FALSE(states.has_value());
    ASSERT_EQ(states.error(), simplex::Phase1Error::INFEASIBLE);
  }
}

TEST(PrimalPhase1Tests, CatalogProblemsWithLinearlyDependentRows) {
  const auto problems = faker::catalog<double>()
                            .problem_type(faker::ProblemType::StandardLP)
                            .solution_type(faker::SolutionType::FEASIBLE)
                            .linearly_dependent_rows(true)
                            .all();

  for (auto instance : problems) {
    problem::StandardLP problem(instance.problem);

    auto result = simplex::primal_phase1(problem);

    ASSERT_TRUE(result.has_value());

    auto old_solution = highs::solve(highs::from_milp(problem::MILP(problem)));

    auto new_problem = problem::remove_rows(problem, result->redundant_rows);

    ASSERT_EQ(linalg::rank(Matrix(new_problem.matrix)), new_problem.matrix.rows());

    auto new_solution =
        highs::solve(highs::from_milp(problem::MILP(new_problem)));

    ASSERT_EQ(old_solution.status, HighsModelStatus::kOptimal);
    ASSERT_EQ(new_solution.status, HighsModelStatus::kOptimal);

    ASSERT_DOUBLE_EQ(old_solution.objective, new_solution.objective);
  }
}
