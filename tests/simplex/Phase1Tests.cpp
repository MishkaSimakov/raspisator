#include <gtest/gtest.h>

#include "faker/Catalog.h"
#include "linear/BigInteger.h"
#include "linear/simplex/Feasibility.h"
#include "linear/simplex/init/primal/Phase1.h"

TEST(PrimalPhase1Tests, CatalogProblems) {
  const auto problems = faker::catalog<Rational>()
                            .problem_type(faker::ProblemType::StandardLP)
                            .solution_type(faker::SolutionType::FEASIBLE)
                            .linearly_dependent_rows(false)
                            .all();

  for (auto instance : problems) {
    problem::StandardLP problem(instance.problem);

    auto states = simplex::primal_phase1(problem);

    ASSERT_TRUE(states.has_value());
    ASSERT_TRUE(simplex::is_primal_feasible(problem, *states));
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
  const auto problems = faker::catalog<Rational>()
                            .problem_type(faker::ProblemType::StandardLP)
                            .solution_type(faker::SolutionType::FEASIBLE)
                            .linearly_dependent_rows(true)
                            .all();

  for (auto instance : problems) {
    problem::StandardLP problem(instance.problem);

    auto states = simplex::primal_phase1(problem);

    ASSERT_FALSE(states.has_value());
    ASSERT_EQ(states.error(), simplex::Phase1Error::LINEARLY_DEPENDENT_ROWS);
  }
}
