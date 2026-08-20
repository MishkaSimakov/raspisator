#include <gtest/gtest.h>

#include "faker/Catalog.h"
#include "field/BigInteger.h"
#include "presolve/obfuscators/AddInfiniteBounds.h"
#include "presolve/obfuscators/AddLinearlyDependentRows.h"
#include "presolve/obfuscators/ShuffleRows.h"
#include "presolve/passes/TransformToEqualities.h"
#include "simplex/init/dual/Subproblem.h"
#include "simplex2/Assertions.h"
#include "support/Highs.h"
#include "support/RandomProblem.h"

TEST(DualSubproblemPhase1Tests, CatalogProblems) {
  const auto problems = faker::catalog<Rational>()
                            .problem_type(faker::ProblemType::StandardLP)
                            .solution_type(faker::SolutionType::FEASIBLE)
                            .linearly_dependent_rows(false)
                            .all();

  for (auto instance : problems) {
    if (instance.problem.name != "textbook13") {
      continue;
    }

    problem::StandardLP problem(instance.problem);

    auto result = simplex::subproblem_dual_phase1(problem);

    ASSERT_TRUE(result.has_value());
    ASSERT_DUAL_FEASIBLE(problem, result->states);
  }
}

TEST(DualSubproblemPhase1Tests, RandomBoundedProblems) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 5;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine engine(0);

  for (size_t iteration = 0; iteration < kIterations; ++iteration) {
    std::cout << "#" << iteration << std::endl;

    auto problem =
        random_feasible_problem<Rational>(kSize, kElementMagnitude, engine);

    add_infinite_bounds(problem, engine, true);

    problem = presolve::TransformToEqualities<Rational>().apply(problem);

    problem::StandardLP<Rational> standard_lp(problem);

    auto phase1 = simplex::subproblem_dual_phase1(standard_lp);

    // dual phase 1 must not fail since problem is dual feasible (because it is
    // bounded and primal feasible)
    ASSERT_TRUE(phase1.has_value());
    ASSERT_DUAL_FEASIBLE(standard_lp, phase1->states);
  }
}

TEST(DualSubproblemPhase1Tests, RandomUnboundedProblems) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 5;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine engine(0);

  for (size_t iteration = 0; iteration < kIterations; ++iteration) {
    std::cout << "#" << iteration << std::endl;

    auto problem =
        random_feasible_problem<Rational>(kSize, kElementMagnitude, engine);

    add_infinite_bounds(problem, engine, false);

    problem = presolve::TransformToEqualities<Rational>().apply(problem);

    problem::StandardLP<Rational> standard_lp(problem);

    auto phase1 = simplex::subproblem_dual_phase1(standard_lp);
    auto solution = highs::solve(highs::from_milp(problem::MILP(standard_lp)));

    if (!phase1.has_value()) {
      ASSERT_EQ(phase1.error(),
                simplex::SubproblemPhase1Error::DUAL_INFEASIBLE);
      ASSERT_TRUE(solution.status == HighsModelStatus::kUnboundedOrInfeasible ||
                  solution.status == HighsModelStatus::kUnbounded ||
                  solution.status == HighsModelStatus::kInfeasible);
    } else {
      ASSERT_DUAL_FEASIBLE(standard_lp, phase1->states);
    }
  }
}
