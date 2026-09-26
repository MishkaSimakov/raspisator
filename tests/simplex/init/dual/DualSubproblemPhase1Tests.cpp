#include <gtest/gtest.h>

#include "faker/Catalog.h"
#include "presolve/obfuscators/AddInfiniteBounds.h"
#include "presolve/obfuscators/ShuffleRows.h"
#include "presolve/passes/TransformToEqualities.h"
#include "simplex/init/dual/Subproblem.h"
#include "simplex2/Assertions.h"
#include "support/GMPRational.h"
#include "support/Highs.h"
#include "support/RandomProblem.h"

TEST(DualSubproblemPhase1Tests, CatalogProblems) {
  const auto problems = faker::catalog<GMPRational>()
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

TEST(DualSubproblemPhase1Tests, DISABLED_RandomBoundedProblems) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 5;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine engine(0);

  for (size_t iteration = 0; iteration < kIterations; ++iteration) {
    std::cout << "#" << iteration << std::endl;

    auto problem =
        random_feasible_problem<GMPRational>(kSize, kElementMagnitude, engine);

    add_infinite_bounds(problem, engine, true);

    problem = presolve::TransformToEqualities<GMPRational>().apply(problem);

    problem::StandardLP<GMPRational> standard_lp(problem);

    auto phase1 = simplex::subproblem_dual_phase1(standard_lp);

    // dual phase 1 must not fail since problem is dual feasible (because it is
    // bounded and primal feasible)
    ASSERT_TRUE(phase1.has_value());
    ASSERT_DUAL_FEASIBLE(standard_lp, phase1->states);
  }
}

TEST(DualSubproblemPhase1Tests, DISABLED_RandomUnboundedProblems) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 5;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine engine(0);

  for (size_t iteration = 0; iteration < kIterations; ++iteration) {
    std::cout << "#" << iteration << std::endl;

    auto problem =
        random_feasible_problem<GMPRational>(kSize, kElementMagnitude, engine);

    add_infinite_bounds(problem, engine, false);

    problem = presolve::TransformToEqualities<GMPRational>().apply(problem);

    problem::StandardLP<GMPRational> standard_lp(problem);

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
