#include <gtest/gtest.h>

#include "faker/Catalog.h"
#include "field/BigInteger.h"
#include "presolve/obfuscators/AddLinearlyDependentRows.h"
#include "presolve/obfuscators/ShuffleRows.h"
#include "presolve/passes/TransformToEqualities.h"
#include "simplex/init/dual/Subproblem.h"
#include "support/RandomProblem.h"

TEST(DualSubproblemPhase1Tests, CatalogProblems) {
  const auto problems = faker::catalog<Rational>()
                            .problem_type(faker::ProblemType::StandardLP)
                            .solution_type(faker::SolutionType::FEASIBLE)
                            .linearly_dependent_rows(std::nullopt)
                            .all();

  for (auto instance : problems) {
    problem::StandardLP problem(instance.problem);

    auto result = simplex::subproblem_dual_phase1(problem);

    ASSERT_TRUE(result.has_value());
    ASSERT_EQ(simplex::get_dual_infeasibility_reason(problem, result->states),
              "");
  }
}

TEST(DualSubproblemPhase1Tests, RandomProblems) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 5;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine engine(0);

  for (size_t iteration = 0; iteration < kIterations; ++iteration) {
    std::cout << "#" << iteration << std::endl;

    auto problem =
        random_feasible_problem<Rational>(kSize, kElementMagnitude, engine);

    add_linearly_dependent_rows(problem, engine);
    shuffle_rows(problem, engine);

    problem = presolve::TransformToEqualities<Rational>().apply(problem);

    problem::StandardLP<Rational> standard_lp(problem);

    auto phase1 = simplex::subproblem_dual_phase1(standard_lp);

    ASSERT_TRUE(phase1.has_value());
    ASSERT_TRUE(simplex::is_dual_feasible(standard_lp, phase1->states));
  }
}
