#include <gtest/gtest.h>

#include "faker/Catalog.h"
#include "field/BigInteger.h"
#include "presolve/obfuscators/AddLinearlyDependentRows.h"
#include "presolve/obfuscators/ShuffleRows.h"
#include "presolve/passes/TransformToEqualities.h"
#include "simplex/init/dual/Subproblem.h"
#include "support/RandomProblem.h"

TEST(ReducedCostTests, CatalogProblems) {
  const auto problems = faker::catalog<Rational>()
                            .problem_type(faker::ProblemType::StandardLP)
                            .solution_type(faker::SolutionType::FEASIBLE)
                            .all_variables_bounded(true)
                            .all();

  for (auto instance : problems) {
    problem::StandardLP problem(instance.problem);

    auto result = simplex::try_init_dual_by_reduced_cost(problem);

    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(simplex::is_dual_feasible(problem, *result));
  }
}
