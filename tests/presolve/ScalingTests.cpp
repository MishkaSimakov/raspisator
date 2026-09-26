#include <gtest/gtest.h>

#include <random>
#include <set>

#include "ConstructSparse.h"
#include "PassAssertions.h"
#include "presolve/passes/Scaling.h"
#include "support/Highs.h"
#include "support/ProblemConstructors.h"
#include "support/RandomProblem.h"

TEST(ScalingTests, RandomProblems) {
  std::default_random_engine random;

  for (size_t i = 0; i < 1'000; ++i) {
    auto problem = random_feasible_problem<double>(5, 100, random);
    auto pass = presolve::Scaling<double>();

    ASSERT_PASS_CORRECT(problem, pass);
  }
}
