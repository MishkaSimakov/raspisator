#include <gtest/gtest.h>

#include <random>

#include "support/RandomProblem.h"

auto add_linearly_dependent() {

}

TEST(RemoveLinearlyDependentConstraintsTests, RandomTests) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 10;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine random;

  for (size_t i = 0; i < kIterations; ++i) {
    auto [A, b, c, bounds] = random_feasible_problem(kSize, kElementMagnitude, random);


  }
}
