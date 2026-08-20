#pragma once

#include "linalg/Linalg.h"
#include "linalg/Random.h"
#include "linalg/Stack.h"
#include "utils/Random.h"

// Adds infinite bounds to the problem's variables.
// If preserve_boundness is true, bounded problems will remain bounded.
// Otherwise, a problem can become unbounded after this method.
template <typename Field, typename Gen>
  requires std::uniform_random_bit_generator<Gen>
void add_infinite_bounds(problem::MILP<Field>& problem, Gen& random,
                         bool preserve_boundness) {
  const auto [n, d] = problem.matrix.shape();

  for (size_t i = 0; i < d; ++i) {
    if (rnd::bernoulli(0.25, random) && (!preserve_boundness || problem.cost[i] >= 0)) {
      problem.var_bounds[i].lower = std::nullopt;
    }

    if (rnd::bernoulli(0.25, random) && (!preserve_boundness || problem.cost[i] <= 0)) {
      problem.var_bounds[i].upper = std::nullopt;
    }
  }

  problem.implied_var_bounds = problem.var_bounds;
}
