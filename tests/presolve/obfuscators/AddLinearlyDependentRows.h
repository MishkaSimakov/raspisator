#pragma once

#include "linalg/Linalg.h"
#include "linalg/Random.h"
#include "linalg/Stack.h"

template <typename Field, typename Gen>
  requires std::uniform_random_bit_generator<Gen>
void add_linearly_dependent_rows(problem::MILP<Field>& problem, Gen& random) {
  const auto [n, d] = problem.matrix.shape();

  std::uniform_int_distribution<int> value_distribution(-5, 5);

  const auto multiplier =
      linalg::random::dense<Field>(n, n, random, value_distribution);

  // add linearly dependent constraints and their bounds
  auto new_matrix = Matrix(problem.matrix);
  new_matrix = linalg::vstack(new_matrix, multiplier * new_matrix);
  problem.matrix = CSCMatrix(new_matrix);

  problem.rhs_bounds.resize(2 * n, Bound<Field>{0, 0});

  for (size_t row = 0; row < n; ++row) {
    for (size_t col = 0; col < n; ++col) {
      problem.rhs_bounds[n + row] +=
          multiplier[row, col] * problem.rhs_bounds[col];
    }
  }

  // widen rhs bounds, problem should still remain feasible
  std::uniform_int_distribution<int> bound_widening(0, 5);
  std::uniform_int_distribution<int> coin(0, 1);

  for (size_t row = 0; row < 2 * n; ++row) {
    if (coin(random) == 1) {
      continue;
    }

    if (problem.rhs_bounds[row].lower) {
      *problem.rhs_bounds[row].lower -= bound_widening(random);
    }
    if (problem.rhs_bounds[row].upper) {
      *problem.rhs_bounds[row].upper += bound_widening(random);
    }
  }

  problem.row_names.resize(2 * n);
}
