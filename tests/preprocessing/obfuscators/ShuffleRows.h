#pragma once

#include <algorithm>
#include <numeric>
#include <random>

#include "problem/MILP.h"

template <typename Field>
void shuffle_rows(problem::MILP<Field>& problem, size_t seed = 0) {
  const auto [n, d] = problem.matrix.shape();

  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), 0);

  std::default_random_engine random(seed);
  std::ranges::shuffle(order, random);

  std::vector<Bound<Field>> new_rhs_bounds(n);
  std::vector<std::string> new_row_names(n);

  for (size_t row = 0; row < n; ++row) {
    new_rhs_bounds[order[row]] = problem.rhs_bounds[row];
    new_row_names[order[row]] = problem.row_names[row];
  }

  problem.rhs_bounds = std::move(new_rhs_bounds);
  problem.row_names = std::move(new_row_names);

  // TODO: implement permutation operator for CSCMatrix
  for (size_t col = 0; col < d; ++col) {
    for (auto& [row, value] : problem.matrix.get_column(col)) {
      row = order[row];
    }
  }
}
