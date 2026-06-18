#pragma once

#include <set>

#include "problem/StandardLP.h"

namespace problem {

template <typename Field>
StandardLP<Field> remove_rows(StandardLP<Field> problem,
                              const std::vector<size_t>& rows_to_remove) {
  const auto [n, d] = problem.matrix.shape();

  std::set to_remove_set(rows_to_remove.begin(), rows_to_remove.end());

  std::vector<size_t> rows_mapping(n, n);
  size_t new_rows_count = 0;

  for (size_t row = 0; row < n; ++row) {
    if (!to_remove_set.contains(row)) {
      rows_mapping[row] = new_rows_count++;
    }
  }

  // update constraints matrix
  for (size_t col = 0; col < d; ++col) {
    for (auto& [row, value] : problem.matrix.get_column(col)) {
      row = rows_mapping[row];
    }
  }
  problem.matrix.resize(new_rows_count, d);

  // update rhs
  for (size_t row = 0; row < n; ++row) {
    if (rows_mapping[row] < n) {
      problem.rhs[rows_mapping[row]] = problem.rhs[row];
    }
  }
  problem.rhs.resize(new_rows_count);

  // update row names
  for (size_t row = 0; row < n; ++row) {
    if (rows_mapping[row] < n) {
      problem.row_names[rows_mapping[row]] = problem.row_names[row];
    }
  }
  problem.row_names.resize(new_rows_count);

  return std::move(problem);
}

}  // namespace problem
