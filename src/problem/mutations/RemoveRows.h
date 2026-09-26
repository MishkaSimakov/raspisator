#pragma once

#include <set>

#include "problem/StandardLP.h"

namespace problem {

namespace detail {

inline std::pair<std::vector<size_t>, size_t> get_rows_mapping(
    size_t rows_count, const std::vector<size_t>& rows_to_remove) {
  std::set to_remove_set(rows_to_remove.begin(), rows_to_remove.end());

  std::vector<size_t> rows_mapping(rows_count, rows_count);
  size_t new_rows_count = 0;

  for (size_t row = 0; row < rows_count; ++row) {
    if (!to_remove_set.contains(row)) {
      rows_mapping[row] = new_rows_count++;
    }
  }

  return {rows_mapping, new_rows_count};
}

template <typename Field>
void remove_rows(CoreLP<Field>& problem,
                 const std::vector<size_t>& rows_mapping,
                 size_t new_rows_count) {
  const auto [n, d] = problem.matrix.shape();

  // update constraints matrix
  for (size_t col = 0; col < d; ++col) {
    for (auto& [row, value] : problem.matrix.get_column(col)) {
      row = rows_mapping[row];
    }
  }
  problem.matrix.resize(new_rows_count, d);

  // update row names
  for (size_t row = 0; row < n; ++row) {
    if (rows_mapping[row] < n) {
      problem.row_names[rows_mapping[row]] = problem.row_names[row];
    }
  }
  problem.row_names.resize(new_rows_count);
}

template <typename Field>
void remove_rows(LP<Field>& problem, const std::vector<size_t>& rows_mapping,
                 size_t new_rows_count) {
  const auto [n, d] = problem.matrix.shape();

  remove_rows(static_cast<CoreLP<Field>&>(problem), rows_mapping,
              new_rows_count);

  // update rhs
  for (size_t row = 0; row < n; ++row) {
    if (rows_mapping[row] < n) {
      problem.rhs_bounds[rows_mapping[row]] = problem.rhs_bounds[row];
    }
  }
  problem.rhs_bounds.resize(new_rows_count);
}

template <typename Field>
void remove_rows(StandardLP<Field>& problem,
                 const std::vector<size_t>& rows_mapping,
                 size_t new_rows_count) {
  const auto [n, d] = problem.matrix.shape();

  remove_rows(static_cast<CoreLP<Field>&>(problem), rows_mapping,
              new_rows_count);

  // update rhs
  for (size_t row = 0; row < n; ++row) {
    if (rows_mapping[row] < n) {
      problem.rhs[rows_mapping[row]] = problem.rhs[row];
    }
  }
  problem.rhs.resize(new_rows_count);
}

}  // namespace detail

template <typename Field>
LP<Field> remove_rows(LP<Field> problem,
                      const std::vector<size_t>& rows_to_remove) {
  auto [rows_mapping, new_rows_count] =
      detail::get_rows_mapping(problem.matrix.rows(), rows_to_remove);

  detail::remove_rows(problem, rows_mapping, new_rows_count);

  return std::move(problem);
}

template <typename Field>
MILP<Field> remove_rows(MILP<Field> problem,
                        const std::vector<size_t>& rows_to_remove) {
  auto [rows_mapping, new_rows_count] =
      detail::get_rows_mapping(problem.matrix.rows(), rows_to_remove);

  detail::remove_rows(static_cast<LP<Field>&>(problem), rows_mapping,
                      new_rows_count);

  return std::move(problem);
}

template <typename Field>
StandardLP<Field> remove_rows(StandardLP<Field> problem,
                              const std::vector<size_t>& rows_to_remove) {
  auto [rows_mapping, new_rows_count] =
      detail::get_rows_mapping(problem.matrix.rows(), rows_to_remove);

  detail::remove_rows(problem, rows_mapping, new_rows_count);

  return std::move(problem);
}

template <typename Field>
StandardMILP<Field> remove_rows(StandardMILP<Field> problem,
                                const std::vector<size_t>& rows_to_remove) {
  auto [rows_mapping, new_rows_count] =
      detail::get_rows_mapping(problem.matrix.rows(), rows_to_remove);

  detail::remove_rows(static_cast<StandardLP<Field>&>(problem), rows_mapping,
                      new_rows_count);

  return std::move(problem);
}

}  // namespace problem
