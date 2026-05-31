#pragma once

#include <optional>
#include <vector>

#include "linear/matrix/RowBasis.h"
#include "linear/model/LP.h"
#include "linear/simplex/Math.h"
#include "linear/sparse/CSCMatrix.h"
#include "linear/sparse/LU.h"

namespace simplex {

// Algorithm for finding initial dual feasible point.
// It is fast, but may fail. It is guaranteed to work when all variables have
// both upper and lower bounds.
template <typename Field>
std::optional<std::vector<VariableState>> try_init_dual_by_reduced_cost(
    const CSCMatrix<Field>& A, const Matrix<Field>& b, const Matrix<Field>& c,
    const Bounds<Field>& bounds, const std::vector<size_t>& basic_variables) {
  auto [n, d] = A.shape();

  if (basic_variables.size() != n) {
    throw std::invalid_argument(std::format(
        "Wrong basic variables count: {} != {}", basic_variables.size(), n));
  }

  std::vector<VariableState> states(d);

  linalg::LUPA<Field> lupa(A);
  lupa.set_columns(basic_variables);

  const auto simplex_multipliers =
      lupa.solve_linear_transposed(detail::get_basic_cost(c, basic_variables));

  const auto reduced_costs =
      detail::get_reduced_cost(A, c, simplex_multipliers);

  for (size_t i = 0; i < d; ++i) {
    if (reduced_costs[i, 0] <= Field(0)) {
      if (!bounds[i].lower) {
        return std::nullopt;
      }

      states[i] = VariableState::AT_LOWER;
    } else {
      if (!bounds[i].upper) {
        return std::nullopt;
      }

      states[i] = VariableState::AT_UPPER;
    }
  }

  for (const size_t basic_var : basic_variables) {
    states[basic_var] = VariableState::BASIC;
  }

  return states;
}

template <typename Field>
std::optional<std::vector<VariableState>> try_init_dual_by_reduced_cost(
    const CSCMatrix<Field>& A, const Matrix<Field>& b, const Matrix<Field>& c,
    const Bounds<Field>& bounds) {
  auto [n, d] = A.shape();

  auto basic_variables =
      linalg::get_row_basis(linalg::transposed(linalg::to_dense(A)));

  return try_init_dual_by_reduced_cost(A, b, c, bounds, basic_variables);
}

}  // namespace simplex
