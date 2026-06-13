#pragma once

#include <optional>
#include <vector>

#include "linalg/CSCMatrix.h"
#include "linalg/RowBasis.h"
#include "linalg/Transpose.h"
#include "linalg/lu/LUPA.h"
#include "simplex/SimplexMath.h"

namespace simplex {

// Algorithm for finding initial dual feasible point.
// It is fast, but may fail. It is guaranteed to work when all variables have
// both upper and lower bounds.
template <typename Field>
std::optional<std::vector<VariableState>> try_init_dual_by_reduced_cost(
    const CSCMatrix<Field>& A, const Vector<Field>& b, const Vector<Field>& c,
    const std::vector<Bound<Field>>& bounds,
    const std::vector<size_t>& basic_variables) {
  auto [n, d] = A.shape();

  if (basic_variables.size() != n) {
    throw std::invalid_argument(std::format(
        "Wrong basic variables count: {} != {}", basic_variables.size(), n));
  }

  std::vector<VariableState> states(d);

  linalg::LUPA<Field> lupa(A);
  lupa.set_columns(basic_variables);

  const Vector pi = lupa.solve_linear_transposed(Vector(c[basic_variables]));
  const Vector reduced_costs = c - linalg::transpose(A) * pi;

  for (size_t i = 0; i < d; ++i) {
    if (reduced_costs[i] <= Field(0)) {
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
    const CSCMatrix<Field>& A, const Vector<Field>& b, const Vector<Field>& c,
    const std::vector<Bound<Field>>& bounds) {
  auto [n, d] = A.shape();

  auto basic_variables = linalg::get_row_basis(Matrix(linalg::transpose(A)));

  return try_init_dual_by_reduced_cost(A, b, c, bounds, basic_variables);
}

template <typename Field>
std::optional<std::vector<VariableState>> try_init_dual_by_reduced_cost(
    const problem::StandardLP<Field>& problem) {
  return try_init_dual_by_reduced_cost(problem.matrix, problem.rhs,
                                       problem.cost, problem.var_bounds);
}

}  // namespace simplex
