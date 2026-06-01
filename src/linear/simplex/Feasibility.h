#pragma once

#include "linear/model/LP.h"
#include "linear/simplex/Math.h"
#include "linear/sparse/CSCMatrix.h"
#include "linear/sparse/LU.h"

namespace simplex {

// This function does not check whether matrix formed by basic columns is
// invertible.
template <typename Field>
bool is_dual_feasible(const CSCMatrix<Field>& A, const Matrix<Field>& b,
                      const Matrix<Field>& c, const Bounds<Field>& bounds,
                      const std::vector<VariableState>& states) {
  auto [n, d] = A.shape();

  std::vector<size_t> basic_variables;
  for (size_t i = 0; i < states.size(); ++i) {
    if (states[i] == VariableState::BASIC) {
      basic_variables.push_back(i);
    }
  }

  if (basic_variables.size() != n) {
    return false;
  }

  linalg::LUPA<Field> lupa(A);
  lupa.set_columns(basic_variables);

  const auto simplex_multipliers =
      lupa.solve_linear_transposed(detail::get_basic_cost(c, basic_variables));

  const auto reduced_costs =
      detail::get_reduced_cost(A, c, simplex_multipliers);

  for (size_t i = 0; i < states.size(); ++i) {
    if ((states[i] == VariableState::AT_LOWER &&
         (!bounds[i].lower ||
          FieldTraits<Field>::is_strictly_positive(reduced_costs[i, 0]))) ||
        (states[i] == VariableState::AT_UPPER &&
         (!bounds[i].upper ||
          FieldTraits<Field>::is_strictly_negative(reduced_costs[i, 0])))) {
      return false;
    }
  }

  return true;
}

// This function does not check whether matrix formed by basic columns is
// invertible.
template <typename Field>
bool is_primal_feasible(const CSCMatrix<Field>& A, const Matrix<Field>& b,
                        const Matrix<Field>& c, const Bounds<Field>& bounds,
                        const std::vector<VariableState>& states) {
  auto [n, d] = A.shape();

  std::vector<size_t> basic_variables;
  for (size_t i = 0; i < states.size(); ++i) {
    if (states[i] == VariableState::BASIC) {
      basic_variables.push_back(i);
    }
  }

  if (basic_variables.size() != n) {
    return false;
  }

  auto rhs = detail::get_adjusted_rhs(A, b, bounds, states);

  linalg::LUPA<Field> lupa(A);
  lupa.set_columns(basic_variables);
  auto basic_point = lupa.solve_linear(rhs);

  for (size_t i = 0; i < n; ++i) {
    if (!bounds[basic_variables[i]].contains(basic_point[i, 0])) {
      return false;
    }
  }

  return true;
}

}  // namespace simplex
