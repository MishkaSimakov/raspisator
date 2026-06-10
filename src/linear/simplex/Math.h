#pragma once

#include <vector>

#include "linalg/Linalg.h"
#include "linear/model/Bound.h"
#include "linear/model/LP.h"

namespace simplex::detail {

template <typename Field>
Matrix<Field> get_adjusted_rhs(const CSCMatrix<Field>& matrix,
                               const Vector<Field>& rhs,
                               const Bounds<Field>& bounds,
                               const std::vector<VariableState>& states) {
  auto [n, d] = matrix.shape();

  Vector<Field> result(rhs);

  for (size_t col = 0; col < d; ++col) {
    if (states[col] == VariableState::AT_LOWER) {
      for (const auto& [row, value] : matrix.get_column(col)) {
        result[row] -= value * *bounds[col].lower;
      }
    } else if (states[col] == VariableState::AT_UPPER) {
      for (const auto& [row, value] : matrix.get_column(col)) {
        result[row] -= value * *bounds[col].upper;
      }
    }
  }

  return result;
}

template <typename Field>
static Vector<Field> get_point_from_basis(
    const Bounds<Field>& bounds, const std::vector<VariableState>& states,
    const std::vector<size_t>& basic_vars, const Vector<Field>& basic_point) {
  Vector<Field> result(states.size());

  for (size_t i = 0; i < basic_vars.size(); ++i) {
    result[basic_vars[i]] = basic_point[i, 0];
  }
  for (size_t i = 0; i < states.size(); ++i) {
    if (states[i] == VariableState::AT_LOWER) {
      result[i] = *bounds[i].lower;
    } else if (states[i] == VariableState::AT_UPPER) {
      result[i] = *bounds[i].upper;
    } else if (states[i] == VariableState::NONBASIC_FREE) {
      result[i] = 0;
    }
  }

  return result;
}

template <typename Field>
Field get_objective(const Vector<Field>& cost, const Bounds<Field>& bounds,
                    const std::vector<VariableState>& states,
                    const std::vector<size_t>& basic_vars,
                    const Vector<Field>& basic_point) {
  KahanSum<Field> objective;

  for (size_t col = 0; col < states.size(); ++col) {
    // TODO: change this to switch, so that new VariableStates can be handled
    if (states[col] == VariableState::AT_LOWER) {
      objective.add(cost[col] * *bounds[col].lower);
    } else if (states[col] == VariableState::AT_UPPER) {
      objective.add(cost[col] * *bounds[col].upper);
    }
  }

  for (size_t i = 0; i < basic_vars.size(); ++i) {
    objective.add(cost[basic_vars[i]] * basic_point[i]);
  }

  return objective.sum();
}

template <typename Field>
Vector<Field> get_basic_cost(const Vector<Field>& cost,
                             const std::vector<size_t>& basic_vars) {
  Vector<Field> result(basic_vars.size());

  for (size_t i = 0; i < basic_vars.size(); ++i) {
    result[i] = cost[basic_vars[i]];
  }

  return result;
}

// TODO: this can be simplified if I implement sparse matrix arithmetics
// (c - A.transposed() * pi, where .transposed is expression template)
// @simplex_multipliers is pi = A_B^-1 c_B
template <typename Field>
Vector<Field> get_reduced_cost(const CSCMatrix<Field>& A,
                               const Vector<Field>& c,
                               const Vector<Field>& simplex_multipliers) {
  auto [n, d] = A.shape();

  Vector<Field> result(d);
  for (size_t i = 0; i < d; ++i) {
    result[i] = c[i];

    for (const auto& [row, value] : A.get_column(i)) {
      result[i] -= value * simplex_multipliers[row];
    }
  }

  return result;
}

}  // namespace simplex::detail
