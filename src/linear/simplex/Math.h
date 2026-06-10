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

}  // namespace simplex::detail
