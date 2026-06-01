#pragma once

#include <variant>

#include "linear/simplex/Types.h"

namespace simplex::detail {

struct ChangeBasis {
  size_t entering_variable;
  size_t leaving_index;
  VariableState new_state;
};

struct ToggleBound {
  size_t variable;
  VariableState new_state;  // should be either AT_UPPER or AT_LOWER
};

struct Unbounded {};

template <typename Field>
std::variant<ChangeBasis, ToggleBound, Unbounded> primal_ratio_test(
    const CSCMatrix<Field>& A, size_t entering_var, Field entering_reduced_cost,
    const IterationState<Field>& state, Field pivot_tolerance) {
  using std::abs;

  const auto [n, d] = A.shape();

  Matrix<Field> column(n, 1, 0);
  for (auto [row, value] : A.get_column(entering_var)) {
    column[row, 0] = value;
  }

  column = state.lupa.solve_linear(std::move(column));

  // theta is maximum entering variable change so that the point would not
  // become primal infeasible for variable i,
  // change = | new_value - old_value |
  const auto get_variable_theta =
      [&](const size_t i, const Field epsilon = 0) -> std::optional<Field> {
    if (abs(column[i, 0]) < pivot_tolerance) {
      return std::nullopt;
    }

    // x_i = \alpha + s \beta x_j, where
    // x_j is entering variable,
    // s = -sign(entering reduced cost)
    const Field alpha = state.basic_point[i, 0];
    const Field beta = entering_reduced_cost > 0 ? -column[i, 0] : column[i, 0];

    const auto bound = (*state.bounds)[state.basic_variables[i]];

    if (beta > 0 && bound.upper) {
      if (alpha > *bound.upper) {
        return epsilon / beta;
      }

      return (*bound.upper - alpha + epsilon) / beta;
    }
    if (beta < 0 && bound.lower) {
      if (alpha < *bound.lower) {
        return -epsilon / beta;
      }

      return (*bound.lower - alpha - epsilon) / beta;
    }

    return std::nullopt;
  };

  // Harris' ratio test
  // There are 2 steps:
  // 1. Determine theta_max
  // 2. Filter variables using theta_max and determine theta_chosen
  Minimum<Field> min_theta_bound;

  for (size_t i = 0; i < n; ++i) {
    min_theta_bound.record(
        get_variable_theta(i, FieldTraits<Field>::tolerance));
  }

  std::optional<size_t> leaving_id = std::nullopt;

  const auto entering_state = state.variables_states[entering_var];
  const auto entering_bound = (*state.bounds)[entering_var];

  if (min_theta_bound.has_value()) {
    const Field theta_max = *min_theta_bound;
    // logging::log_value(*theta_max.min(), "theta_max.txt");

    ArgMaximum<Field> max_pivot;

    for (size_t i = 0; i < n; ++i) {
      auto current_theta = get_variable_theta(i);

      if (current_theta && *current_theta <= theta_max) {
        max_pivot.record(i, abs(column[i, 0]));
      }
    }

    // logging::log_value(*max_pivot.max(), "max_pivot.txt");

    const Field leaving_theta = *get_variable_theta(max_pivot->index);

    // logging::log_value(leaving_theta, "leaving_theta.txt");
    // logging::log_value(state.basic_point[*max_pivot.argmax(), 0],
    // "leaving_value.txt");

    Field new_entering_value;

    switch (entering_state) {
      case VariableState::AT_LOWER:
        new_entering_value = *entering_bound.lower + leaving_theta;
        break;
      case VariableState::AT_UPPER:
        new_entering_value = *entering_bound.upper - leaving_theta;
        break;
      case VariableState::NONBASIC_FREE:
        new_entering_value =
            entering_reduced_cost > 0 ? leaving_theta : -leaving_theta;
        break;
      default:
        throw std::runtime_error("Unexpected variable state.");
    }

    if (entering_bound.contains(new_entering_value)) {
      leaving_id = max_pivot->index;
    }
  }

  if (!leaving_id) {
    if (entering_state == VariableState::AT_LOWER && entering_bound.upper) {
      return ToggleBound{entering_var, VariableState::AT_UPPER};
    }
    if (entering_state == VariableState::AT_UPPER && entering_bound.lower) {
      return ToggleBound{entering_var, VariableState::AT_LOWER};
    }

    return Unbounded{};
  }

  const auto leaving_state = entering_reduced_cost * column[*leaving_id, 0] < 0
                                 ? VariableState::AT_UPPER
                                 : VariableState::AT_LOWER;

  return ChangeBasis{entering_var, *leaving_id, leaving_state};
}

}  // namespace simplex::detail
