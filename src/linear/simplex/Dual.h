#pragma once

#include <optional>
#include <random>

#include "Types.h"

namespace simplex {

template <typename Field>
class DualLeavingVariable {
  constexpr static size_t kAnticyclingIterations = 5;

  std::mt19937 random_engine_;

  // largest violation variable selection
  // turned out to be prone to cycling when coupled with strong branching
  std::optional<LeavingVariable> get_dual_leaving_variable_dantzig(
      const IterationState<Field>& state) {
    ArgMaximum<BoundViolation<Field>, std::less<BoundViolation<Field>>>
        max_violation;

    for (size_t i = 0; i < state.basic_variables.size(); ++i) {
      auto violation = (*state.bounds)[state.basic_variables[i]].get_violation(
          state.basic_point[i, 0]);

      max_violation.record(i, violation);
    }

    if (!max_violation.max()) {
      return std::nullopt;
    }

    switch (max_violation.max()->type) {
      case BoundViolationType::NONE:
        return std::nullopt;
      case BoundViolationType::VIOLATE_LOWER_BOUND:
        return LeavingVariable{*max_violation.argmax(),
                               VariableState::AT_LOWER};
      case BoundViolationType::VIOLATE_UPPER_BOUND:
        return LeavingVariable{*max_violation.argmax(),
                               VariableState::AT_UPPER};
      default:
        std::unreachable();
    }
  }

  // Select a random boundary violating variable. This strategy is used when
  // potential cycling is detected.
  std::optional<LeavingVariable> get_dual_leaving_variable_randomly(
      const IterationState<Field>& state) {
    std::vector<LeavingVariable> result;

    for (size_t i = 0; i < state.basic_variables.size(); ++i) {
      auto violation = (*state.bounds)[state.basic_variables[i]].get_violation(
          state.basic_point[i, 0]);

      if (violation.type == BoundViolationType::VIOLATE_LOWER_BOUND) {
        result.emplace_back(i, VariableState::AT_LOWER);
      } else if (violation.type == BoundViolationType::VIOLATE_UPPER_BOUND) {
        result.emplace_back(i, VariableState::AT_UPPER);
      }
    }

    if (result.empty()) {
      return std::nullopt;
    }

    std::uniform_int_distribution<size_t> index_dist(0, result.size() - 1);
    size_t index = index_dist(random_engine_);

    return result.at(index);
  }

 public:
  std::optional<LeavingVariable> get(const IterationState<Field>& state) {
    if (state.last_cycling_iteration &&
        state.iteration_index <
            *state.last_cycling_iteration + kAnticyclingIterations) {
      return get_dual_leaving_variable_randomly(state);
    }

    return get_dual_leaving_variable_dantzig(state);
  }
};

}  // namespace simplex
