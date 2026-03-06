#pragma once

#include <optional>
#include <random>

#include "Types.h"

namespace simplex {

template <typename Field>
class PrimalEnteringVariable {
  constexpr static size_t kAnticyclingIterations = 5;

  std::mt19937 random_engine_;

  std::optional<size_t> get_primal_entering_variable(
      const IterationState<Field>& state) {
    auto [n, d] = state.problem_shape();

    ArgMaximum<Field> max_cost;

    for (size_t i = 0; i < d; ++i) {
      if (state.variables_states[i] == VariableState::BASIC) {
        continue;
      }

      Field cost = state.reduced_cost[i, 0];
      if (state.variables_states[i] == VariableState::AT_UPPER) {
        cost *= -1;
      }

      if (FieldTraits<Field>::is_strictly_positive(cost)) {
        max_cost.record(i, cost);
      }
    }

    return max_cost.argmax();
  }

  std::optional<size_t> get_primal_entering_variable_randomly(
      const IterationState<Field>& state) {
    std::vector<size_t> result;

    auto [n, d] = state.problem_shape();

    for (size_t i = 0; i < d; ++i) {
      if (state.variables_states[i] == VariableState::BASIC) {
        continue;
      }

      Field cost = state.reduced_cost[i, 0];
      if (state.variables_states[i] == VariableState::AT_UPPER) {
        cost *= -1;
      }

      if (FieldTraits<Field>::is_strictly_positive(cost)) {
        result.push_back(i);
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
  std::optional<size_t> get(const IterationState<Field>& state) {
    if (state.last_cycling_iteration &&
        state.iteration_index <
            *state.last_cycling_iteration + kAnticyclingIterations) {
      return get_primal_entering_variable_randomly(state);
    }

    return get_primal_entering_variable(state);
  }
};

}  // namespace simplex
