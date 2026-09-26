#pragma once

#include <optional>
#include <random>

#include "Pricing.h"
#include "simplex/CyclingDetector.h"
#include "simplex/StateView.h"
#include "utils/Accumulators.h"

namespace simplex {

template <typename Field>
class DualDantzigPricing final : public DualPricing<Field> {
  CyclingDetector<Field> cycling_;

  std::default_random_engine random_;

  std::optional<LeavingVariable> dantzig_pricing(StateView<Field> simplex) {
    ArgMaximum<BoundViolation<Field>> max_violation;

    for (size_t i = 0; i < simplex.basic_vars.size(); ++i) {
      max_violation.record(
          i, simplex.problem.var_bounds[simplex.basic_vars[i]].get_violation(
                 simplex.basic_point[i]));
    }

    if (!max_violation.has_value()) {
      return std::nullopt;
    }

    switch (max_violation->max.type) {
      case BoundViolationType::NONE:
        return std::nullopt;
      case BoundViolationType::VIOLATE_LOWER_BOUND:
        return LeavingVariable{max_violation->index, VariableState::AT_LOWER};
      case BoundViolationType::VIOLATE_UPPER_BOUND:
        return LeavingVariable{max_violation->index, VariableState::AT_UPPER};
      default:
        std::unreachable();
    }
  }

  // Select a random boundary violating variable. This strategy is used when
  // potential cycling is detected.
  std::optional<LeavingVariable> random_pricing(StateView<Field> state) {
    std::vector<LeavingVariable> result;

    for (size_t i = 0; i < state.basic_vars.size(); ++i) {
      auto violation =
          state.problem.var_bounds[state.basic_vars[i]].get_violation(
              state.basic_point[i]);

      if (violation.type == BoundViolationType::VIOLATE_LOWER_BOUND) {
        result.emplace_back(i, VariableState::AT_LOWER);
      } else if (violation.type == BoundViolationType::VIOLATE_UPPER_BOUND) {
        result.emplace_back(i, VariableState::AT_UPPER);
      }
    }

    if (result.empty()) {
      return std::nullopt;
    }

    return result.at(random_() % result.size());
  }

 public:
  std::optional<LeavingVariable> get_dual_leaving(
      StateView<Field> simplex) override {
    if (!simplex.intentional_repeat &&
        cycling_.record(simplex.iteration, simplex.states, simplex.objective) ==
            CyclingState::HAS_CYCLING) {
      return random_pricing(simplex);
    }

    return dantzig_pricing(simplex);
  }
};

}  // namespace simplex
