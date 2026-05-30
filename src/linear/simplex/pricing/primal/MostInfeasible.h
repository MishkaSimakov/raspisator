#pragma once

#include <random>

#include "Pricing.h"
#include "linear/simplex/CyclingDetector.h"
#include "utils/Accumulators.h"

namespace simplex {

template <typename Field>
class PrimalMostInfeasible final : public PrimalPricing<Field> {
  // most infeasible pricing is prone to cycling, therefore cycling detection is
  // used. When cycling is detected, strategy is swapped to random infeasible
  // pricing.
  CyclingDetector<Field> cycling_;

  std::default_random_engine random_;

  static bool is_feasible(VariableState state, Field reduced_cost,
                          Field tolerance) {
    using std::abs;

    switch (state) {
      case VariableState::BASIC:
        return true;
      case VariableState::AT_LOWER:
        return reduced_cost <= tolerance;
      case VariableState::AT_UPPER:
        return reduced_cost >= -tolerance;
      case VariableState::NONBASIC_FREE:
        return abs(reduced_cost) <= tolerance;
      default:
        throw std::runtime_error("Unknown variable state.");
    }
  }

  std::optional<size_t> random_pricing(detail::State<Field> simplex) {
    using std::abs;

    // TODO: Reservoir sampling
    // choose random among top k by reduced cost
    std::vector<std::pair<double, size_t>> costs;

    for (size_t i = 0; i < simplex.reduced_cost.size(); ++i) {
      if (!is_feasible(simplex.states[i], simplex.reduced_cost[i],
                       simplex.tolerance.feasibility)) {
        costs.emplace_back(abs(simplex.reduced_cost[i]), i);
      }
    }

    if (costs.empty()) {
      return std::nullopt;
    }

    std::ranges::sort(costs, {}, [](auto p) { return -p.first; });

    const size_t count = std::min(5uz, costs.size());
    const size_t index = random_() % count;

    return costs[index].second;
  }

  std::optional<size_t> most_infeasible_pricing(detail::State<Field> simplex) {
    ArgMaximum<Field> max_cost;

    for (size_t i = 0; i < simplex.reduced_cost.size(); ++i) {
      if (!is_feasible(simplex.states[i], simplex.reduced_cost[i],
                       simplex.tolerance.feasibility)) {
        max_cost.record(i, abs(simplex.reduced_cost[i]));
      }
    }

    return max_cost.has_value() ? std::optional{max_cost->index} : std::nullopt;
  }

 public:
  std::optional<size_t> get_primal_entering(
      detail::State<Field> simplex) override {
    if (cycling_.record(simplex.iteration, simplex.states, simplex.objective) ==
        CyclingState::HAS_CYCLING) {
      return random_pricing(simplex);
    }

    return most_infeasible_pricing(simplex);
  }
};

}  // namespace simplex
