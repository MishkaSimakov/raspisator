#pragma once

#include <optional>

#include "Pricing.h"
#include "simplex/StateView.h"

namespace simplex {

template <typename Field>
class PrimalSteepestEdge final : public PrimalPricing<Field> {
  // weights_[i] stores squared length of the edge corresponding to variable i
  // if i is non-basic variable
  Vector<Field> weights_;

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

  void recalculate_weights(const problem::StandardLP<Field>& problem,
                           linalg::LUPA<Field>& lupa,
                           const std::vector<VariableState>& var_states) {
    for (size_t i = 0; i < var_states.size(); ++i) {
      if (var_states[i] == VariableState::BASIC) {
        continue;
      }

      const Vector edge =
          lupa.solve_linear(problem.matrix.get_column_as_matrix(i));

      weights_[i] = 1;
      for (size_t j = 0; j < edge.size(); ++j) {
        weights_[i] += edge[j] * edge[j];
      }
    }
  }

  CyclingDetector<Field> cycling_;

 public:
  void init(const problem::StandardLP<Field>& problem,
            linalg::LUPA<Field>& lupa,
            const std::vector<VariableState>& var_states,
            const std::vector<size_t>& basic_vars) override {
    const auto [n, d] = problem.matrix.shape();

    weights_.resize(d);
    recalculate_weights(problem, lupa, var_states);
  }

  std::optional<size_t> get_primal_entering(
      StateView<Field> simplex, const Vector<Field>& reduced_cost) override {
    auto status =
        cycling_.record(simplex.iteration, simplex.states, simplex.objective);

    if (status == CyclingState::HAS_CYCLING) {
      std::cout << "cycling!!!!" << std::endl;
    }

    ArgMaximum<Field> max_cost;

    for (size_t i = 0; i < reduced_cost.size(); ++i) {
      if (!simplex.problem.var_bounds[i].is_fixed() &&
          !is_feasible(simplex.states[i], reduced_cost[i],
                       simplex.tolerance.feasibility)) {
        max_cost.record(i, reduced_cost[i] * reduced_cost[i] / weights_[i]);
      }
    }

    return max_cost.has_value() ? std::optional{max_cost->index} : std::nullopt;
  }

  void move(ChangeBasisMove<Field> move, StateView<Field> simplex,
            const Vector<Field>& pivot_row,
            const Vector<Field>& pivot_col) override {
    // w = pivot_col
    const Vector v = simplex.lupa.solve_linear_transposed(pivot_col);

    // very important, without this line error grows quickly
    weights_[move.entering_variable] = 1 + linalg::dot(pivot_col, pivot_col);

    weights_[simplex.basic_vars[move.leaving_index]] =
        weights_[move.entering_variable] /
        (pivot_row[move.entering_variable] * pivot_row[move.entering_variable]);

    for (size_t i = 0; i < simplex.problem.matrix.cols(); ++i) {
      if (simplex.states[i] == VariableState::BASIC ||
          i == move.entering_variable) {
        continue;
      }

      const Field new_pivot_row_value =
          pivot_row[i] / pivot_row[move.entering_variable];

      if (new_pivot_row_value == 0) {
        continue;
      }

      weights_[i] += new_pivot_row_value * new_pivot_row_value *
                     weights_[move.entering_variable];

      for (const auto& [row, value] : simplex.problem.matrix.get_column(i)) {
        weights_[i] -= 2 * new_pivot_row_value * value * v[row];
      }
    }
  }

  void move(ToggleBoundMove<Field> move, StateView<Field> simplex) override {}

  void post_refactorization(const problem::StandardLP<Field>& problem,
                            linalg::LUPA<Field>& lupa,
                            const std::vector<VariableState>& var_states,
                            const std::vector<size_t>& basic_vars) override {
    recalculate_weights(problem, lupa, var_states);
  }
};

}  // namespace simplex
