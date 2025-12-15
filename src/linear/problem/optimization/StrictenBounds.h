#pragma once

#include <unordered_map>

#include "BaseOptimizer.h"

template <typename Field>
class StrictenBounds final : public BaseOptimizer<Field> {
  std::pair<Field, Field> get_bounds(const MILPProblem<Field>& problem,
                                     const Expression<Field>& expression,
                                     const std::string& target) {
    Field lower = -expression.get_shift();
    Field upper = -expression.get_shift();

    for (const auto& [var, coef] : expression.get_variables()) {
      if (var == target) {
        continue;
      }

      const auto& info = problem.get_variable_info(var);

      if (coef > 0) {
        lower -= coef * info.upper_bound;
        upper -= coef * info.lower_bound;
      } else {
        lower -= coef * info.lower_bound;
        upper -= coef * info.upper_bound;
      }
    }

    return {lower, upper};
  }

 public:
  StrictenBounds() = default;

  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    for (const auto& constraint : problem.constraints) {
      if (constraint.type != ConstraintType::EQUAL_ZERO) {
        continue;
      }

      for (const auto& [var, coef] : constraint.expr.get_variables()) {
        auto [lower, upper] = get_bounds(problem, constraint.expr, var);

        auto& info = problem.get_variable_info(var);

        if (coef > 0) {
          info.lower_bound = std::max(lower / coef, info.lower_bound);
          info.upper_bound = std::min(upper / coef, info.upper_bound);
        } else {
          info.lower_bound = std::max(upper / coef, info.lower_bound);
          info.upper_bound = std::min(lower / coef, info.upper_bound);
        }
      }
    }

    return problem;
  }

  Matrix<Field> inverse(const Matrix<Field>& point) override { return point; }
};
