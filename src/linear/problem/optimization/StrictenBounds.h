#pragma once

#include <unordered_map>

#include "BaseOptimizer.h"

template <typename Field>
class StrictenBounds final : public BaseOptimizer<Field> {
  Bound<Field> get_bound(const MILPProblem<Field>& problem,
                         const Expression<Field>& expression) {
    Bound<Field> result(expression.get_shift(), expression.get_shift());

    for (const auto& [var, coef] : expression.get_variables()) {
      result += coef * problem.get_variable_info(var).bound;
    }

    return result;
  }

 public:
  StrictenBounds() = default;

  // Note: This optimizer does not change infinite bounds!
  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    for (const auto& constraint : problem.constraints) {
      if (constraint.type != ConstraintType::EQUAL_ZERO) {
        continue;
      }

      auto constraint_bound = get_bound(problem, constraint.expr);

      for (const auto& [var, coef] : constraint.expr.get_variables()) {
        auto new_bound = -constraint_bound / coef;
        auto& curr_bound = problem.get_variable_info(var).bound;

        if (curr_bound.lower) {
          curr_bound.stricten_lower(new_bound.lower);
        }
        if (curr_bound.upper) {
          curr_bound.stricten_upper(new_bound.upper);
        }
      }
    }

    for (auto& info : problem.variables) {
      if (info.type == VariableType::INTEGER) {
        if (info.bound.lower) {
          info.bound.lower = -FieldTraits<Field>::floor(-*info.bound.lower);
        }

        if (info.bound.upper) {
          info.bound.upper = FieldTraits<Field>::floor(*info.bound.upper);
        }
      }

      if (info.bound.is_infeasible()) {
        throw std::runtime_error("Problem is trivially infeasible.");
      }
    }

    return problem;
  }

  Matrix<Field> inverse(const Matrix<Field>& point) override { return point; }
};
