#pragma once

#include <unordered_map>

#include "BaseOptimizer.h"

template <typename Field>
class StrictenBounds final : public BaseOptimizer<Field> {
  Bound<Field> get_bound(const MILPProblem<Field>& problem,
                         const Expression<Field>& expression,
                         const std::string& target) {
    Bound<Field> result(expression.get_shift(), expression.get_shift());
    Field target_coef = 0;

    for (const auto& [var, coef] : expression.get_variables()) {
      if (var != target) {
        result += coef * problem.get_variable_info(var).bound;
      } else {
        target_coef = coef;
      }
    }

    return -result / target_coef;
  }

 public:
  StrictenBounds() = default;

  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    for (const auto& constraint : problem.constraints) {
      if (constraint.type != ConstraintType::EQUAL_ZERO) {
        continue;
      }

      for (const auto& [var, coef] : constraint.expr.get_variables()) {
        problem.get_variable_info(var).bound ^=
            get_bound(problem, constraint.expr, var);
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
