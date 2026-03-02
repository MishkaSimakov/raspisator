#pragma once

#include "BaseOptimizer.h"
#include "linear/problem/MILPProblem.h"

template <typename Field>
class RemoveOneVariableConstraints final : public BaseOptimizer<Field> {
 public:
  RemoveOneVariableConstraints() = default;

  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    for (auto itr = problem.constraints.begin();
         itr != problem.constraints.end();) {
      if (itr->expr.get_variables().size() != 1) {
        ++itr;
        continue;
      }

      auto [var, coef] = *itr->expr.get_variables().begin();
      auto& info = problem.get_variable_info(var);

      Field value = -itr->expr.get_shift() / coef;

      if (itr->type == ConstraintType::EQUAL_ZERO) {
        if (!info.bound.is_inside(value)) {
          throw std::runtime_error("Problem is trivially unfeasible.");
        }

        info.bound.lower = value;
        info.bound.upper = value;
      } else {
        if (coef > 0) {
          info.bound.stricten_upper(value);
        } else {
          info.bound.stricten_lower(value);
        }
      }

      itr = problem.constraints.erase(itr);
    }

    return problem;
  }

  Matrix<Field> inverse(const Matrix<Field>& point) override { return point; }
};
