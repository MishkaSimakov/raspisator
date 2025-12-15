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
      auto& info = problem.get_variable(var);

      Field value = -itr->expr.get_shift() / coef;

      if (itr->type == ConstraintType::EQUAL_ZERO) {
        if (FieldTraits<Field>::is_strictly_negative(value -
                                                     info.lower_bound) ||
            FieldTraits<Field>::is_strictly_negative(info.upper_bound -
                                                     value)) {
          throw std::runtime_error("Problem is trivially unfeasible.");
        }

        info.lower_bound = value;
        info.upper_bound = value;
      } else {
        if (coef > 0) {
          info.upper_bound = std::min(info.upper_bound, value);
        } else {
          info.lower_bound = std::max(info.lower_bound, value);
        }
      }

      itr = problem.constraints.erase(itr);
    }

    return problem;
  }

  Matrix<Field> inverse(const Matrix<Field>& point) override { return point; }
};
