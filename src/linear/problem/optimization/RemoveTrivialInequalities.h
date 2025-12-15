#pragma once

#include "BaseOptimizer.h"
#include "linear/problem/MILPProblem.h"

template <typename Field>
class RemoveTrivialInequalities final : public BaseOptimizer<Field> {
 public:
  RemoveTrivialInequalities() = default;

  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    for (auto itr = problem.constraints.begin();
         itr != problem.constraints.end();) {
      if (itr->type == ConstraintType::LESS_OR_EQUAL &&
          itr->lhs.get_variables().size() == 1) {
        auto [var, coef] = *itr->lhs.get_variables().begin();
        auto& info = problem.get_variable(var);

        if (coef > 0) {
          info.upper_bound = std::min(info.upper_bound, itr->rhs / coef);
        } else {
          info.lower_bound = std::max(info.lower_bound, itr->rhs / coef);
        }

        itr = problem.constraints.erase(itr);
      } else {
        ++itr;
      }
    }

    return problem;
  }

  Matrix<Field> inverse(const Matrix<Field>& point) override { return point; }
};
