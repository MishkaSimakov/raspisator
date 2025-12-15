#pragma once

#include "BaseOptimizer.h"
#include "linear/problem/MILPProblem.h"

template <typename Field>
class RemoveConstantConstraints final : public BaseOptimizer<Field> {
 public:
  RemoveConstantConstraints() = default;

  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    for (auto itr = problem.constraints.begin();
         itr != problem.constraints.end();) {
      switch (itr->evaluate()) {
        case EvaluationResult::TRUE:
          itr = problem.constraints.erase(itr);
          break;
        case EvaluationResult::FALSE:
          throw std::runtime_error(
              "Problem is trivially unfeasible. It contains a constant "
              "constraint that evaluates to false.");
        case EvaluationResult::DONT_KNOW:
          ++itr;
          break;
        default:
          std::unreachable();
      }
    }

    return problem;
  }
};
