#pragma once

#include "BaseOptimizer.h"
#include "linear/problem/MILPProblem.h"

template <typename Field>
class ObjectivePerturbation final : public BaseOptimizer<Field> {
 public:
  ObjectivePerturbation() = default;

  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    std::vector<std::string> perturbed;

    for (const auto& variable : problem.variables) {
      if (problem.objective.get_variables().contains(variable.name)) {
        continue;
      }

      if (!variable.bound.lower ||
          FieldTraits<Field>::is_nonzero(*variable.bound.lower)) {
        continue;
      }

      if (variable.bound.is_fixed()) {
        continue;
      }

      // remove variables with big upper bounds so that coefficients are not
      // too small
      if (!variable.bound.upper || variable.bound.upper > 500) {
        continue;
      }

      problem.objective -= static_cast<Field>(1) / *variable.bound.upper *
                           problem.get_variable(variable.name);
      perturbed.push_back(variable.name);
    }

    for (const std::string& name : perturbed) {
      problem.objective.get_variables().at(name) /=
          static_cast<Field>(2 * perturbed.size());
    }

    return problem;
  }

  Matrix<Field> inverse(const Matrix<Field>& point) override { return point; }
};
