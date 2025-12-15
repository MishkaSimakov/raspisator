#pragma once

#include "BaseOptimizer.h"
#include "linear/problem/MILPProblem.h"

template <typename Field>
class TransformToEqualities final : public BaseOptimizer<Field> {
  Field get_slack_upper_bound(const MILPProblem<Field>& problem,
                              const Constraint<Field>& constraint) const {
    Field result = constraint.rhs;

    for (auto [var, coef] : constraint.lhs.get_variables()) {
      const auto& info = problem.variables.at(var);

      if (coef > 0) {
        result -= coef * info.lower_bound;
      } else {
        result -= coef * info.upper_bound;
      }
    }

    return result;
  }

  size_t slack_variables_count(const MILPProblem<Field>& problem) {
    size_t result = 0;

    for (const auto& info : problem.variables | std::views::values) {
      if (info.type == VariableType::SLACK) {
        ++result;
      }
    }

    return result;
  }

 public:
  TransformToEqualities() = default;

  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    for (auto& constraint : problem.constraints) {
      if (constraint.type == ConstraintType::EQUAL) {
        continue;
      }

      auto var_name = std::format("slack({})", slack_variables_count(problem));
      Field upper_bound = get_slack_upper_bound(problem, constraint);
      auto slack_var =
          problem.new_variable(var_name, VariableType::SLACK, 0, upper_bound);

      constraint.lhs += slack_var;
      constraint.type = ConstraintType::EQUAL;
    }

    return problem;
  }
};
