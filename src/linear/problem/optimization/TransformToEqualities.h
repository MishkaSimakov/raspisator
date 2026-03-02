#pragma once

#include "BaseOptimizer.h"
#include "linear/problem/MILPProblem.h"

template <typename Field>
class TransformToEqualities final : public BaseOptimizer<Field> {
  size_t slack_count_{0};

  size_t slack_variables_count(const MILPProblem<Field>& problem) {
    size_t result = 0;

    for (const auto& info : problem.variables) {
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
      if (constraint.type == ConstraintType::EQUAL_ZERO) {
        continue;
      }

      ++slack_count_;

      auto var_name = std::format("slack({})", slack_variables_count(problem));
      auto slack_var =
          problem.new_variable(var_name, VariableType::SLACK, 0, std::nullopt);

      constraint.expr += slack_var;
      constraint.type = ConstraintType::EQUAL_ZERO;
    }

    return problem;
  }

  Matrix<Field> inverse(const Matrix<Field>& point) override {
    return point[{0, point.get_height() - slack_count_}, 0];
  }
};
