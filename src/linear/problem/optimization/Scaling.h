#pragma once

#include <cmath>

#include "BaseOptimizer.h"
#include "linear/matrix/Matrix.h"
#include "linear/problem/MILPProblem.h"
#include "utils/Accumulators.h"

template <typename Field>
class Scaling final : public BaseOptimizer<Field> {
  static_assert(std::is_floating_point_v<Field>,
                "Currently works only for standard types.");

  std::vector<Field> variables_scale_factors_;

  static Field round_scale_factor(Field scale_factor) {
    double power = std::round(std::log2(scale_factor));
    return std::exp2(-power);
  }

  static Field get_scale_factor(const Constraint<Field>& constraint) {
    GeometricAverage<Field> scale_factor;

    for (Field coef : constraint.lhs.get_variables() | std::views::values) {
      scale_factor.record(FieldTraits<Field>::abs(coef));
    }

    return round_scale_factor(scale_factor.average());
  }

  static Field get_scale_factor(const MILPProblem<Field>& problem,
                                const std::string& variable) {
    GeometricAverage<Field> scale_factor;

    for (const auto& constraint : problem.constraints) {
      auto& variables = constraint.lhs.get_variables();
      auto itr = variables.find(variable);

      if (itr != variables.end()) {
        scale_factor.record(FieldTraits<Field>::abs(itr->second));
      }
    }

    return round_scale_factor(scale_factor.average());
  }

 public:
  Scaling() = default;

  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    variables_scale_factors_.resize(problem.variables.size());

    // scale rows
    for (auto& constraint : problem.constraints) {
      Field scale_factor = get_scale_factor(constraint);

      for (Field& coef : constraint.lhs.get_variables() | std::views::values) {
        coef *= scale_factor;
      }
      constraint.rhs *= scale_factor;
    }

    // scale columns
    for (size_t i = 0; i < problem.variables.size(); ++i) {
      auto& info = problem.variables[i];
      Field scale_factor = get_scale_factor(problem, info.name);

      variables_scale_factors_[i] = scale_factor;

      for (auto& constraint : problem.constraints) {
        auto& variables = constraint.lhs.get_variables();
        auto itr = variables.find(info.name);

        if (itr != variables.end()) {
          itr->second *= scale_factor;
        }
      }

      info.lower_bound /= scale_factor;
      info.upper_bound /= scale_factor;

      auto& objective_variables = problem.objective.get_variables();
      auto itr = objective_variables.find(info.name);

      if (itr != objective_variables.end()) {
        itr->second *= scale_factor;
      }
    }

    return problem;
  }

  // Well scaled if returned value < 2
  // https://pure.iiasa.ac.at/id/eprint/4172/7/WP-94-037.pdf
  double get_scaling_quality(const MILPProblem<Field>& problem) const {
    Minimum<Field> min;
    Maximum<Field> max;

    for (const auto& constraint : problem.constraints) {
      for (const auto& [var, coef] : constraint.lhs.get_variables()) {
        Field abs = FieldTraits<Field>::abs(coef);

        min.record(abs);
        max.record(abs);
      }
    }

    return std::log10(*max.max() / *min.min());
  }

  Matrix<Field> inverse(const Matrix<Field>& point) override {
    auto result = point;

    for (size_t i = 0; i < result.get_height(); ++i) {
      result[i, 0] *= variables_scale_factors_[i];
    }

    return result;
  }
};
