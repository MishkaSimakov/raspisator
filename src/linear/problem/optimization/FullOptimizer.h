#pragma once

#include "BaseOptimizer.h"
#include "ObjectivePerturbation.h"
#include "RemoveConstantAndUnusedVariables.h"
#include "RemoveConstantConstraints.h"
#include "RemoveLinearlyDependentConstraints.h"
#include "RemoveOneVariableConstraints.h"
#include "Scaling.h"
#include "StrictenBounds.h"
#include "Substitution.h"
#include "TransformToEqualities.h"

template <typename Field>
class FullOptimizer final : public BaseOptimizer<Field> {
  bool log_;

  std::vector<std::unique_ptr<BaseOptimizer<Field>>> optimizers_chain_;

  template <typename T>
  MILPProblem<Field> apply(const MILPProblem<Field>& problem) {
    auto& optimizer = optimizers_chain_.emplace_back(std::make_unique<T>());

    auto result = optimizer->apply(problem);

    if (log_) {
      std::println("{}: {} variables, {} constraints", typeid(T).name(),
                   static_cast<ssize_t>(result.variables.size()) -
                       static_cast<ssize_t>(problem.variables.size()),
                   static_cast<ssize_t>(result.constraints.size()) -
                       static_cast<ssize_t>(problem.constraints.size()));
    }

    return result;
  }

  MILPProblem<Field> apply_main_loop(MILPProblem<Field> problem) {
    problem = apply<RemoveConstantConstraints<Field>>(problem);
    problem = apply<RemoveLinearlyDependentConstraints<Field>>(problem);
    problem = apply<RemoveConstantAndUnusedVariables<Field>>(problem);
    problem = apply<RemoveOneVariableConstraints<Field>>(problem);
    problem = apply<StrictenBounds<Field>>(problem);
    problem = apply<Substitution<Field>>(problem);

    return problem;
  }

  MILPProblem<Field> scale(MILPProblem<Field> problem) {
    Field initial_quality = Scaling<Field>().get_scaling_quality(problem);

    auto result = apply<Scaling<Field>>(problem);

    if (log_) {
      Field new_quality = Scaling<Field>().get_scaling_quality(result);
      std::println("scaling quality: {} -> {}", initial_quality, new_quality);
    }

    return result;
  }

 public:
  explicit FullOptimizer(bool log = false) : log_(log) {}

  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    problem = apply_main_loop(problem);
    problem = apply_main_loop(problem);

    problem = apply<TransformToEqualities<Field>>(problem);

    for (size_t i = 0; i < 10; ++i) {
      problem = apply_main_loop(problem);
    }

    problem = scale(problem);
    problem = apply<ObjectivePerturbation<Field>>(problem);

    return problem;
  }

  Matrix<Field> inverse(const Matrix<Field>& point) override {
    auto result = point;

    for (std::unique_ptr<BaseOptimizer<Field>>& optimizer :
         optimizers_chain_ | std::views::reverse) {
      result = optimizer->inverse(result);
    }

    return result;
  }
};
