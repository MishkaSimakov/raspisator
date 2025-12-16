#pragma once

#include "BaseOptimizer.h"
#include "linear/problem/MILPProblem.h"

template <typename Field>
class RemoveConstantAndUnusedVariables final : public BaseOptimizer<Field> {
  std::unordered_map<size_t, Field> constants_;

  void replace_in_constraints(MILPProblem<Field>& problem,
                              const std::string& name, Field value) {
    for (auto& constraint : problem.constraints) {
      auto& variables = constraint.expr.get_variables();
      auto itr = variables.find(name);

      if (itr != variables.end()) {
        constraint.expr += itr->second * value;
        variables.erase(itr);
      }
    }
  }

  void replace_in_objective(MILPProblem<Field>& problem,
                            const std::string& name, Field value) {
    auto& variables = problem.objective.get_variables();
    auto itr = variables.find(name);

    if (itr != variables.end()) {
      problem.objective += itr->second * value;
      variables.erase(itr);
    }
  }

  bool is_unused(const MILPProblem<Field>& problem,
                 const std::string& variable) {
    for (auto& constraint : problem.constraints) {
      if (constraint.expr.get_variables().contains(variable)) {
        return false;
      }
    }

    if (problem.objective.get_variables().contains(variable)) {
      return false;
    }

    return true;
  }

 public:
  RemoveConstantAndUnusedVariables() = default;

  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    size_t i = 0;

    for (auto itr = problem.variables.begin();
         itr != problem.variables.end();) {
      if (is_unused(problem, itr->name)) {
        // unused variable may be equal to any value within its bounds
        Field value = itr->upper_bound;
        constants_.emplace(i, value);

        itr = problem.variables.erase(itr);
      } else if (FieldTraits<Field>::is_nonzero(itr->upper_bound -
                                                itr->lower_bound)) {
        ++itr;
      } else {
        Field value = itr->upper_bound;
        constants_.emplace(i, value);

        replace_in_constraints(problem, itr->name, value);
        replace_in_objective(problem, itr->name, value);

        itr = problem.variables.erase(itr);
      }

      ++i;
    }

    return problem;
  }

  Matrix<Field> inverse(const Matrix<Field>& point) override {
    size_t d = point.get_height();

    Matrix<Field> result(d + constants_.size(), 1);

    size_t j = 0;
    for (size_t i = 0; i < d; ++i) {
      if (constants_.contains(i)) {
        result[i, 0] = constants_.at(i);
      } else {
        result[i, 0] = point[j, 0];
        ++j;
      }
    }

    return result;
  }
};
