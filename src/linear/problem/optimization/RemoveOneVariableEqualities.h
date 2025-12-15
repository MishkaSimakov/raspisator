#pragma once

#include "BaseOptimizer.h"
#include "linear/problem/MILPProblem.h"

template <typename Field>
class RemoveOneVariableEqualities final : public BaseOptimizer<Field> {
  void replace_in_constraints(MILPProblem<Field>& problem,
                              const std::string& name, Field value) {
    for (auto& constraint : problem.constraints) {
      auto& variables = constraint.lhs.get_variables();
      auto itr = constraint.lhs.get_by_name(name);

      if (itr != variables.end()) {
        constraint.rhs -= itr->second * value;
        variables.erase(itr);
      }
    }
  }

  void replace_in_objective(MILPProblem<Field>& problem,
                            const std::string& name, Field value) {
    auto& variables = problem.objective.get_variables();
    auto itr = problem.objective.get_by_name(name);

    if (itr != variables.end()) {
      problem.objective -= itr->second * value;
      variables.erase(itr);
    }
  }

 public:
  RemoveOneVariableEqualities() = default;

  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    for (auto itr = problem.constraints.begin();
         itr != problem.constraints.end();) {
      if (itr->type == ConstraintType::EQUAL &&
          itr->lhs.get_variables().size() == 1) {
        Field value = itr->rhs / itr->lhs.get_variables().begin()->second;
        auto& info =
            problem.variables.at(itr->lhs.get_variables().begin()->first);

        if (FieldTraits<Field>::is_strictly_negative(value -
                                                     info.lower_bound) ||
            FieldTraits<Field>::is_strictly_negative(info.upper_bound -
                                                     value)) {
          throw std::runtime_error("Problem is trivially unfeasible.");
        }

        info.lower_bound = value;
        info.upper_bound = value;

        itr = problem.constraints.erase(itr);
      } else {
        ++itr;
      }
    }

    return problem;
  }
};
