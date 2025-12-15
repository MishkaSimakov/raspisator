#pragma once

#include <iostream>
#include <ranges>
#include <unordered_map>
#include <vector>

#include "VariableType.h"
#include "grammar/Constraint.h"
#include "grammar/Expression.h"
#include "grammar/Variable.h"

template <typename Field>
struct VariableInfo {
  std::string name;
  VariableType type;

  Field lower_bound;
  Field upper_bound;
};

template <typename Field>
struct MILPProblem {
  std::unordered_map<std::string, VariableInfo<Field>> variables;
  std::unordered_map<std::string, Field> constants;

  std::vector<Constraint<Field>> constraints;
  Expression<Field> objective;

  //
  MILPProblem() = default;

  std::unordered_map<std::string, size_t> enumerate_variables() const {
    std::unordered_map<std::string, size_t> result;

    size_t current_index = 0;
    for (const std::string& name : variables | std::views::keys) {
      result.emplace(name, current_index);
      ++current_index;
    }

    return result;
  }

  Variable<Field> new_variable(std::string name, VariableType type,
                               Field lower_bound, Field upper_bound) {
    variables.emplace(
        name, VariableInfo<Field>{name, type, lower_bound, upper_bound});

    return Variable<Field>{name};
  }

  void add_constraint(Constraint<Field> constraint) {
    constraints.push_back(std::move(constraint));
  }

  void set_objective(Expression<Field> objective) {
    this->objective = std::move(objective);
  }

  double get_sparsity() const {
    size_t nz_count = 0;

    for (const auto& constraint : constraints) {
      nz_count += constraint.lhs.get_variables().size();
    }

    return static_cast<double>(nz_count) /
           static_cast<double>(constraints.size() * variables.size());
  }
};

template <typename Field>
std::ostream& operator<<(std::ostream& os, const MILPProblem<Field>& problem) {
  std::println(os, "problem with {} constraints", problem.constraints.size());

  os << "max " << problem.objective << "\n";
  os << "such that:\n";

  for (const Constraint<Field>& constraint : problem.constraints) {
    os << constraint << "\n";
  }

  for (const VariableInfo<Field>& info :
       problem.variables | std::views::values) {
    std::println(os, "{} <= {} <= {}", info.lower_bound, info.name,
                 info.upper_bound);
  }

  return os;
}
