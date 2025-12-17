#pragma once

#include <iostream>
#include <ranges>
#include <unordered_map>
#include <vector>

#include "VariableType.h"
#include "grammar/Constraint.h"
#include "grammar/Expression.h"
#include "grammar/Variable.h"
#include "linear/matrix/Matrix.h"
#include "utils/Accumulators.h"

template <typename Field>
struct VariableInfo {
  std::string name;
  VariableType type;

  Field lower_bound;
  Field upper_bound;
};

template <typename Field>
struct MILPProblem {
  std::vector<VariableInfo<Field>> variables;

  std::vector<Constraint<Field>> constraints;
  Expression<Field> objective;

  //
  MILPProblem() = default;

  std::unordered_map<std::string, size_t> enumerate_variables() const {
    std::unordered_map<std::string, size_t> result;

    for (size_t i = 0; i < variables.size(); ++i) {
      result.emplace(variables[i].name, i);
    }

    return result;
  }

  VariableInfo<Field>& get_variable_info(const std::string& name) {
    for (VariableInfo<Field>& info : variables) {
      if (info.name == name) {
        return info;
      }
    }

    throw std::runtime_error(std::format("No variables with name {}.", name));
  }

  const VariableInfo<Field>& get_variable_info(const std::string& name) const {
    for (const VariableInfo<Field>& info : variables) {
      if (info.name == name) {
        return info;
      }
    }

    throw std::runtime_error(std::format("No variables with name {}.", name));
  }

  Variable<Field> get_variable(const std::string& name) const {
    return Variable<Field>{name};
  }

  Field extract_variable(const Variable<Field>& variable,
                         const Matrix<Field>& point) {
    std::optional<size_t> index;

    for (size_t i = 0; i < variables.size(); ++i) {
      if (variables[i].name == variable.get_name()) {
        index = i;
        break;
      }
    }

    if (!index) {
      throw std::runtime_error(std::format(
          "No variable with name {} in the problem.", variable.get_name()));
    }

    return point[*index, 0];
  }

  Variable<Field> new_variable(std::string name, VariableType type,
                               Field lower_bound, Field upper_bound) {
    for (const VariableInfo<Field>& info : variables) {
      if (info.name == name) {
        throw std::runtime_error(
            std::format("There is already a variable with name {}.", name));
      }
    }

    variables.emplace_back(name, type, lower_bound, upper_bound);

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

  Field average_boundary_gap() const {
    ArithmeticMean<Field> gap;

    for (const auto& variable : variables) {
      gap.record(variable.upper_bound - variable.lower_bound);
    }

    return gap.mean();
  }
};

template <typename Field>
std::ostream& operator<<(std::ostream& os, const MILPProblem<Field>& problem) {
  std::println(os, "problem with {} constraints and {} variables",
               problem.constraints.size(), problem.variables.size());

  os << "max " << problem.objective << "\n";
  os << "such that:\n";

  for (const Constraint<Field>& constraint : problem.constraints) {
    os << constraint << "\n";
  }

  for (const VariableInfo<Field>& info : problem.variables) {
    std::println(os, "{} <= {} <= {}", info.lower_bound, info.name,
                 info.upper_bound);
  }

  return os;
}
