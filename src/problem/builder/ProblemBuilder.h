#pragma once

#include <iostream>
#include <ranges>
#include <unordered_map>
#include <vector>

#include "VariableType.h"
#include "grammar/Constraint.h"
#include "grammar/Expression.h"
#include "grammar/Variable.h"
#include "linalg/Linalg.h"
#include "linear/model/Bound.h"
#include "problem/MILP.h"
#include "utils/Accumulators.h"

namespace problem {

template <typename Field>
struct VariableInfo {
  std::string name;
  VariableType type;

  Bound<Field> bound;
};

template <typename Field>
class ProblemBuilder {
  std::vector<VariableInfo<Field>> variables_;

  std::vector<Constraint<Field>> constraints_;
  Expression<Field> objective_;

 public:
  ProblemBuilder() = default;

  std::unordered_map<std::string, size_t> enumerate_variables() const {
    std::unordered_map<std::string, size_t> result;

    for (size_t i = 0; i < variables_.size(); ++i) {
      result.emplace(variables_[i].name, i);
    }

    return result;
  }

  VariableInfo<Field>& get_variable_info(const std::string& name) {
    for (VariableInfo<Field>& info : variables_) {
      if (info.name == name) {
        return info;
      }
    }

    throw std::runtime_error(std::format("No variables with name {}.", name));
  }

  const VariableInfo<Field>& get_variable_info(const std::string& name) const {
    for (const VariableInfo<Field>& info : variables_) {
      if (info.name == name) {
        return info;
      }
    }

    throw std::runtime_error(std::format("No variables with name {}.", name));
  }

  Variable<Field> get_variable(const std::string& name) const {
    return Variable<Field>{name};
  }

  Variable<Field> new_variable(std::string name, VariableType type,
                               Bound<Field> bound) {
    for (const VariableInfo<Field>& info : variables_) {
      if (info.name == name) {
        throw std::runtime_error(
            std::format("There is already a variable with name {}.", name));
      }
    }

    variables_.emplace_back(name, type, bound);

    return Variable<Field>{name};
  }

  Variable<Field> new_variable(std::string name, VariableType type,
                               std::optional<Field> lower_bound,
                               std::optional<Field> upper_bound) {
    return new_variable(name, type, Bound<Field>(lower_bound, upper_bound));
  }

  void add_constraint(Constraint<Field> constraint) {
    constraints_.push_back(std::move(constraint));
  }

  void set_objective(Expression<Field> objective) {
    this->objective_ = std::move(objective);
  }

  double get_sparsity() const {
    size_t nz_count = 0;

    for (const auto& constraint : constraints_) {
      nz_count += constraint.lhs.get_variables().size();
    }

    return static_cast<double>(nz_count) /
           static_cast<double>(constraints_.size() * variables_.size());
  }

  Field average_boundary_gap() const {
    ArithmeticMean<Field> result;

    for (const auto& variable : variables_) {
      if (variable.bound.upper && variable.bound.lower) {
        result.record(*variable.bound.upper - *variable.bound.lower);
      }
    }

    return *result;
  }

  explicit operator MILP<Field>() const {
    const size_t n = constraints_.size();
    const size_t d = variables_.size();

    const auto enumeration = enumerate_variables();

    MILP<Field> result;

    // cost
    result.cost.resize(d);
    for (const auto& [var, coef] : objective_.get_variables()) {
      result.cost[enumeration.at(var)] = coef;
    }

    // rhs bounds
    result.rhs_bounds.resize(n);
    for (size_t row = 0; row < n; ++row) {
      const Field rhs = -constraints_[row].expr.get_shift();

      if (constraints_[row].type == ConstraintType::EQUAL_ZERO) {
        result.rhs_bounds[row] = {rhs, rhs};
      } else {
        result.rhs_bounds[row] = {std::nullopt, rhs};
      }
    }

    // constraints matrix
    result.matrix.resize(n, 0);

    for (size_t col = 0; col < d; ++col) {
      result.matrix.add_column();

      for (size_t row = 0; row < n; ++row) {
        const auto& vars = constraints_[row].expr.get_variables();

        if (vars.contains(variables_[col].name)) {
          result.matrix.push_to_last_column(row, vars.at(variables_[col].name));
        }
      }
    }

    // var info
    result.var_bounds.resize(d);
    result.is_integer.resize(d);
    result.var_names.resize(d);

    for (const auto& info : variables_) {
      const size_t index = enumeration.at(info.name);

      result.var_bounds[index] = info.bound;
      result.is_integer[index] = info.type == VariableType::INTEGER;
      result.var_names[index] = info.name;
    }

    result.implied_is_integer = result.is_integer;
    result.implied_var_bounds = result.var_bounds;

    return result;
  }
};

template <typename Field>
std::ostream& operator<<(std::ostream& os, const MILPProblem<Field>& problem) {
  std::println(os, "problem with {} constraints and {} variables",
               problem.constraints_.size(), problem.variables_.size());

  os << "max " << problem.objective_ << "\n";
  os << "such that:\n";

  for (const Constraint<Field>& constraint : problem.constraints_) {
    os << constraint << "\n";
  }

  for (const VariableInfo<Field>& info : problem.variables_) {
    std::println(os, "{} ∈ {}", info.name, info.bound);
  }

  return os;
}

}  // namespace problem
