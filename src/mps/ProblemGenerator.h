#pragma once

#include <cmath>

#include "Types.h"
#include "linear/problem/MILPProblem.h"

namespace mps {

template <typename Field>
class ProblemGenerator {
  static Bound<Field> get_bound(const Variable<Field>& variable) {
    if (!variable.lower_specified && !variable.upper_specified) {
      return variable.is_integer ? Bound<Field>{0, 1}
                                 : Bound<Field>{0, std::nullopt};
    }

    if (!variable.lower_specified && variable.upper_specified) {
      return variable.bound.upper && *variable.bound.upper < 0
                 ? Bound<Field>(std::nullopt, variable.bound.upper)
                 : Bound<Field>(0, variable.bound.upper);
    }

    if (variable.lower_specified && !variable.upper_specified) {
      return Bound<Field>(variable.bound.lower, std::nullopt);
    }

    // both bounds are specified
    return variable.bound;
  }

  static size_t get_objective_row(const MPSParsingState<Field>& state) {
    for (size_t i = 0; i < state.rows.size(); ++i) {
      if (state.rows[i].type == RowSense::FREE) {
        return i;
      }
    }

    throw std::runtime_error(
        "Objective row was not specified in the MPS file.");
  }

  static Expression<Field> get_row_expr(
      const MPSParsingState<Field>& state, size_t row,
      const std::vector<::Variable<Field>>& variables) {
    Expression<Field> result;

    for (size_t i = 0; i < state.cols.size(); ++i) {
      auto itr = state.cols[i].values.find(row);

      if (itr != state.cols[i].values.end()) {
        result += itr->second * variables[i];
      }
    }

    result -= state.rows[row].rhs.value_or(0);

    return result;
  }

 public:
  static MILPProblem<Field> generate(const MPSParsingState<Field>& state) {
    using std::abs;

    MILPProblem<Field> result;
    std::vector<::Variable<Field>> variables;

    for (const Variable<Field>& var : state.cols) {
      const auto variable_type =
          var.is_integer ? VariableType::INTEGER : VariableType::REAL;

      variables.push_back(
          result.new_variable(var.name, variable_type, get_bound(var)));
    }

    const size_t objective_row_index = get_objective_row(state);

    auto objective_expr = get_row_expr(state, objective_row_index, variables);

    if (state.objective == ObjectiveType::MAXIMIZE) {
      std::cerr << "MPS objective is MAXIMIZE, negating objective value."
                << std::endl;

      objective_expr *= -1;
    }

    result.set_objective(objective_expr);

    // process constraints
    for (size_t i = 0; i < state.rows.size(); ++i) {
      if (state.rows[i].type == RowSense::FREE) {
        continue;
      }

      const auto row = get_row_expr(state, i, variables);

      if (!state.rows[i].range) {
        if (state.rows[i].type == RowSense::LESS_THAN) {
          result.add_constraint(row <= Expression<Field>{0});
        } else if (state.rows[i].type == RowSense::GREATER_THAN) {
          result.add_constraint(row >= Expression<Field>{0});
        } else {
          result.add_constraint(row == Expression<Field>{0});
        }
      } else {
        Field upper = 0;
        Field lower = 0;

        const Field range = *state.rows[i].range;

        if (state.rows[i].type == RowSense::LESS_THAN) {
          upper = 0;
          lower = -abs(range);
        } else if (state.rows[i].type == RowSense::GREATER_THAN) {
          upper = abs(range);
          lower = 0;
        } else if (range < 0) {
          upper = 0;
          lower = range;
        } else {
          upper = range;
          lower = 0;
        }

        result.add_constraint(row <= Expression{upper});
        result.add_constraint(row >= Expression{lower});
      }
    }

    return result;
  }
};

}  // namespace mps
