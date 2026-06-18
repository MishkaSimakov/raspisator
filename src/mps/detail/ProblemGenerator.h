#pragma once

#include <cmath>
#include <iostream>

#include "Types.h"
#include "problem/MILP.h"

namespace mps::detail {

template <typename Field>
class ProblemGenerator {
  static Bound<Field> get_var_bound(const Variable<Field>& variable) {
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

  static Bound<Field> get_rhs_bound(const Row<Field>& row) {
    using std::abs;

    if (row.type == RowSense::FREE) {
      return Bound<Field>{std::nullopt, std::nullopt};
    }

    const Field rhs = row.rhs.value_or(0);

    if (!row.range) {
      switch (row.type) {
        case RowSense::LESS_THAN:
          return Bound<Field>{std::nullopt, rhs};
        case RowSense::GREATER_THAN:
          return Bound<Field>{rhs, std::nullopt};
        case RowSense::EQUAL:
          return Bound<Field>{rhs, rhs};
        default:
          throw std::runtime_error("Unknown row type.");
      }
    }

    const Field range = *row.range;

    switch (row.type) {
      case RowSense::LESS_THAN:
        return Bound<Field>{rhs - abs(range), rhs};
      case RowSense::GREATER_THAN:
        return Bound<Field>{rhs, rhs + abs(range)};
      case RowSense::EQUAL:
        return range < 0 ? Bound<Field>{rhs + range, rhs}
                         : Bound<Field>{rhs, rhs + range};
      default:
        throw std::runtime_error("Unknown row type.");
    }
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

 public:
  static problem::MILP<Field> generate(const MPSParsingState<Field>& state) {
    using std::abs;

    problem::MILP<Field> result;

    result.name = state.problem_name;

    // create variables
    result.is_integer.resize(state.cols.size());
    result.var_bounds.resize(state.cols.size());
    result.var_names.resize(state.cols.size());

    for (size_t i = 0; i < state.cols.size(); ++i) {
      result.var_names[i] = state.cols[i].name;
      result.var_bounds[i] = get_var_bound(state.cols[i]);
      result.is_integer[i] = state.cols[i].is_integer;
    }

    // fill in cost information
    const size_t objective_row_index = get_objective_row(state);

    result.cost_name = state.rows[objective_row_index].name;

    double cost_multiplier = 1;
    if (state.objective == ObjectiveType::MINIMIZE) {
      // MPS objective is MINIMIZE, but problem::MILP assumes maximization
      // problem so we need to negate objective coefficients.

      cost_multiplier = -1;
    }

    result.cost.resize(state.cols.size());
    for (size_t i = 0; i < state.cols.size(); ++i) {
      auto itr = state.cols[i].values.find(objective_row_index);

      if (itr != state.cols[i].values.end()) {
        result.cost[i] = itr->second * cost_multiplier;
      }
    }

    // fill in constraints matrix
    // exclude objective row from constraints
    result.matrix.resize(state.rows.size() - 1, 0);

    for (size_t col = 0; col < state.cols.size(); ++col) {
      result.matrix.add_column();

      for (const auto& [row, value] : state.cols[col].values) {
        if (row < objective_row_index) {
          result.matrix.push_to_last_column(row, value);
        } else if (row > objective_row_index) {
          result.matrix.push_to_last_column(row - 1, value);
        }
      }
    }

    // constraints names
    result.row_names.resize(state.rows.size() - 1);

    for (size_t row = 0; row < state.rows.size(); ++row) {
      if (row < objective_row_index) {
        result.row_names[row] = state.rows[row].name;
      } else if (row > objective_row_index) {
        result.row_names[row - 1] = state.rows[row].name;
      }
    }

    // fill in rhs bounds
    result.rhs_bounds.resize(state.rows.size() - 1);

    for (size_t row = 0; row < state.rows.size(); ++row) {
      if (row < objective_row_index) {
        result.rhs_bounds[row] = get_rhs_bound(state.rows[row]);
      } else if (row > objective_row_index) {
        result.rhs_bounds[row - 1] = get_rhs_bound(state.rows[row]);
      }
    }

    // copy implied bounds and integrality
    result.implied_var_bounds = result.var_bounds;
    result.implied_is_integer = result.is_integer;

    return result;
  }
};

}  // namespace mps::detail
