#pragma once

#include <expected>
#include <vector>

#include "linalg/Linalg.h"
#include "simplex/Simplex.h"
#include "simplex/SimplexMath.h"

namespace simplex {

struct Phase1Result {
  std::vector<VariableState> states;
  std::vector<size_t> redundant_rows;

  // Phase 1 uses simplex to find primal feasible basis. This is internal
  // simplex method iterations count.
  size_t iterations_count;
};

enum class Phase1Error {
  INFEASIBLE,
  LINEARLY_DEPENDENT_ROWS,
  REACHED_ITERATIONS_LIMIT,
};

// Note: algorithm is taken from:
// https://people.orie.cornell.edu/dpw/orie6300/Lectures/lec12.pdf
template <typename Field>
std::expected<Phase1Result, Phase1Error> primal_phase1(
    const problem::StandardLP<Field>& problem, Config<Field> config = {}) {
  using std::abs;

  // save tolerances, they will be needed later
  const auto tolerance = config.tolerance;

  const auto [n, d] = problem.matrix.shape();

  std::vector<VariableState> states(d);

  auto new_problem = problem;

  // slack variable for each row (if row has any)
  std::vector<std::optional<size_t>> slacks(n);

  for (size_t col = 0; col < d; ++col) {
    // if variable has singleton column, tentatively make it basic
    if (problem.matrix.get_column(col).size() == 1) {
      const auto [row, coef] = problem.matrix.get_column(col).front();

      if (abs(coef) > tolerance.pivot && !slacks[row]) {
        slacks[row] = col;
        states[col] = VariableState::BASIC;

        continue;
      }
    }

    if (problem.var_bounds[col].lower) {
      states[col] = VariableState::AT_LOWER;
    } else if (problem.var_bounds[col].upper) {
      states[col] = VariableState::AT_UPPER;
    } else {
      states[col] = VariableState::NONBASIC_FREE;
    }
  }

  // add additional slack variable to each constraint that needs it
  auto rhs = detail::get_adjusted_rhs(problem, states);

  for (size_t row = 0; row < n; ++row) {
    if (slacks[row]) {
      const size_t col = *slacks[row];
      const Field coef = problem.matrix.get_column(col).front().second;

      // if constraint is feasible, then no need to add new slacks
      if (problem.var_bounds[col].contains(rhs[row] / coef,
                                           tolerance.feasibility)) {
        continue;
      }

      // otherwise set existing slack to one of its bounds
      if (problem.var_bounds[col].lower) {
        states[col] = VariableState::AT_LOWER;
        rhs[row] -= coef * *problem.var_bounds[col].lower;
      } else if (problem.var_bounds[col].upper) {
        states[col] = VariableState::AT_UPPER;
        rhs[row] -= coef * *problem.var_bounds[col].upper;
      } else {
        states[col] = VariableState::NONBASIC_FREE;
      }
    }

    // add a new slack variable to account for constraint infeasibility
    new_problem.var_names.push_back(std::format("phase1_slack_{}", row));
    new_problem.matrix.add_column();

    if (rhs[row] > 0) {
      new_problem.matrix.push_to_last_column(row, 1);
    } else {
      new_problem.matrix.push_to_last_column(row, -1);
    }
  }

  // add information about new slack variables
  const size_t new_d = new_problem.matrix.cols();
  new_problem.cost = Vector<Field>::zeros(new_d);

  for (size_t i = d; i < new_d; ++i) {
    new_problem.cost[i] = -1;
  }

  new_problem.var_bounds.resize(new_d, Bound<Field>{0, std::nullopt});
  states.resize(new_d, VariableState::BASIC);

  // solve phase 1 problem
  auto helper = Simplex<Field>(std::move(config));
  helper.set_problem(new_problem);

  const auto result = helper.primal(states);

  if (result.status == Status::ITERATIONS_LIMIT) {
    return std::unexpected{Phase1Error::REACHED_ITERATIONS_LIMIT};
  }

  if (result.status != Status::OPTIMAL) {
    throw std::runtime_error("Something went wrong in primal implementation.");
  }

  // Case 1
  if (*result.objective < -tolerance.feasibility) {
    return std::unexpected{Phase1Error::INFEASIBLE};
  }

  std::vector<size_t> redundant_rows;

  // Case 2: try to eliminate artificial variables from basic variables (if
  // there are any) using pivot operation
  for (size_t basic_index = 0; basic_index < n; ++basic_index) {
    const size_t i = helper.get_basic_vars()[basic_index];

    if (i < d) {
      continue;
    }

    const auto row = helper.get_tableau_row(basic_index);

    // try to find replacement for i among non-artificial variables
    ArgMaximum<Field> max_pivot;

    for (size_t j = 0; j < d; ++j) {
      if (helper.get_states()[j] == VariableState::BASIC) {
        continue;
      }

      KahanSum<Field> coef;

      for (const auto [index, value] : problem.matrix.get_column(j)) {
        coef.add(row[index] * value);
      }

      max_pivot.record(j, abs(coef.sum()));
    }

    if (!max_pivot.has_value() || max_pivot->max <= tolerance.pivot) {
      // Row associated with the current slack variable is linearly dependent.
      redundant_rows.push_back(new_problem.matrix.get_column(i).front().first);
      continue;
    }

    helper.change_basis(basic_index, max_pivot->index, VariableState::AT_LOWER);
  }

  states = helper.get_states();
  states.resize(d);

  return Phase1Result{
      .states = std::move(states),
      .redundant_rows = std::move(redundant_rows),
      .iterations_count = result.iterations_count,
  };
}

}  // namespace simplex
