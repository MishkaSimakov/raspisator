#pragma once

#include <expected>
#include <vector>

#include "linalg/CSCMatrix.h"
#include "simplex/SimplexMath.h"
#include "simplex/Simplex.h"

namespace simplex {

enum class Phase1Error {
  INFEASIBLE,
  LINEARLY_DEPENDENT_ROWS,
};

// Returns primal feasible basis or std::nullopt if problem is infeasible.
// Note: algorithm is taken from
// https://people.orie.cornell.edu/dpw/orie6300/Lectures/lec12.pdf
template <typename Field>
std::expected<std::vector<VariableState>, Phase1Error> primal_phase1(
    const problem::StandardLP<Field>& problem, Config<Field> config = {}) {
  using std::abs;

  const auto [n, old_d] = problem.matrix.shape();
  const size_t new_d = old_d + n;

  std::vector<VariableState> states(new_d);

  auto new_problem = problem;
  new_problem.var_bounds.resize(new_d);
  new_problem.var_names.resize(new_d);

  for (size_t i = 0; i < old_d; ++i) {
    if (problem.var_bounds[i].lower) {
      states[i] = VariableState::AT_LOWER;
    } else if (problem.var_bounds[i].upper) {
      states[i] = VariableState::AT_UPPER;
    } else {
      states[i] = VariableState::NONBASIC_FREE;
    }
  }

  for (size_t i = old_d; i < new_d; ++i) {
    states[i] = VariableState::BASIC;
    new_problem.var_bounds[i] = {0, std::nullopt};
    new_problem.var_names[i] = "phase1_slack" + std::to_string(i);
  }

  // add additional slack variable to each constraint
  const auto rhs = detail::get_adjusted_rhs(problem, states);

  for (size_t i = 0; i < n; ++i) {
    new_problem.matrix.add_column();

    if (rhs[i, 0] > 0) {
      new_problem.matrix.push_to_last_column(i, 1);
    } else {
      new_problem.matrix.push_to_last_column(i, -1);
    }
  }

  new_problem.cost = Vector<Field>::zeros(new_d);

  for (size_t i = old_d; i < new_d; ++i) {
    new_problem.cost[i] = -1;
  }

  // save pivot tolerance, it will be needed later
  const Field pivot_tolerance = config.tolerance.pivot;

  auto helper = Simplex<Field>(std::move(config));
  helper.set_problem(new_problem);

  const auto result = helper.primal(states);

  if (!result.is_feasible()) {
    throw std::runtime_error("Something went wrong in primal implementation.");
  }

  auto solution = std::get<FiniteLPSolution<Field>>(result.solution);

  // Case 1
  if (FieldTraits<Field>::is_strictly_negative(solution.value)) {
    return std::unexpected{Phase1Error::INFEASIBLE};
  }

  // Case 2: try to eliminate artificial variables from basic variables (if
  // there are any) using pivot operation
  for (size_t basic_index = 0; basic_index < n; ++basic_index) {
    const size_t i = helper.get_basic_vars()[basic_index];

    if (i < old_d) {
      continue;
    }

    const auto row = helper.get_tableau_row(basic_index);

    // try to find replacement for i among non-artificial variables
    ArgMaximum<Field> max_pivot;

    for (size_t j = 0; j < old_d; ++j) {
      if (solution.variables[j] == VariableState::BASIC) {
        continue;
      }

      KahanSum<Field> coef;

      for (const auto [index, value] : problem.matrix.get_column(j)) {
        coef.add(row[index] * value);
      }

      max_pivot.record(j, abs(coef.sum()));
    }

    if (!max_pivot.has_value() || max_pivot->max <= pivot_tolerance) {
      // Problem contains linearly dependent rows.
      return std::unexpected{Phase1Error::LINEARLY_DEPENDENT_ROWS};
    }

    solution.variables[max_pivot->index] = VariableState::BASIC;
    helper.change_basis(basic_index, max_pivot->index, VariableState::AT_LOWER);
  }

  solution.variables.resize(old_d);
  return solution.variables;
}

}  // namespace simplex
