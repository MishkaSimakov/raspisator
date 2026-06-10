#pragma once

#include <expected>
#include <vector>

#include "linalg/CSCMatrix.h"
#include "linear/model/LP.h"
#include "linear/simplex/Math.h"
#include "linear/simplex/pricing/primal/MostInfeasible.h"

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
    const CSCMatrix<Field>& A, const Vector<Field>& b, const Vector<Field>& c,
    const Bounds<Field>& bounds) {
  const auto [n, old_d] = A.shape();
  const size_t new_d = old_d + n;

  std::vector<VariableState> states(new_d);
  Bounds<Field> new_bounds(new_d);

  for (size_t i = 0; i < old_d; ++i) {
    if (bounds[i].lower) {
      states[i] = VariableState::AT_LOWER;
    } else if (bounds[i].upper) {
      states[i] = VariableState::AT_UPPER;
    } else {
      states[i] = VariableState::NONBASIC_FREE;
    }

    new_bounds[i] = bounds[i];
  }

  for (size_t i = old_d; i < new_d; ++i) {
    states[i] = VariableState::BASIC;
    new_bounds[i] = {0, std::nullopt};
  }

  // add additional slack variable to each constraint
  const auto rhs = detail::get_adjusted_rhs(A, b, bounds, states);

  auto new_A = A;
  for (size_t i = 0; i < n; ++i) {
    new_A.add_column();

    if (rhs[i, 0] > 0) {
      new_A.push_to_last_column(i, 1);
    } else {
      new_A.push_to_last_column(i, -1);
    }
  }

  auto new_c = Vector<Field>(new_d);

  for (size_t i = old_d; i < new_d; ++i) {
    new_c[i] = -1;
  }

  auto helper = Simplex(new_A, b, new_c);
  const auto result = helper.primal(new_bounds, states);

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
    bool found = false;

    for (size_t j = 0; j < old_d; ++j) {
      if (solution.variables[j] == VariableState::BASIC) {
        continue;
      }

      Field coef = 0;

      for (const auto [index, value] : A.get_column(j)) {
        coef += row[index] * value;
      }

      if (FieldTraits<Field>::is_nonzero(coef)) {
        // change i -> j in basis
        solution.variables[j] = VariableState::BASIC;
        helper.change_basis(basic_index, j, VariableState::AT_LOWER);
        found = true;

        break;
      }
    }

    if (!found) {
      // Problem contains linearly dependent rows.
      return std::unexpected{Phase1Error::LINEARLY_DEPENDENT_ROWS};
    }
  }

  solution.variables.resize(old_d);
  return solution.variables;
}

}  // namespace simplex
