#pragma once

#include <format>

#include "linalg/CSCMatrix.h"
#include "linalg/lu/LUPA.h"
#include "problem/StandardLP.h"
#include "simplex/SimplexMath.h"

namespace simplex {

// Returns empty string if basis is dual feasible, otherwise returns string
// description of the dual infeasibility reason.
// Note: This function does not check whether a matrix formed by basic columns
// is invertible.
template <typename Field>
std::string get_dual_infeasibility_reason(
    const CSCMatrix<Field>& A, const Vector<Field>& b, const Vector<Field>& c,
    const std::vector<Bound<Field>>& bounds,
    const std::vector<VariableState>& states,
    Field tolerance = FieldTraits<Field>::tolerance) {
  const auto [n, d] = A.shape();

  std::vector<size_t> basic_variables;
  for (size_t i = 0; i < states.size(); ++i) {
    if (states[i] == VariableState::BASIC) {
      basic_variables.push_back(i);
    }
  }

  if (basic_variables.size() != n) {
    return "Wrong basic variables count.";
  }

  linalg::LUPA<Field> lupa(A);
  lupa.set_columns(basic_variables);

  const Vector pi = lupa.solve_linear_transposed(Vector(c[basic_variables]));
  const Vector reduced_costs = c - linalg::transpose(A) * pi;

  for (size_t i = 0; i < states.size(); ++i) {
    if (reduced_costs[i] > tolerance) {
      // infeasible if variable value can be increased
      if (!bounds[i].upper ||
          states[i] == VariableState::AT_LOWER && !bounds[i].is_fixed()) {
        return std::format(
            "Variable {} has reduced cost {} > 0 and can be increased.", i,
            reduced_costs[i]);
      }
    }

    if (reduced_costs[i] < -tolerance) {
      // infeasible if variable value can be decreased
      if (!bounds[i].lower ||
          states[i] == VariableState::AT_UPPER && !bounds[i].is_fixed()) {
        return std::format(
            "Variable {} has reduced cost {} < 0 and can be decreased.", i,
            reduced_costs[i]);
      }
    }
  }

  // the point is dual feasible
  return "";
}

template <typename Field>
bool is_dual_feasible(const problem::StandardLP<Field>& problem,
                      const std::vector<VariableState>& states,
                      Field tolerance = FieldTraits<Field>::tolerance) {
  return get_dual_infeasibility_reason(problem.matrix, problem.rhs,
                                       problem.cost, problem.var_bounds, states,
                                       tolerance)
      .empty();
}

template <typename Field>
std::string get_dual_infeasibility_reason(
    const problem::StandardLP<Field>& problem,
    const std::vector<VariableState>& states,
    Field tolerance = FieldTraits<Field>::tolerance) {
  return get_dual_infeasibility_reason(problem.matrix, problem.rhs,
                                       problem.cost, problem.var_bounds, states,
                                       tolerance);
}

// Returns empty string if basis is primal feasible, otherwise returns string
// description of the primal infeasibility reason.
// Note: This function does not check whether a matrix formed by basic columns
// is invertible.
template <typename Field>
std::string get_primal_infeasibility_reason(
    const CSCMatrix<Field>& A, const Vector<Field>& b, const Vector<Field>& c,
    const std::vector<Bound<Field>>& bounds,
    const std::vector<VariableState>& states,
    Field tolerance = FieldTraits<Field>::tolerance) {
  const auto [n, d] = A.shape();

  std::vector<size_t> basic_variables;
  for (size_t i = 0; i < states.size(); ++i) {
    if (states[i] == VariableState::BASIC) {
      basic_variables.push_back(i);
    }
  }

  if (basic_variables.size() != n) {
    return "Wrong basic variables count.";
  }

  auto rhs = detail::get_adjusted_rhs(A, b, bounds, states);

  linalg::LUPA<Field> lupa(A);
  lupa.set_columns(basic_variables);
  auto basic_point = lupa.solve_linear(rhs);

  for (size_t i = 0; i < n; ++i) {
    if (!bounds[basic_variables[i]].contains(basic_point[i], tolerance)) {
      return std::format("Variable {} is outside of its bound: {} \\notin {}.",
                         basic_variables[i], basic_point[i],
                         bounds[basic_variables[i]]);
    }
  }

  // the point is primal feasible
  return "";
}

template <typename Field>
bool is_primal_feasible(const problem::StandardLP<Field>& problem,
                        const std::vector<VariableState>& states,
                        Field tolerance = FieldTraits<Field>::tolerance) {
  return get_primal_infeasibility_reason(problem.matrix, problem.rhs,
                                         problem.cost, problem.var_bounds,
                                         states, tolerance)
      .empty();
}

template <typename Field>
std::string get_primal_infeasibility_reason(
    const problem::StandardLP<Field>& problem,
    const std::vector<VariableState>& states,
    Field tolerance = FieldTraits<Field>::tolerance) {
  return get_primal_infeasibility_reason(problem.matrix, problem.rhs,
                                         problem.cost, problem.var_bounds,
                                         states, tolerance);
}

}  // namespace simplex
