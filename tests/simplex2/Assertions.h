#pragma once

#include <gtest/gtest.h>

#include "linalg/Linalg.h"
#include "linalg/lu/FullPivotingLU.h"
#include "linalg/lu/Solve.h"
#include "problem/StandardLP.h"
#include "simplex/Feasibility.h"
#include "simplex/Result.h"
#include "simplex/Simplex.h"

// Validates that simplex position is optimal, feasible, and consistent with
// problem. Checks: Ax = b, bounds, nonbasic variable positions, optimality
// conditions.
template <typename Field>
void validate_simplex_solution(const problem::StandardLP<Field>& problem,
                               const simplex::Simplex<Field>& simplex) {
  const auto [n, d] = problem.matrix.shape();

  const Vector point = simplex.get_point();

  // Residual check: A*x == b
  {
    const Vector<Field> residue = problem.matrix * point - problem.rhs;
    for (size_t i = 0; i < n; ++i) {
      ASSERT_TRUE(!FieldTraits<Field>::is_nonzero(residue[i]))
          << "Constraint " << i << " violated: residue = " << residue[i];
    }
  }

  // Bound feasibility and nonbasic position
  for (size_t i = 0; i < d; ++i) {
    const auto& bound = problem.var_bounds[i];

    if (bound.lower) {
      ASSERT_TRUE(
          !FieldTraits<Field>::is_strictly_positive(*bound.lower - point[i]))
          << "Lower bound violated for variable " << i;
    }
    if (bound.upper) {
      ASSERT_TRUE(
          !FieldTraits<Field>::is_strictly_positive(point[i] - *bound.upper))
          << "Upper bound violated for variable " << i;
    }

    if (simplex.get_states()[i] == simplex::VariableState::AT_LOWER) {
      ASSERT_TRUE(bound.lower &&
                  !FieldTraits<Field>::is_nonzero(*bound.lower - point[i]))
          << "Variable " << i << " is AT_LOWER but not at lower bound";
    } else if (simplex.get_states()[i] == simplex::VariableState::AT_UPPER) {
      ASSERT_TRUE(bound.upper &&
                  !FieldTraits<Field>::is_nonzero(*bound.upper - point[i]))
          << "Variable " << i << " is AT_UPPER but not at upper bound";
    }
  }

  // Optimality: reduced costs must have correct sign for nonbasic variables
  const auto basic_vars = simplex.get_basic_vars();

  auto [P, Q, ls, us] =
      linalg::FullPivotingLU<Field>(n).get(problem.matrix, basic_vars);

  const Vector<Field> pi = linalg::solve_linear_transposed(
      Vector<Field>(problem.cost[basic_vars]), P, Q, ls, us);

  const Vector<Field> reduced_cost =
      problem.cost - linalg::transpose(problem.matrix) * pi;

  for (size_t i = 0; i < d; ++i) {
    if (simplex.get_states()[i] == simplex::VariableState::BASIC) {
      continue;
    }

    ASSERT_TRUE((simplex.get_states()[i] == simplex::VariableState::AT_LOWER &&
                 !FieldTraits<Field>::is_strictly_positive(reduced_cost[i])) ||
                (simplex.get_states()[i] == simplex::VariableState::AT_UPPER &&
                 !FieldTraits<Field>::is_strictly_negative(reduced_cost[i])))
        << "Optimality condition violated for variable " << i
        << " (state=" << static_cast<int>(simplex.get_states()[i])
        << ", reduced_cost=" << reduced_cost[i] << ")";
  }
}

#define ASSERT_PRIMAL_FEASIBLE(problem, states) \
  ASSERT_EQ(simplex::get_primal_infeasibility_reason(problem, states), "");

#define ASSERT_DUAL_FEASIBLE(problem, states) \
  ASSERT_EQ(simplex::get_dual_infeasibility_reason(problem, states), "");
