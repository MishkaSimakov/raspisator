#pragma once

#include "problem/MILP.h"

template <typename Field>
problem::MILP<Field> feasible_from_matrix(CSCMatrix<Field> matrix) {
  const auto [n, d] = matrix.shape();

  problem::MILP<Field> problem;

  problem.cost = Vector<Field>(d);
  if (n > 0) {
    problem.cost[0] = 1;
  }

  problem.cost_offset = 0;

  // generate feasible point and construct bounds from it
  const auto feasible_point = Vector<Field>::ones(d);

  Vector rhs = matrix * feasible_point;

  problem.var_bounds.resize(d);
  for (size_t i = 0; i < d; ++i) {
    problem.var_bounds[i] =
        Bound<Field>{feasible_point[i] - 1, feasible_point[i] + 1};
  }

  problem.rhs_bounds.resize(n);
  for (size_t i = 0; i < n; ++i) {
    problem.rhs_bounds[i] = Bound<Field>{rhs[i] - 1, rhs[i] + 1};
  }

  problem.implied_var_bounds = problem.var_bounds;

  problem.name = "test";
  problem.cost_name = "cost";

  problem.var_names.resize(d);
  problem.row_names.resize(n);

  problem.proven_infeasible = false;
  problem.proven_unbounded = false;

  problem.matrix = std::move(matrix);

  problem.is_integer.resize(d, false);
  problem.implied_is_integer.resize(d, false);

  return problem;
}
