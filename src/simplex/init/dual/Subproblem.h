#pragma once

#include <expected>
#include <vector>

#include "ReducedCost.h"
#include "linalg/Linalg.h"
#include "simplex/Simplex.h"
#include "simplex/SimplexMath.h"

// subproblem approach from 4.2.1 of "The Dual Simplex Method, Techniques for a
// fast and stable implementation"
namespace simplex {

struct SubproblemPhase1Result {
  std::vector<VariableState> states;

  // Phase 1 uses simplex to find primal feasible basis. This is internal
  // simplex method iterations count.
  size_t iterations_count;
};

enum class SubproblemPhase1Error {
  DUAL_INFEASIBLE,
  REACHED_ITERATIONS_LIMIT,
};

template <typename Field>
std::expected<SubproblemPhase1Result, SubproblemPhase1Error>
subproblem_dual_phase1(const problem::StandardLP<Field>& problem,
                       Config<Field> config = {}) {
  using std::abs;

  // save tolerances, they will be needed later
  const auto tolerance = config.tolerance;

  const auto [n, d] = problem.matrix.shape();

  auto new_problem = problem;

  for (size_t i = 0; i < d; ++i) {
    const auto bound = problem.var_bounds[i];

    if (bound.lower && bound.upper) {
      // boxed variable
      new_problem.var_bounds[i] = {0, 0};
    } else if (!bound.lower && bound.upper) {
      // J_u
      new_problem.var_bounds[i] = {-1, 0};
    } else if (bound.lower && !bound.upper) {
      // J_l
      new_problem.var_bounds[i] = {0, 1};
    } else {
      // J_f
      new_problem.var_bounds[i] = {-1, 1};
    }
  }

  new_problem.rhs = Vector<Field>::zeros(n);

  // all variables in new_problem are boxed, so this method can't fail
  auto init_states = try_init_dual_by_reduced_cost(new_problem);
  assert(init_states.has_value());

  // solve phase 1 problem
  auto helper = Simplex<Field>(std::move(config));
  helper.set_problem(new_problem);

  const auto result = helper.dual(*init_states);

  if (result.status == Status::ITERATIONS_LIMIT) {
    return std::unexpected{SubproblemPhase1Error::REACHED_ITERATIONS_LIMIT};
  }

  if (result.status != Status::OPTIMAL) {
    throw std::runtime_error("Something went wrong in dual implementation.");
  }

  // Case 1, problem is dual infeasible
  if (*result.objective > tolerance.feasibility) {
    return std::unexpected{SubproblemPhase1Error::DUAL_INFEASIBLE};
  }

  return SubproblemPhase1Result{
      .states = helper.get_states(),
      .iterations_count = result.iterations_count,
  };
}

}  // namespace simplex
