#pragma once

#include <gtest/gtest.h>

#include "presolve/Pass.h"
#include "problem/MILP.h"
#include "support/Highs.h"

inline void assert_pass_correct(const problem::MILP<double>& problem,
                                presolve::Pass<double>& pass) {
  // solve without preprocessing
  auto solution = highs::solve(highs::from_milp(problem));

  // preprocess
  auto new_problem = pass.apply(problem);
  new_problem.validate();

  // solve new problem
  auto new_solution = highs::solve(highs::from_milp(new_problem));

  ASSERT_EQ(solution.status, new_solution.status);

  if (solution.status == HighsModelStatus::kOptimal) {
    // Objective must match
    new_solution = highs::inverse_pass(problem, pass, new_solution);
    ASSERT_NEAR(solution.objective, new_solution.objective, 1e-6);
  }
}

#define ASSERT_PASS_CORRECT(problem, pass) \
  ASSERT_NO_FATAL_FAILURE(assert_pass_correct(problem, pass));
