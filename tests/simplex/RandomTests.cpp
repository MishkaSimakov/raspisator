#include <gtest/gtest.h>

#include <random>

#include "Assertions.h"
#include "field/BigInteger.h"
#include "linalg/Matrix.h"
#include "linalg/Random.h"
#include "presolve/passes/TransformToEqualities.h"
#include "simplex/Simplex.h"
#include "simplex/init/dual/ReducedCost.h"
#include "simplex/init/primal/Phase1.h"
#include "support/Highs.h"
#include "support/RandomProblem.h"

TEST(RandomSimplexMethodTests, SimpleRandomMatrixDual) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 10;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine engine(0);

  for (size_t iteration = 0; iteration < kIterations; ++iteration) {
    std::cout << "#" << iteration << std::endl;

    auto problem =
        random_feasible_problem<Rational>(kSize, kElementMagnitude, engine);
    problem = presolve::TransformToEqualities<Rational>().apply(problem);

    problem::StandardLP<Rational> standard_lp(problem);

    auto solver = simplex::Simplex<Rational>();

    solver.set_problem(standard_lp);

    // all variables have all bounds -> this method is guaranteed to find dual
    // feasible point
    auto states = simplex::try_init_dual_by_reduced_cost(standard_lp);

    ASSERT_TRUE(states.has_value());

    auto run_result = solver.dual(*states);

    // check solution
    ASSERT_TRUE(run_result.is_feasible());

    auto finite_solution =
        std::get<simplex::FiniteLPSolution<Rational>>(run_result.solution);

    ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(
        Matrix(problem.matrix), standard_lp.rhs, problem.cost,
        problem.var_bounds, finite_solution));
  }
}

TEST(RandomSimplexMethodTests, SimpleRandomMatrixPrimal) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 10;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine engine(0);

  for (size_t iteration = 0; iteration < kIterations; ++iteration) {
    std::cout << "#" << iteration << std::endl;

    auto problem =
        random_feasible_problem<Rational>(kSize, kElementMagnitude, engine);
    problem = presolve::TransformToEqualities<Rational>().apply(problem);

    problem::StandardLP<Rational> standard_lp(problem);

    auto states = simplex::primal_phase1(standard_lp);

    ASSERT_TRUE(states.has_value());

    auto solver = simplex::Simplex<Rational>();
    solver.set_problem(standard_lp);

    auto run_result = solver.primal(*states);

    // check solution
    ASSERT_TRUE(run_result.is_feasible());

    auto finite_solution =
        std::get<simplex::FiniteLPSolution<Rational>>(run_result.solution);

    ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(
        Matrix(problem.matrix), standard_lp.rhs, problem.cost,
        problem.var_bounds, finite_solution));
  }
}
