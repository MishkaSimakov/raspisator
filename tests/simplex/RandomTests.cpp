#include <gtest/gtest.h>

#include <random>

#include "Assertions.h"
#include "linalg/Matrix.h"
#include "linalg/Random.h"
#include "linear/BigInteger.h"
#include "linear/simplex/Simplex.h"
#include "linear/simplex/init/dual/ReducedCost.h"
#include "linear/simplex/init/primal/Phase1.h"
#include "presolve/passes/TransformToEqualities.h"
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

    // calculate solution
    Vector<Rational> b(problem.matrix.rows());
    for (size_t i = 0; i < problem.matrix.rows(); ++i) {
      b[i] = *problem.rhs_bounds[i].lower;
    }

    auto solver = simplex::Simplex(problem.matrix, b, problem.cost);

    // all variables have all bounds -> this method is guaranteed to find dual
    // feasible point
    auto states = simplex::try_init_dual_by_reduced_cost(
        problem.matrix, b, problem.cost, Bounds(problem.var_bounds));

    ASSERT_TRUE(states.has_value());

    auto run_result = solver.dual(Bounds(problem.var_bounds), *states);

    // check solution
    ASSERT_TRUE(run_result.is_feasible());

    auto finite_solution =
        std::get<FiniteLPSolution<Rational>>(run_result.solution);

    ASSERT_NO_FATAL_FAILURE(
        validate_simplex_solution(Matrix(problem.matrix), b, problem.cost,
                                  Bounds(problem.var_bounds), finite_solution));
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

    // calculate solution
    Vector<Rational> b(problem.matrix.rows());
    for (size_t i = 0; i < problem.matrix.rows(); ++i) {
      b[i] = *problem.rhs_bounds[i].lower;
    }

    auto solver = simplex::Simplex(problem.matrix, b, problem.cost);

    auto states = simplex::primal_phase1(problem.matrix, b, problem.cost,
                                         Bounds(problem.var_bounds));

    ASSERT_TRUE(states.has_value());

    auto run_result = solver.primal(Bounds(problem.var_bounds), *states);

    // check solution
    ASSERT_TRUE(run_result.is_feasible());

    auto finite_solution =
        std::get<FiniteLPSolution<Rational>>(run_result.solution);

    ASSERT_NO_FATAL_FAILURE(
        validate_simplex_solution(Matrix(problem.matrix), b, problem.cost,
                                  Bounds(problem.var_bounds), finite_solution));
  }
}
