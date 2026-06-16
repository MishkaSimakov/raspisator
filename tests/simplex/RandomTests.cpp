#include <gtest/gtest.h>

#include <random>

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

    auto result = solver.dual(*states);

    // check solution
    ASSERT_EQ(result.status, simplex::Status::OPTIMAL);

    ASSERT_TRUE(simplex::is_primal_feasible(standard_lp, solver.get_states()));
    ASSERT_TRUE(simplex::is_dual_feasible(standard_lp, solver.get_states()));
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

    auto result = solver.primal(*states);

    // check solution
    ASSERT_EQ(result.status, simplex::Status::OPTIMAL);

    ASSERT_TRUE(simplex::is_primal_feasible(standard_lp, solver.get_states()));
    ASSERT_TRUE(simplex::is_dual_feasible(standard_lp, solver.get_states()));
  }
}
