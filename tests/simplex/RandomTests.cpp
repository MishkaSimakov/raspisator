#include <gtest/gtest.h>

#include <random>

#include "Assertions.h"
#include "linear/BigInteger.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/Random.h"
#include "linear/simplex/Simplex.h"
#include "support/RandomProblem.h"

TEST(RandomSimplexMethodTests, SimpleRandomMatrixDual) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 10;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine engine(0);

  for (size_t iteration = 0; iteration < kIterations; ++iteration) {
    std::cout << "#" << iteration << std::endl;

    auto [A, b, c, bounds] = random_feasible_problem(kSize, kElementMagnitude, engine);

    // calculate solution
    auto solver = simplex::Simplex(CSCMatrix(A), b, c);

    // all variables have all bounds -> this method is guaranteed to find dual
    // feasible point
    auto states = solver.try_get_dual_feasible(bounds);

    ASSERT_TRUE(states.has_value());

    auto run_result = solver.dual(bounds, *states);

    // check solution
    ASSERT_TRUE(run_result.is_feasible());

    auto finite_solution =
        std::get<FiniteLPSolution<Rational>>(run_result.solution);

    ASSERT_NO_FATAL_FAILURE(
        validate_simplex_solution(A, b, c, bounds, finite_solution));
  }
}

TEST(RandomSimplexMethodTests, SimpleRandomMatrixPrimal) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 10;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine engine(0);

  for (size_t iteration = 0; iteration < kIterations; ++iteration) {
    std::cout << "#" << iteration << std::endl;

    auto [A, b, c, bounds] = random_feasible_problem(kSize, kElementMagnitude, engine);

    // calculate solution
    auto solver = simplex::Simplex(CSCMatrix(A), b, c);

    auto states = solver.get_primal_feasible(bounds);

    ASSERT_TRUE(states.has_value());

    auto run_result = solver.primal(bounds, *states);

    // check solution
    ASSERT_TRUE(run_result.is_feasible());

    auto finite_solution =
        std::get<FiniteLPSolution<Rational>>(run_result.solution);

    ASSERT_NO_FATAL_FAILURE(
        validate_simplex_solution(A, b, c, bounds, finite_solution));
  }
}
