#include <gtest/gtest.h>

#include <random>

#include "Assertions.h"
#include "linalg/CSCMatrix.h"
#include "linalg/Matrix.h"
#include "linalg/Random.h"
#include "linalg/Stack.h"
#include "linear/BigInteger.h"
#include "linear/simplex/Simplex.h"

using namespace linalg;

auto random_problem(size_t size, size_t magnitude,
                    std::default_random_engine& engine) {
  std::uniform_int_distribution<int> height_distribution(1, size);
  std::uniform_int_distribution<int> width_increase_distribution(1, size);
  std::uniform_int_distribution<int> value_distribution(-magnitude, magnitude);

  // generate an LP-problem
  size_t n = height_distribution(engine);
  size_t d = n + width_increase_distribution(engine);

  auto A_basic =
      random::dense_invertible<Rational>(n, engine, value_distribution);
  auto A_nonbasic =
      random::dense<Rational>(n, d - n, engine, value_distribution);

  Vector c = random::dense<Rational>(d, 1, engine, value_distribution);

  Vector<Rational> point(d);
  Bounds<Rational> bounds(d);

  for (size_t i = 0; i < d; ++i) {
    int first = value_distribution(engine);
    int second = value_distribution(engine);

    if (first > second) {
      std::swap(first, second);
    }

    bounds[i] = Bound<Rational>(first, second);

    point[i, 0] = std::uniform_int_distribution<int>(first, second)(engine);
  }

  Matrix A = hstack(A_basic, A_nonbasic);
  Vector b = A * point;

  return std::make_tuple(std::move(A), std::move(b), std::move(c),
                         std::move(bounds));
}

TEST(RandomSimplexMethodTests, SimpleRandomMatrixDual) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 10;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine engine(0);

  for (size_t iteration = 0; iteration < kIterations; ++iteration) {
    std::cout << "#" << iteration << std::endl;

    auto [A, b, c, bounds] = random_problem(kSize, kElementMagnitude, engine);

    // calculate solution
    auto solver = simplex::Simplex(CSCMatrix<Rational>(A), b, c);

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

    auto [A, b, c, bounds] = random_problem(kSize, kElementMagnitude, engine);

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
