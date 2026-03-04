#include <gtest/gtest.h>

#include <random>

#include "Assertions.h"
#include "linear/BigInteger.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/Random.h"
#include "linear/simplex/Simplex.h"

auto random_problem(size_t size, size_t magnitude,
                    std::default_random_engine engine) {
  std::uniform_int_distribution<int> height_distribution(1, size);
  std::uniform_int_distribution<int> width_increase_distribution(1, size);
  std::uniform_int_distribution<int> elements_distribution(-magnitude,
                                                           magnitude);

  auto elements_generator = [&elements_distribution, &engine] {
    return elements_distribution(engine);
  };

  // generate an LP-problem
  size_t n = height_distribution(engine);
  size_t d = n + width_increase_distribution(engine);

  auto A_basic = linalg::random_invertible<Rational>(n, elements_generator);
  auto A_nonbasic = linalg::random<Rational>(n, d - n, elements_generator);

  auto c = linalg::random<Rational>(1, d, elements_generator);

  Matrix<Rational> point(d, 1);

  Bounds<Rational> bounds(d);

  for (size_t i = 0; i < d; ++i) {
    int first = elements_generator();
    int second = elements_generator();

    if (first > second) {
      std::swap(first, second);
    }

    bounds[i] = Bound<Rational>(first, second);

    point[i, 0] = std::uniform_int_distribution<int>(first, second)(engine);
  }

  auto A = linalg::hstack(A_basic, A_nonbasic);
  auto b = A * point;

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

    auto [A, b, c, bounds] = random_problem(kSize, kElementMagnitude, engine);

    // calculate solution
    auto solver = simplex::Simplex(CSCMatrix(A), b, c);

    auto states = solver.try_get_primal_feasible(bounds);

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
