#include <gtest/gtest.h>

#include <random>

#include "../../src/linear/simplex/BoundedSimplexMethod.h"
#include "Assertions.h"
#include "linear/BigInteger.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/Random.h"

TEST(RandomSimplexMethodTests, SimpleRandomMatrix) {
  constexpr size_t kIterations = 1'000;
  constexpr size_t kSize = 10;
  constexpr int kElementMagnitude = 10;

  std::default_random_engine engine(0);
  std::uniform_int_distribution<int> height_distribution(1, kSize);
  std::uniform_int_distribution<int> width_increase_distribution(1, kSize);
  std::uniform_int_distribution<int> elements_distribution(-kElementMagnitude,
                                                           kElementMagnitude);

  std::uniform_int_distribution<int> coin_distribution(0, 1);

  auto elements_generator = [&elements_distribution, &engine] {
    return elements_distribution(engine);
  };

  for (size_t iteration = 0; iteration < kIterations; ++iteration) {
    std::cout << "#" << iteration << std::endl;

    // generate an LP-problem
    size_t n = height_distribution(engine);
    size_t d = n + width_increase_distribution(engine);

    auto A_basic = linalg::random_invertible<Rational>(n, elements_generator);
    auto A_nonbasic = linalg::random<Rational>(n, d - n, elements_generator);

    auto c = linalg::random<Rational>(1, d, elements_generator);

    Matrix<Rational> point(d, 1);

    std::vector<Rational> lower(d);
    std::vector<Rational> upper(d);
    for (size_t i = 0; i < d; ++i) {
      int first = elements_generator();
      int second = elements_generator();

      if (first > second) {
        std::swap(first, second);
      }

      lower[i] = first;
      upper[i] = second;

      point[i, 0] = std::uniform_int_distribution<int>(first, second)(engine);
    }

    auto A = linalg::hstack(A_basic, A_nonbasic);
    auto b = A * point;

    std::vector<size_t> basic_variables(n);
    std::iota(basic_variables.begin(), basic_variables.end(), 0);

    // calculate solution
    auto solver = simplex::BoundedSimplexMethod(CSCMatrix(A), b, c);
    solver.setup_warm_start(basic_variables);
    auto solution = solver.dual(lower, upper);

    // check solution
    auto finite_solution = std::get<FiniteLPSolution<Rational>>(solution);

    ASSERT_NO_FATAL_FAILURE(
        validate_simplex_solution(A, b, c, lower, upper, finite_solution));
  }
}
