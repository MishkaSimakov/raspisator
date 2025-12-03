#include <gtest/gtest.h>

#include <random>

#include "linear/BigInteger.h"
#include "linear/BoundedSimplexMethod.h"
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
    // generate an LP-problem
    size_t n = height_distribution(engine);
    size_t d = n + width_increase_distribution(engine);

    auto A_basic = linalg::random_invertible<Rational>(n, elements_generator);
    auto A_nonbasic = linalg::random<Rational>(n, d - n, elements_generator);

    auto c = linalg::random<Rational>(1, d, elements_generator);

    Matrix<Rational> point(d, 1);
    std::vector<VariableState> variables(d);

    std::vector<Rational> lower(d);
    std::vector<Rational> upper(d);
    for (size_t i = 0; i < n; ++i) {
      int first = elements_generator();
      int second = elements_generator();

      if (first > second) {
        std::swap(first, second);
      }

      lower[i] = first;
      upper[i] = second;

      variables[i] = VariableState::BASIC;

      point[i, 0] = std::uniform_int_distribution<int>(first, second)(engine);
    }

    for (size_t i = n; i < d; ++i) {
      if (coin_distribution(engine) == 0) {
        variables[i] = VariableState::AT_LOWER;
        point[i, 0] = lower[i];
      } else {
        variables[i] = VariableState::AT_UPPER;
        point[i, 0] = upper[i];
      }
    }

    auto A = linalg::hstack(A_basic, A_nonbasic);
    auto b = A * point;

    // calculate solution
    auto solution = SimplexMethod(CSCMatrix(A), b, c, lower, upper)
                        .solve_from(point, variables);

    // check solution
    auto finite_solution = std::get<FiniteLPSolution<Rational>>(solution);

    ASSERT_EQ(A * finite_solution.point, b);

    for (size_t i = 0; i < d; ++i) {
      ASSERT_TRUE((lower[i] <= finite_solution.point[i, 0]));
      ASSERT_TRUE((finite_solution.point[i, 0] <= upper[i]));

      if (std::ranges::find(finite_solution.basic_variables, i) ==
          finite_solution.basic_variables.end()) {
        ASSERT_TRUE((lower[i] == finite_solution.point[i, 0] ||
                     finite_solution.point[i, 0] == upper[i]));
      }
    }

    // TODO: check reduced costs
  }
}
