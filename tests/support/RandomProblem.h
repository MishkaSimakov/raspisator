#pragma once

#include "linear/BigInteger.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/Random.h"
#include "linear/simplex/Simplex.h"
#include "support/RandomProblem.h"

inline auto random_feasible_problem(size_t size, size_t magnitude,
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
