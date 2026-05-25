#pragma once

#include "linear/BigInteger.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/Random.h"
#include "linear/simplex/Simplex.h"
#include "problem/MILP.h"
#include "support/RandomProblem.h"

template <typename Field>
problem::MILP<Field> random_feasible_problem(
    size_t size, int magnitude, std::default_random_engine engine) {
  std::uniform_int_distribution<size_t> height_distribution(1, size);
  std::uniform_int_distribution<size_t> width_increase_distribution(1, size);
  std::uniform_int_distribution<int> elements_distribution(-magnitude,
                                                           magnitude);

  auto elements_generator = [&elements_distribution, &engine] -> int {
    return elements_distribution(engine);
  };

  // generate an LP-problem
  size_t n = height_distribution(engine);
  size_t d = n + width_increase_distribution(engine);

  auto A_basic = linalg::random_invertible<Field>(n, elements_generator);
  auto A_nonbasic = linalg::random<Field>(n, d - n, elements_generator);

  std::vector<Field> c(d);
  for (size_t i = 0; i < d; ++i) {
    c[i] = elements_generator();
  }

  Matrix<Field> point(d, 1);

  std::vector<Bound<Field>> bounds(d);

  for (size_t i = 0; i < d; ++i) {
    int first = elements_generator();
    int second = elements_generator();

    if (first > second) {
      std::swap(first, second);
    }

    bounds[i] = Bound<Field>(first, second);

    point[i, 0] = std::uniform_int_distribution<int>(first, second)(engine);
  }

  auto A = linalg::hstack(A_basic, A_nonbasic);
  auto b = A * point;

  problem::MILP<Field> result;

  result.matrix = CSCMatrix(A);
  result.cost = c;

  result.var_bounds = bounds;

  result.rhs_bounds.resize(n);
  std::uniform_int_distribution<int> bound_range(0, 10);
  for (size_t i = 0; i < n; ++i) {
    result.rhs_bounds[i] = Bound<Field>{b[i, 0] - bound_range(engine),
                                        b[i, 0] + bound_range(engine)};
  }

  result.implied_var_bounds = result.var_bounds;

  result.var_names.resize(d);
  result.row_names.resize(n);

  result.is_integer.resize(d, false);
  result.implied_is_integer = result.is_integer;

  return result;
}
