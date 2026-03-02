#include <gtest/gtest.h>

#include "Assertions.h"
#include "linear/BigInteger.h"
#include "linear/model/LP.h"
#include "linear/simplex/BoundedSimplexMethod.h"

using Field = double;

#include "simplex_core_dump_1.h"
#include "simplex_core_dump_2.h"

TEST(BigTests, PotentiallyCyclingProblem2) {
  using namespace SimplexDump_1;

  auto solver = simplex::BoundedSimplexMethod(CSCMatrix(A), b, c);

  Bounds bounds(lower, upper);
  auto solution = solver.dual(bounds, states).solution;

  auto finite_solution = std::get<FiniteLPSolution<Field>>(solution);

  ASSERT_NO_FATAL_FAILURE(
      validate_simplex_solution(A, b, c, bounds, finite_solution));
}

TEST(BigTests, PotentiallyCyclingProblem3) {
  using namespace SimplexDump_2;

  auto solver = simplex::BoundedSimplexMethod(CSCMatrix(A), b, c);

  Bounds bounds(lower, upper);
  auto solution = solver.dual(bounds, states).solution;

  auto finite_solution = std::get<FiniteLPSolution<Field>>(solution);

  ASSERT_NO_FATAL_FAILURE(
      validate_simplex_solution(A, b, c, bounds, finite_solution));
}
