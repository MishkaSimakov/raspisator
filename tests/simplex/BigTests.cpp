#include <gtest/gtest.h>

#include "../../src/linear/simplex/BoundedSimplexMethod.h"
#include "Assertions.h"
#include "linear/BigInteger.h"
#include "linear/model/LP.h"

using Field = double;

#include "simplex_core_dump_0.h"
#include "simplex_core_dump_1.h"

TEST(BigTests, PotentiallyCyclingProblem1) {
  using namespace SimplexDump_0;

  auto solver = simplex::BoundedSimplexMethod(CSCMatrix(A), b, c);
  solver.setup_warm_start(states);
  auto solution = solver.dual(lower, upper);

  auto finite_solution = std::get<FiniteLPSolution<Field>>(solution);

  ASSERT_NO_FATAL_FAILURE(
      validate_simplex_solution(A, b, c, lower, upper, finite_solution));
}

TEST(BigTests, PotentiallyCyclingProblem2) {
  using namespace SimplexDump_1;

  auto solver = simplex::BoundedSimplexMethod(CSCMatrix(A), b, c);
  solver.setup_warm_start(states);
  auto solution = solver.dual(lower, upper);

  auto finite_solution = std::get<FiniteLPSolution<Field>>(solution);

  ASSERT_NO_FATAL_FAILURE(
      validate_simplex_solution(A, b, c, lower, upper, finite_solution));
}
