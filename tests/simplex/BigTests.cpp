#include <gtest/gtest.h>

#include "Assertions.h"
#include "linear/BigInteger.h"
#include "linear/model/LP.h"
#include "linear/simplex/BoundedSimplexMethod.h"

using Field = double;

#include "../../playground/simplex_core_dump_0.h"
#include "../../playground/simplex_core_dump_5.h"
#include "simplex_core_dump_1.h"
#include "simplex_core_dump_2.h"

TEST(BigTests, PotentiallyCyclingProblem1) {
  using namespace SimplexDump_0;

  auto solver = simplex::BoundedSimplexMethod(CSCMatrix(A), b, c);
  auto solution = solver.dual(lower, upper, states).solution;

  auto finite_solution = std::get<FiniteLPSolution<Field>>(solution);

  ASSERT_NO_FATAL_FAILURE(
      validate_simplex_solution(A, b, c, lower, upper, finite_solution));
}

TEST(BigTests, PotentiallyCyclingProblem2) {
  using namespace SimplexDump_1;

  auto solver = simplex::BoundedSimplexMethod(CSCMatrix(A), b, c);
  auto solution = solver.dual(lower, upper, states).solution;

  auto finite_solution = std::get<FiniteLPSolution<Field>>(solution);

  ASSERT_NO_FATAL_FAILURE(
      validate_simplex_solution(A, b, c, lower, upper, finite_solution));
}

TEST(BigTests, PotentiallyCyclingProblem3) {
  using namespace SimplexDump_2;

  auto solver = simplex::BoundedSimplexMethod(CSCMatrix(A), b, c);
  auto solution = solver.dual(lower, upper, states).solution;

  auto finite_solution = std::get<FiniteLPSolution<Field>>(solution);

  ASSERT_NO_FATAL_FAILURE(
      validate_simplex_solution(A, b, c, lower, upper, finite_solution));
}

TEST(BigTests, SingularBasicMatrix) {
  using namespace SimplexDump_5;

  // auto solver = simplex::BoundedSimplexMethod(CSCMatrix(A), b, c);
  // solver.setup_warm_start(states);
  // auto solution = solver.dual(lower, upper).solution;
  //
  // auto finite_solution = std::get<FiniteLPSolution<Field>>(solution);
  //
  // ASSERT_NO_FATAL_FAILURE(
  //     validate_simplex_solution(A, b, c, lower, upper, finite_solution));
}
