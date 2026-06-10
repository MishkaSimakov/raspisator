#include <gtest/gtest.h>

#include "ConstructSparse.h"
#include "linalg/CSCMatrix.h"
#include "linear/BigInteger.h"
#include "linear/model/LP.h"
#include "linear/simplex/Simplex.h"
#include "linear/simplex/init/primal/Phase1.h"

TEST(FreeVariablesTests, SmallTest) {
  // \max x
  // x + y = 5
  // x free
  // y \in [1, 2]
  auto A = sparse<Rational>({
      {1, 1},
  });
  Vector<Rational> b = {5};
  Vector<Rational> c = {1, 0};

  Bounds<Rational> bounds(2);
  bounds[0] = {std::nullopt, std::nullopt};
  bounds[1] = {1, 2};

  auto basis = simplex::primal_phase1(A, b, c, bounds);

  ASSERT_TRUE(basis.has_value());

  auto simplex = simplex::Simplex(
      A, b, c,
      {.primal_pricing =
           std::make_unique<simplex::PrimalMostInfeasible<Rational>>()});

  const auto result = simplex.primal(bounds, *basis);

  ASSERT_TRUE(result.is_feasible());
  const auto solution = std::get<FiniteLPSolution<Rational>>(result.solution);

  const Vector<Rational> expected = {4, 1};

  ASSERT_EQ(solution.point, expected);
}

TEST(FreeVariablesTests, AllFree) {
  // \max x
  // x + y = 5
  // x free
  // y free
  auto A = sparse<Rational>({
      {1, 1},
  });
  Vector<Rational> b = {5};
  Vector<Rational> c = {1, 0};

  Bounds<Rational> bounds(2);
  bounds[0] = {std::nullopt, std::nullopt};
  bounds[1] = {std::nullopt, std::nullopt};

  auto simplex = simplex::Simplex(
      A, b, c,
      {.primal_pricing =
           std::make_unique<simplex::PrimalMostInfeasible<Rational>>()});

  auto basis = simplex::primal_phase1(A, b, c, bounds);

  ASSERT_TRUE(basis.has_value());

  const auto result = simplex.primal(bounds, *basis);

  ASSERT_TRUE(std::holds_alternative<Unbounded>(result.solution));
}
