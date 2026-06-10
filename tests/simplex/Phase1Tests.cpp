#include <gtest/gtest.h>

#include "ConstructSparse.h"
#include "linear/BigInteger.h"
#include "linear/simplex/Feasibility.h"
#include "linear/simplex/init/primal/Phase1.h"

// Simple two-variable equality: x1 + x2 = 3, x1,x2 >= 0.
TEST(PrimalPhase1Tests, SimpleFeasible) {
  auto A = sparse<Rational>({{1, 1}});
  Vector<Rational> b = {3};
  Vector<Rational> c = {0, 0};

  Bounds<Rational> bounds(2);
  bounds[0] = {0, std::nullopt};
  bounds[1] = {0, std::nullopt};

  auto result = simplex::primal_phase1(A, b, c, bounds);

  ASSERT_TRUE(result.has_value());
  ASSERT_TRUE(simplex::is_primal_feasible(A, b, c, bounds, *result));
}

// One row has positive adjusted RHS, the other negative.
// The two branches of artificial-variable sign selection are both exercised.
// System: x1 = 2, -x2 = -3 (i.e. x2 = 3), unique feasible point (2, 3).
TEST(PrimalPhase1Tests, NegativeRHSFeasible) {
  auto A = sparse<Rational>({{1, 0}, {0, -1}});
  Vector<Rational> b = {2, -3};
  Vector<Rational> c = {0, 0};

  Bounds<Rational> bounds(2);
  bounds[0] = {0, std::nullopt};
  bounds[1] = {0, std::nullopt};

  auto result = simplex::primal_phase1(A, b, c, bounds);

  ASSERT_TRUE(result.has_value());
  ASSERT_TRUE(simplex::is_primal_feasible(A, b, c, bounds, *result));
}

// Variables have both lower and upper bounds.
// System: x1 + x2 = 5, x1 in [2,4], x2 in [1,3]. Feasible at, e.g., (2,3).
TEST(PrimalPhase1Tests, FeasibleWithUpperBounds) {
  auto A = sparse<Rational>({{1, 1}});
  Vector<Rational> b = {5};
  Vector<Rational> c = {0, 0};

  Bounds<Rational> bounds(2);
  bounds[0] = {2, 4};
  bounds[1] = {1, 3};

  auto result = simplex::primal_phase1(A, b, c, bounds);

  ASSERT_TRUE(result.has_value());
  ASSERT_TRUE(simplex::is_primal_feasible(A, b, c, bounds, *result));
}

// System: x1 = 5, but x1 in [0,4]. No feasible point.
TEST(PrimalPhase1Tests, Infeasible) {
  auto A = sparse<Rational>({{1, 0}, {0, 1}});
  Vector<Rational> b = {5, 3};
  Vector<Rational> c = {0, 0};

  Bounds<Rational> bounds(2);
  bounds[0] = {0, 4};
  bounds[1] = {0, 2};

  auto result = simplex::primal_phase1(A, b, c, bounds);

  ASSERT_FALSE(result.has_value());
  ASSERT_EQ(result.error(), simplex::Phase1Error::INFEASIBLE);
}

// Rows are identical and RHS is consistent ([1,0]x=1 twice).
// The system has solutions (x1=1 works), but the rows are linearly dependent,
// so no basis of full column rank exists for the original matrix.
TEST(PrimalPhase1Tests, LinearlyDependentRowsConsistent) {
  auto A = sparse<Rational>({{1, 0}, {1, 0}});
  Vector<Rational> b = {1, 1};
  Vector<Rational> c = {0, 0};

  Bounds<Rational> bounds(2);
  bounds[0] = {0, std::nullopt};
  bounds[1] = {0, std::nullopt};

  auto result = simplex::primal_phase1(A, b, c, bounds);

  ASSERT_FALSE(result.has_value());
  ASSERT_EQ(result.error(), simplex::Phase1Error::LINEARLY_DEPENDENT_ROWS);
}

// Rows are identical but RHS is inconsistent ([1,0]x=1 and [1,0]x=2).
// The system is infeasible, so INFEASIBLE is detected before LD rows.
TEST(PrimalPhase1Tests, LinearlyDependentRowsInconsistent) {
  auto A = sparse<Rational>({{1, 0}, {1, 0}});
  Vector<Rational> b = {1, 2};
  Vector<Rational> c = {0, 0};

  Bounds<Rational> bounds(2);
  bounds[0] = {0, std::nullopt};
  bounds[1] = {0, std::nullopt};

  auto result = simplex::primal_phase1(A, b, c, bounds);

  ASSERT_FALSE(result.has_value());
  ASSERT_EQ(result.error(), simplex::Phase1Error::INFEASIBLE);
}
