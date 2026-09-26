#include <gtest/gtest.h>

#include <sstream>

#include "mps/Format.h"
#include "mps/MPS.h"
#include "problem/MILP.h"

using namespace mps;

static problem::MILP<double> parse(std::string_view text) {
  std::stringstream ss{std::string(text)};
  return read<double>(ss, Format::FREE);
}

TEST(MPSTests, MinimalProblem) {
  const auto problem = parse(
      "NAME minimal\n"
      "ROWS\n"
      " N obj\n"
      "COLUMNS\n"
      "   x1 obj 3\n"
      "RHS\n"
      "   RHS obj 0\n"
      "ENDATA");

  ASSERT_EQ(problem.matrix.shape(), (std::pair{0, 1}));

  ASSERT_EQ(problem.var_names.size(), 1);
  ASSERT_EQ(problem.var_names[0], "x1");
  ASSERT_EQ(problem.row_names.size(), 0);
}

TEST(MPSTests, MultipleVariablesAndConstraints) {
  const auto problem = parse(
      "NAME multi\n"
      "ROWS\n"
      " N obj\n"
      " L c1\n"
      " G c2\n"
      "COLUMNS\n"
      "   x1 obj 1  c1 2\n"
      "   x2 obj 4  c2 5\n"
      "RHS\n"
      "   RHS c1 10  c2 3\n"
      "ENDATA");

  ASSERT_EQ(problem.matrix.shape(), (std::pair{2, 2}));

  ASSERT_EQ(problem.var_names.size(), 2);
  ASSERT_EQ(problem.row_names.size(), 2);
}

// ---------------------------------------------------------------------------
// Constraint type tests
// ---------------------------------------------------------------------------

TEST(MPSTests, LessThanConstraint) {
  const auto problem = parse(
      "NAME lt\n"
      "ROWS\n"
      " N obj\n"
      " L c1\n"
      "COLUMNS\n"
      "   x1 obj 0  c1 1\n"
      "RHS\n"
      "   RHS c1 5\n"
      "ENDATA");

  ASSERT_EQ(problem.matrix.shape(), (std::pair{1, 1}));
  ASSERT_EQ(problem.rhs_bounds[0], (Bound<double>{std::nullopt, 5.0}));
}

TEST(MPSTests, GreaterThanConstraint) {
  const auto problem = parse(
      "NAME gt\n"
      "ROWS\n"
      " N obj\n"
      " G c1\n"
      "COLUMNS\n"
      "   x1 obj 0  c1 1\n"
      "RHS\n"
      "   RHS c1 5\n"
      "ENDATA");

  ASSERT_EQ(problem.matrix.shape(), (std::pair{1, 1}));
  ASSERT_EQ(problem.rhs_bounds[0], (Bound<double>{5.0, std::nullopt}));
}

TEST(MPSTests, EqualityConstraint) {
  const auto problem = parse(
      "NAME eq\n"
      "ROWS\n"
      " N obj\n"
      " E c1\n"
      "COLUMNS\n"
      "   x1 obj 0  c1 1\n"
      "RHS\n"
      "   RHS c1 7\n"
      "ENDATA");

  ASSERT_EQ(problem.matrix.shape(), (std::pair{1, 1}));
  ASSERT_EQ(problem.rhs_bounds[0], (Bound<double>{7.0, 7.0}));
}

TEST(MPSTests, DefaultRealBound) {
  // Without a BOUNDS section a real variable gets [0, +∞).
  const auto problem = parse(
      "NAME defbound\n"
      "ROWS\n"
      " N obj\n"
      "COLUMNS\n"
      "   x1 obj 1\n"
      "RHS\n"
      "   RHS obj 0\n"
      "ENDATA");

  ASSERT_EQ(problem.matrix.shape(), (std::pair{0, 1}));
  ASSERT_FALSE(problem.is_integer[0]);
  ASSERT_EQ(problem.var_bounds[0], (Bound<double>{0.0, std::nullopt}));
}

TEST(MPSTests, UpperBound) {
  const auto problem = parse(
      "NAME upbound\n"
      "ROWS\n"
      " N obj\n"
      "COLUMNS\n"
      "   x1 obj 1\n"
      "RHS\n"
      "   RHS obj 0\n"
      "BOUNDS\n"
      "   UP BND x1 10\n"
      "ENDATA");

  ASSERT_EQ(problem.var_bounds[0], (Bound<double>{0.0, 10.0}));
}

TEST(MPSTests, FixedBound) {
  const auto problem = parse(
      "NAME fxbound\n"
      "ROWS\n"
      " N obj\n"
      "COLUMNS\n"
      "   x1 obj 1\n"
      "RHS\n"
      "   RHS obj 0\n"
      "BOUNDS\n"
      "   FX BND x1 5\n"
      "ENDATA");

  ASSERT_EQ(problem.var_bounds[0], (Bound<double>{5.0, 5.0}));
}

TEST(MPSTests, FreeBound) {
  // FR makes a variable unbounded in both directions.
  const auto problem = parse(
      "NAME frbound\n"
      "ROWS\n"
      " N obj\n"
      "COLUMNS\n"
      "   x1 obj 1\n"
      "RHS\n"
      "   RHS obj 0\n"
      "BOUNDS\n"
      "   FR BND x1\n"
      "ENDATA");

  ASSERT_EQ(problem.var_bounds[0], (Bound<double>{std::nullopt, std::nullopt}));
}

TEST(MPSTests, BinaryVariableBound) {
  // BV sets [0, 1] and marks the variable as integer.
  const auto problem = parse(
      "NAME bvbound\n"
      "ROWS\n"
      " N obj\n"
      "COLUMNS\n"
      "   x1 obj 1\n"
      "RHS\n"
      "   RHS obj 0\n"
      "BOUNDS\n"
      "   BV BND x1\n"
      "ENDATA");

  ASSERT_TRUE(problem.is_integer[0]);
  ASSERT_EQ(problem.var_bounds[0], (Bound<double>{0.0, 1.0}));
}

TEST(MPSTests, IntegerMarker) {
  // Variables inside INTORG/INTEND markers are integer and get [0, 1] default.
  const auto problem = parse(
      "NAME intvar\n"
      "ROWS\n"
      " N obj\n"
      "COLUMNS\n"
      "   x1 obj 1\n"
      "   MARK001 'MARKER' 'INTORG'\n"
      "   x2 obj 2\n"
      "   MARK002 'MARKER' 'INTEND'\n"
      "RHS\n"
      "   RHS obj 0\n"
      "ENDATA");

  ASSERT_FALSE(problem.is_integer[0]);
  ASSERT_TRUE(problem.is_integer[1]);
  ASSERT_EQ(problem.var_bounds[1], (Bound<double>{0.0, 1.0}));
}

TEST(MPSTests, RangeConstraint) {
  // A ranged L row with range r gives a two-sided bound: [rhs - |r|, rhs].
  const auto problem = parse(
      "NAME ranged\n"
      "ROWS\n"
      " N obj\n"
      " L c1\n"
      "COLUMNS\n"
      "   x1 obj 0  c1 1\n"
      "RHS\n"
      "   RHS c1 10\n"
      "RANGES\n"
      "   RNG c1 4\n"
      "ENDATA");

  ASSERT_EQ(problem.matrix.shape(), (std::pair{1, 1}));
  ASSERT_EQ(problem.rhs_bounds[0], (Bound<double>{6.0, 10.0}));
}

TEST(MPSTests, FixedFormatEndToEnd) {
  std::string mps =
      "NAME          FIXTEST\n"
      "ROWS\n"
      " N  NOBJ\n"
      " L  LIM1\n"
      "COLUMNS\n"
      "    X1        NOBJ               1.0   LIM1               2.0\n"
      "RHS\n"
      "    RHS       LIM1              10.0\n"
      "ENDATA";

  std::stringstream ss(mps);
  const auto problem = read<double>(ss, Format::FIXED);

  ASSERT_EQ(problem.matrix.shape(), (std::pair{1, 1}));
  ASSERT_EQ(problem.var_names.size(), 1u);
  ASSERT_EQ(problem.var_names[0], "X1      ");
  ASSERT_EQ(problem.var_bounds[0], (Bound<double>{0.0, std::nullopt}));
}

TEST(MPSTests, MultipleObjectiveRows) {
  // The MPS spec says the first N row is the objective; extra N rows are
  // free rows and are silently skipped during constraint generation.
  const auto problem = parse(
      "NAME multi_n\n"
      "ROWS\n"
      " N obj\n"
      " N extra\n"
      " L c1\n"
      "COLUMNS\n"
      "   x1 obj 3  c1 2\n"
      "   x1 extra 99\n"
      "RHS\n"
      "   RHS c1 10\n"
      "ENDATA");

  ASSERT_EQ(problem.matrix.shape(), (std::pair{2, 1}));
  ASSERT_EQ(problem.var_names.size(), 1u);
}

TEST(MPSTests, DataRowsInObjectSection) {
  // Object section is skipped by the parser. No data rows should be present
  // inside it.
  const auto mps =
      "NAME ranged\n"
      "ROWS\n"
      " N obj\n"
      " L c1\n"
      "COLUMNS\n"
      "   x1 obj 0  c1 1\n"
      "RHS\n"
      "   RHS c1 10\n"
      "RANGES\n"
      "   RNG c1 4\n"
      "OBJECT\n"
      "   RNG2 c1 4\n"
      "ENDATA";

  std::stringstream ss{std::string(mps)};

  ASSERT_ANY_THROW({ read<double>(ss, Format::FREE); });
}

TEST(MPSTests, DuplicateObjectSection) {
  const auto mps =
      "NAME minimal\n"
      "ROWS\n"
      " N obj\n"
      "COLUMNS\n"
      "   x1 obj 3\n"
      "RHS\n"
      "   RHS obj 0\n"
      "OBJECT\n"
      "OBJECT\n"
      "ENDATA";

  std::stringstream ss{std::string(mps)};

  ASSERT_ANY_THROW({ read<double>(ss, Format::FREE); });
}
