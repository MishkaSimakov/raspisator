#include <gtest/gtest.h>

#include <format>
#include <sstream>

#include "linear/problem/VariableType.h"
#include "mps/MPS.h"

using namespace mps;

static MILPProblem<double> parse(std::string_view text) {
  std::stringstream ss{std::string(text)};
  return MPS<double>::read(ss, Format::FREE);
}

static std::string fmt(const Constraint<double>& c) {
  return std::format("{}", c);
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

  ASSERT_EQ(problem.variables.size(), 1u);
  ASSERT_EQ(problem.variables[0].name, "x1");
  ASSERT_EQ(problem.constraints.size(), 0u);
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

  ASSERT_EQ(problem.variables.size(), 2u);
  ASSERT_EQ(problem.constraints.size(), 2u);
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

  ASSERT_EQ(problem.constraints.size(), 1u);
  ASSERT_EQ(fmt(problem.constraints[0]), "x1 + -5 <= 0");
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

  ASSERT_EQ(problem.constraints.size(), 1u);
  ASSERT_EQ(fmt(problem.constraints[0]), "-x1 + 5 <= 0");
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

  ASSERT_EQ(problem.constraints.size(), 1u);
  ASSERT_EQ(fmt(problem.constraints[0]), "x1 + -7 == 0");
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

  ASSERT_EQ(problem.variables[0].type, VariableType::REAL);
  ASSERT_EQ(problem.variables[0].bound, (Bound<double>{0.0, std::nullopt}));
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

  ASSERT_EQ(problem.variables[0].bound, (Bound<double>{0.0, 10.0}));
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

  ASSERT_EQ(problem.variables[0].bound, (Bound<double>{5.0, 5.0}));
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

  ASSERT_EQ(problem.variables[0].bound,
            (Bound<double>{std::nullopt, std::nullopt}));
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

  ASSERT_EQ(problem.variables[0].type, VariableType::INTEGER);
  ASSERT_EQ(problem.variables[0].bound, (Bound<double>{0.0, 1.0}));
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

  ASSERT_EQ(problem.variables[0].type, VariableType::REAL);
  ASSERT_EQ(problem.variables[1].type, VariableType::INTEGER);
  ASSERT_EQ(problem.variables[1].bound, (Bound<double>{0.0, 1.0}));
}

TEST(MPSTests, RangeConstraintProducesTwoConstraints) {
  // A ranged L row with range r gives: rhs - |r| <= expr <= rhs,
  // translated into two constraints in the problem.
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

  ASSERT_EQ(problem.constraints.size(), 2u);
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

  ASSERT_ANY_THROW({ MPS<double>::read(ss, Format::FREE); });
}
