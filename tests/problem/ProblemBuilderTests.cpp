#include <gtest/gtest.h>

#include <sstream>

#include "field/BigInteger.h"
#include "problem/builder/Builder.h"
#include "support/Highs.h"

TEST(ProblemBuilderTests, RemoveConstantConstraints) {
  problem::Builder<Rational> builder;

  // max x
  // x + y <= 10
  // x <= 10
  auto x = builder.new_variable("x", VariableType::INTEGER, 0, 100);
  auto y = builder.new_variable("y", VariableType::INTEGER, 0, 100);

  builder.set_objective(x);

  builder.add_constraint(x + y <= Expression<Rational>{10});
  builder.add_constraint(x <= Expression<Rational>{10});
  builder.add_constraint(Expression<Rational>{5} <= Expression<Rational>{10});
  builder.add_constraint(Expression<Rational>{10} >= Expression<Rational>{2});

  problem::MILP<Rational> problem(builder);

  problem.validate();

  ASSERT_EQ(problem.matrix.rows(), 4);
  ASSERT_EQ(problem.matrix.cols(), 2);

  const auto solution = highs::solve(highs::from_milp(problem));

  ASSERT_EQ(solution.status, HighsModelStatus::kOptimal);
  ASSERT_DOUBLE_EQ(solution.objective, 10);

  ASSERT_DOUBLE_EQ(solution.x[0], 10);
  ASSERT_DOUBLE_EQ(solution.x[1], 0);
}

// Test that operator<< builds successfully
TEST(ProblemBuilderTests, PrintingBuilds) {
  problem::Builder<Rational> builder;

  std::stringstream ss;

  ss << builder;
}
