#include <gtest/gtest.h>

#include "../src/linear/bb/BranchAndBound.h"
#include "linear/BigInteger.h"
#include "linear/SimplexMethod.h"
#include "linear/builder/ProblemBuilder.h"
#include "linear/matrix/Matrix.h"

TEST(ProblemBuilderTests, WithSimplexMethod) {
  ProblemBuilder<Rational> builder;

  // setup constraints
  auto x = builder.new_variable("x", VariableType::REAL);
  auto y = builder.new_variable("y", VariableType::REAL);

  builder.add_constraint(x + y <= Expression<Rational>{10});

  builder.set_objective(x);

  // solve problem
  auto problem = builder.get_problem();

  auto solution = std::get<FiniteLPSolution<Rational>>(
      SimplexMethod(problem.A, problem.b, problem.c).solve());

  // check solution
  Rational x_value = builder.extract_variable(solution.point, x);
  Rational y_value = builder.extract_variable(solution.point, y);

  ASSERT_EQ(solution.value, 10);
  ASSERT_EQ(x_value, 10);
  ASSERT_EQ(y_value, 0);
}

TEST(ProblemBuilderTests, WithBranchAndBound) {
  ProblemBuilder<Rational> builder;

  // setup constraints
  auto x = builder.new_variable("x", VariableType::INTEGER);
  auto y = builder.new_variable("y", VariableType::INTEGER);

  builder.add_constraint(x + y <= Expression{Rational{5} / 3});

  builder.set_objective(x);

  // solve problem
  auto problem = builder.get_problem();

  auto solution = std::get<FiniteMILPSolution<Rational>>(
      BranchAndBound<Rational, SimplexMethod<Rational>>(problem).solve());

  // check solution
  Rational x_value = builder.extract_variable(solution.point, x);
  Rational y_value = builder.extract_variable(solution.point, y);

  ASSERT_EQ(solution.value, 1);
  ASSERT_EQ(x_value, 1);
  ASSERT_EQ(y_value, 0);
}

TEST(ProblemBuilderTests, WithBranchAndBoundBigProblem) {
  ProblemBuilder<Rational> builder;

  // setup constraints
  auto x1 = builder.new_variable("x1", VariableType::INTEGER);
  auto x2 = builder.new_variable("x2", VariableType::INTEGER);
  auto x3 = builder.new_variable("x3", VariableType::INTEGER);
  auto x4 = builder.new_variable("x4", VariableType::INTEGER);

  Expression budget = 35'000 * x1 + 10'000 * x2 + 25'000 * x3 + 90'000 * x4;
  Expression space = 4 * x1 + 2 * x2 + 7 * x3 + 3 * x4;

  builder.add_constraint(budget <= Expression<Rational>{120'000});
  builder.add_constraint(Expression<Rational>{12} >= space);
  builder.add_constraint(x1 + x2 <= Expression<Rational>{1});

  builder.add_constraint(x1 <= Expression<Rational>{1});
  builder.add_constraint(x2 <= Expression<Rational>{1});
  builder.add_constraint(x3 <= Expression<Rational>{1});
  builder.add_constraint(x4 <= Expression<Rational>{1});

  builder.set_objective(300 * x1 + 90 * x2 + 400 * x3 + 150 * x4);

  // solve problem
  auto problem = builder.get_problem();

  auto solution = std::get<FiniteMILPSolution<Rational>>(
      BranchAndBound<Rational, SimplexMethod<Rational>>(problem).solve());

  // check solution
  ASSERT_EQ(solution.value, 700);
  ASSERT_EQ(builder.extract_variable(solution.point, x1), 1);
  ASSERT_EQ(builder.extract_variable(solution.point, x2), 0);
  ASSERT_EQ(builder.extract_variable(solution.point, x3), 1);
  ASSERT_EQ(builder.extract_variable(solution.point, x4), 0);
}
