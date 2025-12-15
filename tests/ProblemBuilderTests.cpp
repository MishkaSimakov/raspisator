#include <gtest/gtest.h>

#include "linear/BigInteger.h"
#include "linear/matrix/Matrix.h"
#include "linear/model/LP.h"
#include "linear/problem/MILPProblem.h"
#include "linear/problem/ToMatrices.h"
#include "linear/problem/optimization/RemoveConstantConstraints.h"
#include "linear/problem/optimization/RemoveLinearlyDependentConstraints.h"
#include "linear/problem/optimization/TransformToEqualities.h"
#include "linear/simplex/BoundedSimplexMethod.h"

TEST(ProblemBuilderTests, RemoveConstantConstraints) {
  MILPProblem<Rational> problem;

  auto x = problem.new_variable("x", VariableType::INTEGER, 0, 100);
  auto y = problem.new_variable("y", VariableType::INTEGER, 0, 100);

  problem.add_constraint(x + y <= Expression<Rational>{10});
  problem.add_constraint(x <= Expression<Rational>{10});
  problem.add_constraint(Expression<Rational>{5} <= Expression<Rational>{10});
  problem.add_constraint(Expression<Rational>{10} >= Expression<Rational>{2});

  auto optimized = RemoveConstantConstraints<Rational>().apply(problem);

  ASSERT_EQ(optimized.constraints.size(), 2);
  ASSERT_EQ(std::format("{}", optimized.constraints[0]), "x + y <= 10");
  ASSERT_EQ(std::format("{}", optimized.constraints[1]), "x <= 10");
}

TEST(ProblemBuilderTests, ThrowsWhenFalseConstantConstraint) {
  MILPProblem<Rational> problem;

  problem.add_constraint(Expression<Rational>{5} <= Expression<Rational>{1});

  ASSERT_ANY_THROW(RemoveConstantConstraints<Rational>().apply(problem));
}

TEST(ProblemBuilderTests, RemoveLinearlyDependent) {
  MILPProblem<Rational> problem;

  auto x = problem.new_variable("x", VariableType::INTEGER, 0, 100);
  auto y = problem.new_variable("y", VariableType::INTEGER, 0, 100);

  problem.add_constraint(x + y == Expression<Rational>{10});
  problem.add_constraint(2 * x + 2 * y == Expression<Rational>{20});

  auto optimized =
      RemoveLinearlyDependentConstraints<Rational>().apply(problem);

  ASSERT_EQ(optimized.constraints.size(), 1);
  ASSERT_TRUE(std::format("{}", optimized.constraints[0]) == "x + y == 10" ||
              std::format("{}", optimized.constraints[0]) == "2*x + 2*y == 20");
}

TEST(ProblemBuilderTests, WithSimplexMethod) {
  MILPProblem<Rational> builder;

  // setup constraints
  auto x = builder.new_variable("x", VariableType::REAL, 0, 10);
  auto y = builder.new_variable("y", VariableType::REAL, 0, 10);

  builder.add_constraint(x + y <= Expression<Rational>{10});

  builder.set_objective(x);

  // solve problem
  auto optimized = TransformToEqualities<Rational>().apply(builder);

  auto matrices = to_matrices(optimized);
  auto basic_vars = linalg::get_row_basis(linalg::transposed(matrices.A));

  auto solver = simplex::BoundedSimplexMethod(CSCMatrix(matrices.A), matrices.b,
                                              matrices.c);
  solver.setup_warm_start(basic_vars);
  auto solution = std::get<FiniteLPSolution<Rational>>(
      solver.dual(matrices.lower, matrices.upper).solution);

  // check solution
  Rational x_value = matrices.extract_variable(x, solution.point);
  Rational y_value = matrices.extract_variable(y, solution.point);

  ASSERT_EQ(solution.value, 10);
  ASSERT_EQ(x_value, 10);
  ASSERT_EQ(y_value, 0);
}

// TEST(ProblemBuilderTests, WithBranchAndBound) {
//   ProblemBuilder<Rational> builder;
//
//   // setup constraints
//   auto x = builder.new_variable("x", VariableType::INTEGER);
//   auto y = builder.new_variable("y", VariableType::INTEGER);
//
//   builder.add_constraint(x + y <= Expression{Rational{5} / 3});
//
//   builder.set_objective(x);
//
//   // solve problem
//   auto problem = builder.get_problem();
//
//   auto solution = std::get<FiniteMILPSolution<Rational>>(
//       BranchAndBound<Rational,
//       BoundedSimplexMethod<Rational>>(problem).solve());
//
//   // check solution
//   Rational x_value = builder.extract_variable(solution.point, x);
//   Rational y_value = builder.extract_variable(solution.point, y);
//
//   ASSERT_EQ(solution.value, 1);
//   ASSERT_EQ(x_value, 1);
//   ASSERT_EQ(y_value, 0);
// }
//
// TEST(ProblemBuilderTests, WithBranchAndBoundBigProblem) {
//   ProblemBuilder<Rational> builder;
//
//   // setup constraints
//   auto x1 = builder.new_variable("x1", VariableType::INTEGER);
//   auto x2 = builder.new_variable("x2", VariableType::INTEGER);
//   auto x3 = builder.new_variable("x3", VariableType::INTEGER);
//   auto x4 = builder.new_variable("x4", VariableType::INTEGER);
//
//   Expression budget = 35'000 * x1 + 10'000 * x2 + 25'000 * x3 + 90'000 * x4;
//   Expression space = 4 * x1 + 2 * x2 + 7 * x3 + 3 * x4;
//
//   builder.add_constraint(budget <= Expression<Rational>{120'000});
//   builder.add_constraint(Expression<Rational>{12} >= space);
//   builder.add_constraint(x1 + x2 <= Expression<Rational>{1});
//
//   builder.add_constraint(x1 <= Expression<Rational>{1});
//   builder.add_constraint(x2 <= Expression<Rational>{1});
//   builder.add_constraint(x3 <= Expression<Rational>{1});
//   builder.add_constraint(x4 <= Expression<Rational>{1});
//
//   builder.set_objective(300 * x1 + 90 * x2 + 400 * x3 + 150 * x4);
//
//   // solve problem
//   auto problem = builder.get_problem();
//
//   auto solution = std::get<FiniteMILPSolution<Rational>>(
//       BranchAndBound<Rational,
//       BoundedSimplexMethod<Rational>>(problem).solve());
//
//   // check solution
//   ASSERT_EQ(solution.value, 700);
//   ASSERT_EQ(builder.extract_variable(solution.point, x1), 1);
//   ASSERT_EQ(builder.extract_variable(solution.point, x2), 0);
//   ASSERT_EQ(builder.extract_variable(solution.point, x3), 1);
//   ASSERT_EQ(builder.extract_variable(solution.point, x4), 0);
// }
