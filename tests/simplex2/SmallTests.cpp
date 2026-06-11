#include <gtest/gtest.h>

#include <variant>

#include "Assertions.h"
#include "ConstructSparse.h"
#include "linalg/Matrix.h"
#include "linear/BigInteger.h"
#include "linear/model/LP.h"
#include "linear/simplex/Feasibility.h"
#include "linear/simplex/Simplex.h"
#include "linear/simplex/init/dual/ReducedCost.h"
#include "problem/StandardMILP.h"

using linalg::CSCMatrix, linalg::Matrix, linalg::Vector;

// Convenience: build a StandardMILP<Rational> from dense matrix representation
static problem::StandardMILP<Rational> make_problem(
    linalg::CSCMatrix<Rational> matrix, Vector<Rational> rhs,
    Vector<Rational> cost, std::vector<Bound<Rational>> bounds) {
  const auto [n, d] = matrix.shape();
  problem::StandardMILP<Rational> p;
  p.matrix = std::move(matrix);
  p.rhs = std::move(rhs);
  p.cost = std::move(cost);
  p.var_bounds = std::move(bounds);
  p.var_names.resize(d);
  p.row_names.resize(n);
  p.is_integer.assign(d, false);
  return p;
}

// Convenience: get dual initial states for a problem
static std::optional<std::vector<VariableState>> dual_init(
    const problem::StandardMILP<Rational>& p) {
  return simplex::try_init_dual_by_reduced_cost(p.matrix, p.rhs, p.cost,
                                                p.var_bounds);
}

// ---- Adapted from tests/simplex/SmallTests.cpp ----

TEST(Simplex2SmallTests, StartingInSolution) {
  auto p = make_problem(sparse<Rational>({{1, -1, 1, 0}, {2, 1, 0, 1}}),
                        Vector<Rational>{1, 3}, Vector<Rational>{2, 1, 1, -1},
                        {{Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}}});

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto states = simplex::try_init_dual_by_reduced_cost(
      p.matrix, p.rhs, p.cost, p.var_bounds, std::vector<size_t>{1, 2});
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);
  auto& sol = std::get<FiniteLPSolution<Rational>>(result.solution);

  EXPECT_EQ(sol.point, (Vector<Rational>{0, 3, 4, 0}));
  EXPECT_EQ(sol.value, Rational{7});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, sol));
}

TEST(Simplex2SmallTests, Simple1) {
  auto p = make_problem(sparse<Rational>({{1, -1, 1, 0}, {2, 1, 0, 1}}),
                        Vector<Rational>{1, 3}, Vector<Rational>{2, 1, 1, -1},
                        {{Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}}});

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);
  auto& sol = std::get<FiniteLPSolution<Rational>>(result.solution);

  EXPECT_EQ(sol.point, (Vector<Rational>{0, 3, 4, 0}));
  EXPECT_EQ(sol.value, Rational{7});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, sol));
}

TEST(Simplex2SmallTests, Simple2) {
  auto p = make_problem(sparse<Rational>({{1, 1, -1, 1}, {1, 14, 10, -10}}),
                        Vector<Rational>{2, 24}, Vector<Rational>{1, 2, 3, -4},
                        {{Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}}});

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);
  auto& sol = std::get<FiniteLPSolution<Rational>>(result.solution);

  EXPECT_EQ(sol.point, (Vector<Rational>{4, 0, 2, 0}));
  EXPECT_EQ(sol.value, Rational{10});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, sol));
}

TEST(Simplex2SmallTests, Simple3) {
  auto p = make_problem(
      sparse<Rational>({{1, 1}}), Vector<Rational>{1}, Vector<Rational>{1, 2},
      {{Rational{0}, Rational{10}}, {Rational{0}, Rational{10}}});

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);
  auto& sol = std::get<FiniteLPSolution<Rational>>(result.solution);

  EXPECT_EQ(sol.point, (Vector<Rational>{0, 1}));
  EXPECT_EQ(sol.value, Rational{2});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, sol));
}

TEST(Simplex2SmallTests, Simple4) {
  auto p = make_problem(
      sparse<Rational>({{1, 1}}), Vector<Rational>{1}, Vector<Rational>{2, 1},
      {{Rational{0}, Rational{10}}, {Rational{0}, Rational{10}}});

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);
  auto& sol = std::get<FiniteLPSolution<Rational>>(result.solution);

  EXPECT_EQ(sol.point, (Vector<Rational>{1, 0}));
  EXPECT_EQ(sol.value, Rational{2});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, sol));
}

TEST(Simplex2SmallTests, TrivialOneDimensional) {
  // min/max x subject to x = 3, 0 <= x <= 10
  auto p = make_problem(sparse<Rational>({{1}}), Vector<Rational>{3},
                        Vector<Rational>{1}, {{Rational{0}, Rational{10}}});

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);
  auto& sol = std::get<FiniteLPSolution<Rational>>(result.solution);

  EXPECT_EQ(sol.point, (Vector<Rational>{3}));
  EXPECT_EQ(sol.value, Rational{3});
}

TEST(Simplex2SmallTests, NonTrivialBounds) {
  auto p = make_problem(sparse<Rational>({{1, -1, 1, 0}, {2, 1, 0, 1}}),
                        Vector<Rational>{1, 3}, Vector<Rational>{2, 1, 1, -1},
                        {{Rational{0}, Rational{1}},
                         {Rational{0}, Rational{3}},
                         {Rational{1}, Rational{10}},
                         {Rational{0}, Rational{10}}});

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);
  auto& sol = std::get<FiniteLPSolution<Rational>>(result.solution);

  EXPECT_EQ(sol.point, (Vector<Rational>{0, 3, 4, 0}));
  EXPECT_EQ(sol.value, Rational{7});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, sol));
}

TEST(Simplex2SmallTests, NonTrivialBounds2) {
  auto p = make_problem(sparse<Rational>({{1, 1, 1, 0, 0, 0}}),
                        Vector<Rational>{1}, Vector<Rational>{1, 0, 0, 1, 1, 1},
                        {{Rational{0}, Rational{1}},
                         {Rational{0}, Rational{1}},
                         {Rational{0}, Rational{1}},
                         {Rational{0}, Rational{1}},
                         {Rational{0}, Rational{1}},
                         {Rational{0}, Rational{1}}});

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);
  auto& sol = std::get<FiniteLPSolution<Rational>>(result.solution);

  EXPECT_EQ(sol.point, (Vector<Rational>{1, 0, 0, 1, 1, 1}));
  EXPECT_EQ(sol.value, Rational{4});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, sol));
}

TEST(Simplex2SmallTests, InfeasibleDetected) {
  // x1 + x2 = 5, but 0 <= x1 <= 1, 0 <= x2 <= 1 (max achievable is 2)
  auto p = make_problem(
      sparse<Rational>({{1, 1}}), Vector<Rational>{5}, Vector<Rational>{1, 1},
      {{Rational{0}, Rational{1}}, {Rational{0}, Rational{1}}});

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);
  EXPECT_TRUE(std::holds_alternative<NoFeasibleElements>(result.solution));
}

TEST(Simplex2SmallTests, UnboundedDetected) {
  // max x2 subject to x1 = 0, x1 in [0,1], x2 in [0, +inf)
  // x2 is decoupled (column of zeros), so problem is unbounded
  auto p = make_problem(
      sparse<Rational>({{1, 0}}), Vector<Rational>{0}, Vector<Rational>{0, 1},
      {{Rational{0}, Rational{1}}, {Rational{0}, std::nullopt}});

  // Primal init: x1 basic (at value 0 = rhs), x2 at lower bound
  std::vector<VariableState> states = {VariableState::BASIC,
                                       VariableState::AT_LOWER};

  ASSERT_TRUE(simplex::is_primal_feasible(p, states));

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto result = solver.primal(states);
  EXPECT_TRUE(std::holds_alternative<Unbounded>(result.solution));
}

TEST(Simplex2SmallTests, PrimalSimplex) {
  // Textbook primal simplex problem: max x1 - x2 + x3
  // subject to  2x1 - x2 + 2x3 + x4              = 4
  //             2x1 - 3x2 +  x3       + x5        = -5 (note negative rhs)
  //              -x1 +  x2 - 2x3            + x6  = -1
  // All vars >= 0 (no upper bound), initial basis = {1, 2, 3} (x2, x3, x4)
  auto p = make_problem(
      sparse<Rational>(
          {{2, -1, 2, 1, 0, 0}, {2, -3, 1, 0, 1, 0}, {-1, 1, -2, 0, 0, 1}}),
      Vector<Rational>{4, -5, -1}, Vector<Rational>{1, -1, 1, 0, 0, 0},
      {{Rational{0}, std::nullopt},
       {Rational{0}, std::nullopt},
       {Rational{0}, std::nullopt},
       {Rational{0}, std::nullopt},
       {Rational{0}, std::nullopt},
       {Rational{0}, std::nullopt}});

  // Primal feasible initial states: x2, x3, x4 basic (indices 1, 2, 3)
  std::vector<VariableState> states = {
      VariableState::AT_LOWER, VariableState::BASIC,    VariableState::BASIC,
      VariableState::BASIC,    VariableState::AT_LOWER, VariableState::AT_LOWER,
  };

  ASSERT_TRUE(simplex::is_primal_feasible(p, states));

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto result = solver.primal(states);
  ASSERT_TRUE(result.is_feasible());

  auto& sol = std::get<FiniteLPSolution<Rational>>(result.solution);
  EXPECT_EQ(sol.point,
            (Vector<Rational>{0, Rational{14} / 5, Rational{17} / 5, 0, 0, 3}));
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, sol));
}

TEST(Simplex2SmallTests, SetProblemTwice) {
  // Verify that calling set_problem twice gives correct results for the new
  // problem (the internal LUPA is re-initialized).
  auto p1 = make_problem(
      sparse<Rational>({{1, 1}}), Vector<Rational>{1}, Vector<Rational>{2, 1},
      {{Rational{0}, Rational{10}}, {Rational{0}, Rational{10}}});

  auto p2 = make_problem(
      sparse<Rational>({{1, 1}}), Vector<Rational>{1}, Vector<Rational>{1, 2},
      {{Rational{0}, Rational{10}}, {Rational{0}, Rational{10}}});

  simplex::Simplex<Rational> solver;

  // Solve p1
  solver.set_problem(p1);
  auto s1 = dual_init(p1);
  ASSERT_TRUE(s1.has_value());
  auto r1 = solver.dual(*s1);
  EXPECT_EQ(std::get<FiniteLPSolution<Rational>>(r1.solution).value,
            Rational{2});

  // Solve p2 with the same solver
  solver.set_problem(p2);
  auto s2 = dual_init(p2);
  ASSERT_TRUE(s2.has_value());
  auto r2 = solver.dual(*s2);
  auto& sol2 = std::get<FiniteLPSolution<Rational>>(r2.solution);
  EXPECT_EQ(sol2.point, (Vector<Rational>{0, 1}));
  EXPECT_EQ(sol2.value, Rational{2});
}
