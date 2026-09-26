#include <gtest/gtest.h>

#include "Assertions.h"
#include "ConstructSparse.h"
#include "problem/StandardMILP.h"
#include "simplex/Feasibility.h"
#include "simplex/Simplex.h"
#include "simplex/init/dual/ReducedCost.h"
#include "support/GMPRational.h"

// Convenience: build a StandardMILP<GMPRational> from dense matrix
// representation
static problem::StandardMILP<GMPRational> make_problem(
    CSCMatrix<GMPRational> matrix, Vector<GMPRational> rhs,
    Vector<GMPRational> cost, std::vector<Bound<GMPRational>> bounds) {
  const auto [n, d] = matrix.shape();
  problem::StandardMILP<GMPRational> p;
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
static std::optional<std::vector<simplex::VariableState>> dual_init(
    const problem::StandardMILP<GMPRational>& p) {
  return simplex::try_init_dual_by_reduced_cost(p.matrix, p.rhs, p.cost,
                                                p.var_bounds);
}

// ---- Adapted from tests/simplex/SmallTests.cpp ----

TEST(Simplex2SmallTests, StartingInSolution) {
  auto p =
      make_problem(sparse<GMPRational>({{1, -1, 1, 0}, {2, 1, 0, 1}}),
                   Vector<GMPRational>{1, 3}, Vector<GMPRational>{2, 1, 1, -1},
                   {{GMPRational{0}, GMPRational{10}},
                    {GMPRational{0}, GMPRational{10}},
                    {GMPRational{0}, GMPRational{10}},
                    {GMPRational{0}, GMPRational{10}}});

  simplex::Simplex<GMPRational> solver;
  solver.set_problem(p);

  auto states = simplex::try_init_dual_by_reduced_cost(
      p.matrix, p.rhs, p.cost, p.var_bounds, std::vector<size_t>{1, 2});
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);

  EXPECT_EQ(solver.get_point(), (Vector<GMPRational>{0, 3, 4, 0}));
  EXPECT_EQ(*result.objective, GMPRational{7});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, solver));
}

TEST(Simplex2SmallTests, Simple1) {
  auto p =
      make_problem(sparse<GMPRational>({{1, -1, 1, 0}, {2, 1, 0, 1}}),
                   Vector<GMPRational>{1, 3}, Vector<GMPRational>{2, 1, 1, -1},
                   {{GMPRational{0}, GMPRational{10}},
                    {GMPRational{0}, GMPRational{10}},
                    {GMPRational{0}, GMPRational{10}},
                    {GMPRational{0}, GMPRational{10}}});

  simplex::Simplex<GMPRational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);

  EXPECT_EQ(solver.get_point(), (Vector<GMPRational>{0, 3, 4, 0}));
  EXPECT_EQ(result.objective, GMPRational{7});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, solver));
}

TEST(Simplex2SmallTests, DISABLED_Simple2) {
  auto p =
      make_problem(sparse<GMPRational>({{1, 1, -1, 1}, {1, 14, 10, -10}}),
                   Vector<GMPRational>{2, 24}, Vector<GMPRational>{1, 2, 3, -4},
                   {{GMPRational{0}, GMPRational{10}},
                    {GMPRational{0}, GMPRational{10}},
                    {GMPRational{0}, GMPRational{10}},
                    {GMPRational{0}, GMPRational{10}}});

  simplex::Simplex<GMPRational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);

  EXPECT_EQ(solver.get_point(), (Vector<GMPRational>{4, 0, 2, 0}));
  EXPECT_EQ(result.objective, GMPRational{10});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, solver));
}

TEST(Simplex2SmallTests, Simple3) {
  auto p = make_problem(
      sparse<GMPRational>({{1, 1}}), Vector<GMPRational>{1},
      Vector<GMPRational>{1, 2},
      {{GMPRational{0}, GMPRational{10}}, {GMPRational{0}, GMPRational{10}}});

  simplex::Simplex<GMPRational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);

  EXPECT_EQ(solver.get_point(), (Vector<GMPRational>{0, 1}));
  EXPECT_EQ(result.objective, GMPRational{2});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, solver));
}

TEST(Simplex2SmallTests, Simple4) {
  auto p = make_problem(
      sparse<GMPRational>({{1, 1}}), Vector<GMPRational>{1},
      Vector<GMPRational>{2, 1},
      {{GMPRational{0}, GMPRational{10}}, {GMPRational{0}, GMPRational{10}}});

  simplex::Simplex<GMPRational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);

  EXPECT_EQ(solver.get_point(), (Vector<GMPRational>{1, 0}));
  EXPECT_EQ(result.objective, GMPRational{2});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, solver));
}

TEST(Simplex2SmallTests, TrivialOneDimensional) {
  // min/max x subject to x = 3, 0 <= x <= 10
  auto p =
      make_problem(sparse<GMPRational>({{1}}), Vector<GMPRational>{3},
                   Vector<GMPRational>{1}, {{GMPRational{0}, GMPRational{10}}});

  simplex::Simplex<GMPRational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);

  EXPECT_EQ(solver.get_point(), (Vector<GMPRational>{3}));
  EXPECT_EQ(result.objective, GMPRational{3});
}

TEST(Simplex2SmallTests, NonTrivialBounds) {
  auto p =
      make_problem(sparse<GMPRational>({{1, -1, 1, 0}, {2, 1, 0, 1}}),
                   Vector<GMPRational>{1, 3}, Vector<GMPRational>{2, 1, 1, -1},
                   {{GMPRational{0}, GMPRational{1}},
                    {GMPRational{0}, GMPRational{3}},
                    {GMPRational{1}, GMPRational{10}},
                    {GMPRational{0}, GMPRational{10}}});

  simplex::Simplex<GMPRational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);

  EXPECT_EQ(solver.get_point(), (Vector<GMPRational>{0, 3, 4, 0}));
  EXPECT_EQ(result.objective, GMPRational{7});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, solver));
}

TEST(Simplex2SmallTests, NonTrivialBounds2) {
  auto p = make_problem(sparse<GMPRational>({{1, 1, 1, 0, 0, 0}}),
                        Vector<GMPRational>{1},
                        Vector<GMPRational>{1, 0, 0, 1, 1, 1},
                        {{GMPRational{0}, GMPRational{1}},
                         {GMPRational{0}, GMPRational{1}},
                         {GMPRational{0}, GMPRational{1}},
                         {GMPRational{0}, GMPRational{1}},
                         {GMPRational{0}, GMPRational{1}},
                         {GMPRational{0}, GMPRational{1}}});

  simplex::Simplex<GMPRational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);

  EXPECT_EQ(solver.get_point(), (Vector<GMPRational>{1, 0, 0, 1, 1, 1}));
  EXPECT_EQ(result.objective, GMPRational{4});
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, solver));
}

TEST(Simplex2SmallTests, DISABLED_InfeasibleDetected) {
  // x1 + x2 = 5, but 0 <= x1 <= 1, 0 <= x2 <= 1 (max achievable is 2)
  auto p = make_problem(
      sparse<GMPRational>({{1, 1}}), Vector<GMPRational>{5},
      Vector<GMPRational>{1, 1},
      {{GMPRational{0}, GMPRational{1}}, {GMPRational{0}, GMPRational{1}}});

  simplex::Simplex<GMPRational> solver;
  solver.set_problem(p);

  auto states = dual_init(p);
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);
  EXPECT_EQ(result.status, simplex::Status::INFEASIBLE);
}

TEST(Simplex2SmallTests, UnboundedDetected) {
  // max x2 subject to x1 = 0, x1 in [0,1], x2 in [0, +inf)
  // x2 is decoupled (column of zeros), so problem is unbounded
  auto p = make_problem(
      sparse<GMPRational>({{1, 0}}), Vector<GMPRational>{0},
      Vector<GMPRational>{0, 1},
      {{GMPRational{0}, GMPRational{1}}, {GMPRational{0}, std::nullopt}});

  // Primal init: x1 basic (at value 0 = rhs), x2 at lower bound
  std::vector states = {simplex::VariableState::BASIC,
                        simplex::VariableState::AT_LOWER};

  ASSERT_TRUE(simplex::is_primal_feasible(p, states));

  simplex::Simplex<GMPRational> solver;
  solver.set_problem(p);

  auto result = solver.primal(states);
  ASSERT_EQ(result.status, simplex::Status::UNBOUNDED);
}

TEST(Simplex2SmallTests, PrimalSimplex) {
  // Textbook primal simplex problem: max x1 - x2 + x3
  // subject to  2x1 - x2 + 2x3 + x4              = 4
  //             2x1 - 3x2 +  x3       + x5        = -5 (note negative rhs)
  //              -x1 +  x2 - 2x3            + x6  = -1
  // All vars >= 0 (no upper bound), initial basis = {1, 2, 3} (x2, x3, x4)
  auto p = make_problem(
      sparse<GMPRational>(
          {{2, -1, 2, 1, 0, 0}, {2, -3, 1, 0, 1, 0}, {-1, 1, -2, 0, 0, 1}}),
      Vector<GMPRational>{4, -5, -1}, Vector<GMPRational>{1, -1, 1, 0, 0, 0},
      {{GMPRational{0}, std::nullopt},
       {GMPRational{0}, std::nullopt},
       {GMPRational{0}, std::nullopt},
       {GMPRational{0}, std::nullopt},
       {GMPRational{0}, std::nullopt},
       {GMPRational{0}, std::nullopt}});

  // Primal feasible initial states: x2, x3, x4 basic (indices 1, 2, 3)
  std::vector states = {
      simplex::VariableState::AT_LOWER, simplex::VariableState::BASIC,
      simplex::VariableState::BASIC,    simplex::VariableState::BASIC,
      simplex::VariableState::AT_LOWER, simplex::VariableState::AT_LOWER,
  };

  ASSERT_TRUE(simplex::is_primal_feasible(p, states));

  simplex::Simplex<GMPRational> solver;
  solver.set_problem(p);

  auto result = solver.primal(states);
  ASSERT_EQ(result.status, simplex::Status::OPTIMAL);

  EXPECT_EQ(solver.get_point(),
            (Vector<GMPRational>{0, GMPRational{14} / 5, GMPRational{17} / 5, 0,
                                 0, 3}));
  ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(p, solver));
}

TEST(Simplex2SmallTests, SetProblemTwice) {
  // Verify that calling set_problem twice gives correct results for the new
  // problem (the internal LUPA is re-initialized).
  auto p1 = make_problem(
      sparse<GMPRational>({{1, 1}}), Vector<GMPRational>{1},
      Vector<GMPRational>{2, 1},
      {{GMPRational{0}, GMPRational{10}}, {GMPRational{0}, GMPRational{10}}});

  auto p2 = make_problem(
      sparse<GMPRational>({{1, 1}}), Vector<GMPRational>{1},
      Vector<GMPRational>{1, 2},
      {{GMPRational{0}, GMPRational{10}}, {GMPRational{0}, GMPRational{10}}});

  simplex::Simplex<GMPRational> solver;

  // Solve p1
  solver.set_problem(p1);
  auto s1 = dual_init(p1);
  ASSERT_TRUE(s1.has_value());
  auto r1 = solver.dual(*s1);
  EXPECT_EQ(r1.objective, GMPRational{2});

  // Solve p2 with the same solver
  solver.set_problem(p2);
  auto s2 = dual_init(p2);
  ASSERT_TRUE(s2.has_value());
  auto r2 = solver.dual(*s2);
  EXPECT_EQ(solver.get_point(), (Vector<GMPRational>{0, 1}));
  EXPECT_EQ(r2.objective, GMPRational{2});
}
