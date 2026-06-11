#include <gtest/gtest.h>

#include <variant>

#include "ConstructSparse.h"
#include "linalg/Matrix.h"
#include "linear/BigInteger.h"
#include "linear/model/LP.h"
#include "linear/simplex/Simplex.h"
#include "linear/simplex/init/dual/ReducedCost.h"
#include "problem/StandardMILP.h"

using linalg::CSCMatrix, linalg::Matrix, linalg::Vector;

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

// After change_basis(k, j, leaving_state):
//   - basic_vars[k] becomes j
//   - var_states[j] becomes BASIC
//   - var_states[old_basic_vars[k]] becomes leaving_state
TEST(Simplex2TraversalTests, ChangeBasisUpdatesVectors) {
  auto p = make_problem(sparse<Rational>({{1, -1, 1, 0}, {2, 1, 0, 1}}),
                        Vector<Rational>{1, 3}, Vector<Rational>{2, 1, 1, -1},
                        {{Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}}});

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto states = simplex::try_init_dual_by_reduced_cost(
      p.matrix, p.rhs, p.cost, Bounds<Rational>(p.var_bounds),
      std::vector<size_t>{1, 2});
  ASSERT_TRUE(states.has_value());

  // Run to a feasible basis
  auto result = solver.dual(*states);
  ASSERT_TRUE(result.is_feasible());

  const auto basic_vars_before = solver.get_basic_vars();
  const auto states_before = solver.get_states();
  const size_t n = basic_vars_before.size();
  const size_t d = states_before.size();

  // Find a nonbasic AT_LOWER variable to bring in
  size_t entering = d;  // sentinel
  for (size_t j = 0; j < d; ++j) {
    if (states_before[j] == VariableState::AT_LOWER) {
      entering = j;
      break;
    }
  }
  ASSERT_NE(entering, d) << "No AT_LOWER nonbasic variable found";

  // Use basis position 0 as the leaving slot
  const size_t leaving_pos = 0;
  const size_t leaving_var = basic_vars_before[leaving_pos];
  ASSERT_NE(leaving_var, entering);

  solver.change_basis(leaving_pos, entering, VariableState::AT_LOWER);

  const auto basic_vars_after = solver.get_basic_vars();
  const auto states_after = solver.get_states();

  // basic_vars[leaving_pos] is now the entering variable
  EXPECT_EQ(basic_vars_after[leaving_pos], entering);

  // Other basis positions unchanged
  for (size_t k = 0; k < n; ++k) {
    if (k != leaving_pos) {
      EXPECT_EQ(basic_vars_after[k], basic_vars_before[k]);
    }
  }

  // States updated correctly
  EXPECT_EQ(states_after[entering], VariableState::BASIC);
  EXPECT_EQ(states_after[leaving_var], VariableState::AT_LOWER);

  // States of other variables unchanged
  for (size_t j = 0; j < d; ++j) {
    if (j != entering && j != leaving_var) {
      EXPECT_EQ(states_after[j], states_before[j]);
    }
  }
}

// change_bound(j, new_state) simply updates var_states_[j].
TEST(Simplex2TraversalTests, ChangeBoundFlipsState) {
  auto p = make_problem(sparse<Rational>({{1, -1, 1, 0}, {2, 1, 0, 1}}),
                        Vector<Rational>{1, 3}, Vector<Rational>{2, 1, 1, -1},
                        {{Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}}});

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto states = simplex::try_init_dual_by_reduced_cost(
      p.matrix, p.rhs, p.cost, Bounds<Rational>(p.var_bounds));
  ASSERT_TRUE(states.has_value());

  auto result = solver.dual(*states);
  ASSERT_TRUE(result.is_feasible());

  const auto states_after_solve = solver.get_states();
  const size_t d = states_after_solve.size();

  // Find an AT_LOWER nonbasic variable
  size_t target = d;
  for (size_t j = 0; j < d; ++j) {
    if (states_after_solve[j] == VariableState::AT_LOWER) {
      target = j;
      break;
    }
  }
  ASSERT_NE(target, d) << "No AT_LOWER variable found after solve";

  solver.change_bound(target, VariableState::AT_UPPER);

  const auto updated_states = solver.get_states();
  EXPECT_EQ(updated_states[target], VariableState::AT_UPPER);

  // All other states unchanged
  for (size_t j = 0; j < d; ++j) {
    if (j != target) {
      EXPECT_EQ(updated_states[j], states_after_solve[j]);
    }
  }
}

// Verify get_tableau_row returns a row of the basis inverse times A.
// After solving, the tableau row for basic variable at position k should
// give B^{-1} A (k-th row). In particular, B^{-1} A[:, basic_vars[k]]
// should equal e_k (k-th unit vector).
TEST(Simplex2TraversalTests, TableauRowBasisIdentity) {
  auto p = make_problem(sparse<Rational>({{1, -1, 1, 0}, {2, 1, 0, 1}}),
                        Vector<Rational>{1, 3}, Vector<Rational>{2, 1, 1, -1},
                        {{Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}},
                         {Rational{0}, Rational{10}}});

  simplex::Simplex<Rational> solver;
  solver.set_problem(p);

  auto states = simplex::try_init_dual_by_reduced_cost(
      p.matrix, p.rhs, p.cost, Bounds<Rational>(p.var_bounds));
  ASSERT_TRUE(states.has_value());

  solver.dual(*states);

  const auto basic_vars = solver.get_basic_vars();
  const size_t n = basic_vars.size();
  const Matrix<Rational> A(p.matrix);

  // For each basis row k: (B^-1 A)_k,basic_vars[j] should be delta(k,j)
  for (size_t k = 0; k < n; ++k) {
    const auto row = solver.get_tableau_row(k);

    for (size_t j = 0; j < n; ++j) {
      // Compute the dot product of tableau row k with column basic_vars[j]
      Rational dot = 0;
      for (const auto& [row_idx, val] : p.matrix.get_column(basic_vars[j])) {
        dot += row[row_idx] * val;
      }
      if (j == k) {
        EXPECT_EQ(dot, Rational{1}) << "Diagonal entry of B^-1 B must be 1 at ("
                                    << k << "," << j << ")";
      } else {
        EXPECT_EQ(dot, Rational{0})
            << "Off-diagonal entry of B^-1 B must be 0 at (" << k << "," << j
            << ")";
      }
    }
  }
}
