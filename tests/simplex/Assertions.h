#pragma once

#include <gtest/gtest.h>

#include <algorithm>

#include "linalg/Matrix.h"
#include "linalg/Transpose.h"
#include "linalg/Vector.h"
#include "linalg/lu/LUPA.h"
#include "linear/model/LP.h"

using linalg::Matrix, linalg::Vector, linalg::CSCMatrix;

template <typename Field>
void validate_simplex_solution(const Matrix<Field>& A, const Vector<Field>& b,
                               const Vector<Field>& c,
                               const Bounds<Field>& bounds,
                               const FiniteLPSolution<Field>& solution) {
  auto [n, d] = A.shape();

  Vector residue = A * solution.point - b;
  for (size_t i = 0; i < n; ++i) {
    ASSERT_TRUE(!FieldTraits<Field>::is_nonzero(residue[i]));
  }

  for (size_t i = 0; i < d; ++i) {
    if (bounds[i].lower) {
      ASSERT_TRUE(!FieldTraits<Field>::is_strictly_positive(*bounds[i].lower -
                                                            solution.point[i]));
    }
    if (bounds[i].upper) {
      ASSERT_TRUE(!FieldTraits<Field>::is_strictly_positive(solution.point[i] -
                                                            *bounds[i].upper));
    }

    for (size_t i = 0; i < d; ++i) {
      if (solution.variables[i] == VariableState::AT_LOWER) {
        ASSERT_TRUE(bounds[i].lower &&
                    !FieldTraits<Field>::is_nonzero(*bounds[i].lower -
                                                    solution.point[i]));
      } else if (solution.variables[i] == VariableState::AT_UPPER) {
        ASSERT_TRUE(bounds[i].upper &&
                    !FieldTraits<Field>::is_nonzero(*bounds[i].upper -
                                                    solution.point[i]));
      }
    }
  }

  // check optimality
  auto sparse_A = CSCMatrix(A);
  auto basic_vars = solution.get_basic_variables();

  auto [P, Q, ls, us] =
      linalg::FullPivotingLU<Field>(n).get(sparse_A, basic_vars);

  auto pi =
      linalg::solve_linear_transposed(Vector(c[basic_vars]), P, Q, ls, us);

  Vector reduced_cost = c - linalg::transpose(sparse_A) * pi;

  for (size_t i = 0; i < d; ++i) {
    if (solution.variables[i] == VariableState::BASIC) {
      continue;
    }

    ASSERT_TRUE((solution.variables[i] == VariableState::AT_LOWER &&
                 !FieldTraits<Field>::is_strictly_positive(reduced_cost[i])) ||
                (solution.variables[i] == VariableState::AT_UPPER &&
                 !FieldTraits<Field>::is_strictly_negative(reduced_cost[i])));
  }
}
