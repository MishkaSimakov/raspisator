#pragma once

#include <gtest/gtest.h>

#include <algorithm>

#include "linear/matrix/Matrix.h"
#include "linear/model/LP.h"
#include "linear/sparse/LU.h"

template <typename Field>
void validate_simplex_solution(const Matrix<Field>& A, const Matrix<Field>& b,
                               const Matrix<Field>& c,
                               const Bounds<Field>& bounds,
                               const FiniteLPSolution<Field>& solution) {
  auto [n, d] = A.shape();

  auto residue = A * solution.point - b;
  for (size_t i = 0; i < n; ++i) {
    ASSERT_TRUE(!FieldTraits<Field>::is_nonzero(residue[i, 0]));
  }

  for (size_t i = 0; i < d; ++i) {
    if (bounds[i].lower) {
      ASSERT_TRUE(!FieldTraits<Field>::is_strictly_positive(
          *bounds[i].lower - solution.point[i, 0]));
    }
    if (bounds[i].upper) {
      ASSERT_TRUE(!FieldTraits<Field>::is_strictly_positive(
          solution.point[i, 0] - *bounds[i].upper));
    }

    for (size_t i = 0; i < d; ++i) {
      if (solution.variables[i] == VariableState::AT_LOWER) {
        ASSERT_TRUE(bounds[i].lower &&
                    !FieldTraits<Field>::is_nonzero(*bounds[i].lower -
                                                    solution.point[i, 0]));
      } else if (solution.variables[i] == VariableState::AT_UPPER) {
        ASSERT_TRUE(bounds[i].upper &&
                    !FieldTraits<Field>::is_nonzero(*bounds[i].upper -
                                                    solution.point[i, 0]));
      }
    }
  }

  // check optimality
  auto sparse_A = CSCMatrix(A);
  auto basic_vars = solution.get_basic_variables();

  auto [P, Q, ls, us] =
      linalg::FullPivotingLU<Field>(n).get(sparse_A, basic_vars);

  Matrix<Field> cb(n, 1);
  for (size_t i = 0; i < n; ++i) {
    cb[i, 0] = c[0, basic_vars[i]];
  }
  cb = linalg::solve_linear_transposed(cb, P, Q, ls, us);

  for (size_t i = 0; i < d; ++i) {
    if (solution.variables[i] == VariableState::BASIC) {
      continue;
    }

    Field reduced_cost = c[0, i];

    for (const auto& [row, value] : sparse_A.get_column(i)) {
      reduced_cost -= value * cb[row, 0];
    }

    ASSERT_TRUE((solution.variables[i] == VariableState::AT_LOWER &&
                 !FieldTraits<Field>::is_strictly_positive(reduced_cost)) ||
                (solution.variables[i] == VariableState::AT_UPPER &&
                 !FieldTraits<Field>::is_strictly_negative(reduced_cost)));
  }
}
