#pragma once

#include <vector>

#include "linear/matrix/Elimination.h"
#include "linear/sparse/Permutation.h"
#include "presolve/Pass.h"
#include "utils/Accumulators.h"

namespace presolve {

template <typename Field>
class RemoveLinearlyDependentConstraints final : public Pass<Field> {
  void gauss_elimination(Matrix<Field>& matrix,
                         std::vector<Bound<Field>>& rhs_bounds) {
    using std::abs;

    const auto [n, d] = matrix.shape();

    auto permutation = linalg::Permutation::id(n);
    size_t current_row = 0;

    for (size_t col = 0; col < d; ++col) {
      ArgMaximum<Field> max_abs;

      for (size_t row = current_row; row < n; ++row) {
        max_abs.record(row, abs(matrix[permutation.apply(row), col]));
      }

      if (!max_abs.has_value() ||
          !FieldTraits<Field>::is_nonzero(max_abs->max)) {
        continue;
      }

      permutation.swap(max_abs->index, current_row);

      for (size_t row = current_row + 1; row < n; ++row) {
        const Field coef = matrix[permutation.apply(row), col];

        if (!FieldTraits<Field>::is_nonzero(coef)) {
          continue;
        }

        rhs_bounds[row] -= rhs_bounds[permutation.apply(current_row)] * coef /
                           matrix[permutation.apply(current_row), col];

        matrix[permutation.apply(row), {0, d}].sub_mul(
            matrix[permutation.apply(current_row), {0, d}],
            coef / matrix[permutation.apply(current_row), col]);
        matrix[permutation.apply(current_row), col] = 0;
      }

      ++current_row;
    }
  }

 public:
  RemoveLinearlyDependentConstraints() = default;

  problem::MILP<Field> apply(problem::MILP<Field> problem) override {
    this->register_apply();

    auto matrix = linalg::to_dense(problem.matrix);
    auto bounds = problem.implied_var_bounds;

    gauss_elimination(matrix, bounds);

    // find linearly dependent rows using matrix after gaussian elimination
    const auto [n, d] = matrix.shape();
    std::vector<size_t> rows_mapping(n, n);
    size_t new_rows = 0;

    for (size_t row = 0; row < n; ++row) {
      bool is_empty = true;

      for (size_t col = 0; col < d; ++col) {
        if (FieldTraits<Field>::is_nonzero(matrix[row, col])) {
          is_empty = false;
          break;
        }
      }

      if (!is_empty) {
        rows_mapping[row] = new_rows;
        ++new_rows;
        continue;
      }

      if (!bounds[row].is_inside(0)) {
        problem.proven_infeasible = true;
        return problem;
      }
    }

    // construct new problem
    for (size_t col = 0; col < d; ++col) {
      for (auto& [row, value] : problem.matrix.get_column(col)) {
        row = rows_mapping[row];
      }
    }
    problem.matrix.resize(new_rows, d);

    for (size_t row = 0; row < n; ++row) {
      if (rows_mapping[row] != n) {
        problem.rhs_bounds[rows_mapping[row]] = problem.rhs_bounds[row];
        problem.row_names[rows_mapping[row]] = problem.row_names[row];
      }
    }

    problem.rhs_bounds.resize(new_rows);
    problem.row_names.resize(new_rows);

    return problem;
  }

  std::vector<Field> inverse(std::vector<Field> solution) const override {
    return solution;
  }
};

}  // namespace presolve
