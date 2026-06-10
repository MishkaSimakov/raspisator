#pragma once

#include <vector>

#include "linalg/Linalg.h"
#include "linalg/Permutation.h"
#include "presolve/Pass.h"
#include "utils/Accumulators.h"

namespace presolve {

template <typename Field>
class RemoveLinearlyDependentEqualities final : public Pass<Field> {
  // Performs the first part of row reduction. Ignores rows with non-zero range.
  void row_reduction(linalg::Matrix<Field>& matrix,
                     std::vector<Bound<Field>>& rhs_bounds) {
    using std::abs;

    const auto [n, d] = matrix.shape();

    auto permutation = linalg::Permutation::id(n);
    size_t current_row = 0;

    for (size_t col = 0; col < d; ++col) {
      ArgMaximum<Field> max_abs;

      for (size_t row = current_row; row < n; ++row) {
        if (rhs_bounds[row].is_fixed()) {
          max_abs.record(row, abs(matrix[permutation.apply(row), col]));
        }
      }

      if (!max_abs.has_value() ||
          !FieldTraits<Field>::is_nonzero(max_abs->max)) {
        continue;
      }

      permutation.swap(max_abs->index, current_row);

      const Field pivot = matrix[permutation.apply(current_row), col];

      for (size_t row = current_row + 1; row < n; ++row) {
        if (!rhs_bounds[permutation.apply(row)].is_fixed()) {
          continue;
        }

        const Field value = matrix[permutation.apply(row), col];

        if (!FieldTraits<Field>::is_nonzero(value)) {
          continue;
        }

        rhs_bounds[permutation.apply(row)] -=
            rhs_bounds[permutation.apply(current_row)] * value / pivot;

        for (size_t j = 0; j < d; ++j) {
          matrix[permutation.apply(row), j] -=
              matrix[permutation.apply(current_row), j] * value / pivot;
        }
        matrix[permutation.apply(row), col] = 0;
      }

      ++current_row;
    }
  }

 public:
  RemoveLinearlyDependentEqualities() = default;

  problem::MILP<Field> apply(problem::MILP<Field> problem) override {
    this->register_apply();

    auto matrix = Matrix(problem.matrix);
    auto bounds = problem.rhs_bounds;

    row_reduction(matrix, bounds);

    // find linearly dependent rows using matrix after row reduction
    const auto [n, d] = matrix.shape();
    std::vector<size_t> rows_mapping(n, n);
    size_t new_rows = 0;

    for (size_t row = 0; row < n; ++row) {
      if (!problem.rhs_bounds[row].is_fixed()) {
        rows_mapping[row] = new_rows++;
        continue;
      }

      bool is_empty = true;

      for (size_t col = 0; col < d; ++col) {
        if (FieldTraits<Field>::is_nonzero(matrix[row, col])) {
          is_empty = false;
          break;
        }
      }

      if (!is_empty) {
        rows_mapping[row] = new_rows++;
        continue;
      }

      if (!bounds[row].contains(0)) {
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

  Vector<Field> inverse(Vector<Field> solution) const override {
    return solution;
  }
};

}  // namespace presolve
