#pragma once

#include <cmath>
#include <vector>

#include "linalg/Linalg.h"
#include "linalg/Permutation.h"
#include "linalg/RRQR.h"
#include "presolve/Pass.h"
#include "utils/Accumulators.h"
#include "utils/Logging.h"

namespace presolve {

template <typename Field>
class RemoveLinearlyDependentEqualities final : public Pass<Field> {
  const Field pivot_tolerance_;
  const Field feasibility_tolerance_;

  // Performs the first part of row reduction.
  void row_reduction(Matrix<Field>& matrix, std::vector<Field>& rhs_bounds) {
    using std::abs;

    const auto [n, d] = matrix.shape();

    std::vector<size_t> permutation(n);
    std::iota(permutation.begin(), permutation.end(), 0);

    size_t current_row = 0;

    for (size_t col = 0; col < d; ++col) {
      ArgMaximum<Field> max_abs;

      for (size_t row = current_row; row < n; ++row) {
        max_abs.record(row, abs(matrix[permutation[row], col]));
      }

      if (!max_abs.has_value() || abs(max_abs->max) <= pivot_tolerance_) {
        continue;
      }

      std::swap(permutation[max_abs->index], permutation[current_row]);

      const Field pivot = matrix[permutation[current_row], col];

      for (size_t row = current_row + 1; row < n; ++row) {
        const Field value = matrix[permutation[row], col];

        if (value == 0) {
          continue;
        }

        rhs_bounds[permutation[row]] -=
            rhs_bounds[permutation[current_row]] * value / pivot;

        for (size_t j = col + 1; j < d; ++j) {
          matrix[permutation[row], j] -=
              matrix[permutation[current_row], j] * value / pivot;
        }
        matrix[permutation[row], col] = 0;
      }

      ++current_row;
    }
  }

  bool is_zero_row(const Matrix<Field>& matrix, size_t row) const {
    for (size_t col = 0; col < matrix.cols(); ++col) {
      if (abs(matrix[row, col]) > pivot_tolerance_) {
        return false;
      }
    }

    return true;
  }

  static Vector<Field> solve(const Matrix<Field>& L, Vector<Field> rhs) {
    const size_t rank = rhs.size();

    for (size_t i = 0; i < rank; ++i) {
      Field r = rhs[i];

      for (size_t j = 0; j < i; ++j) {
        r -= L[i, j] * rhs[j];
      }

      rhs[i] = r / L[i, i];
    }

    return std::move(rhs);
  }

 public:
  explicit RemoveLinearlyDependentEqualities(Field pivot_tolerance = 1e-7,
                                             Field feasibility_tolerance = 1e-7)
      : pivot_tolerance_(pivot_tolerance),
        feasibility_tolerance_(feasibility_tolerance) {}

  problem::MILP<Field> apply(problem::MILP<Field> problem) override {
    using std::abs;

    this->register_apply();

    const auto [n, d] = problem.matrix.shape();

    std::vector<size_t> rows_mapping(n, n);
    size_t equalities_count = 0;

    for (size_t i = 0; i < n; ++i) {
      if (problem.rhs_bounds[i].is_fixed()) {
        rows_mapping[i] = equalities_count++;
      }
    }

    // submatrix formed by equalities
    auto rhs = Vector<Field>::zeros(equalities_count);
    auto matrix = Matrix<Field>::zeros(equalities_count, d);

    for (size_t row = 0; row < n; ++row) {
      if (rows_mapping[row] != n) {
        rhs[rows_mapping[row]] = *problem.rhs_bounds[row].lower;
      }
    }

    for (size_t col = 0; col < d; ++col) {
      for (const auto& [row, value] : problem.matrix.get_column(col)) {
        if (rows_mapping[row] != n) {
          matrix[rows_mapping[row], col] = value;
        }
      }
    }

    auto rrqr_result = linalg::rrqr(matrix, pivot_tolerance_);

    std::vector<size_t> inverse_permutation(equalities_count);
    for (size_t i = 0; i < equalities_count; ++i) {
      inverse_permutation[rrqr_result.permutation[i]] = i;
    }

    // verify that rhs is feasible
    Vector<Field> z(rrqr_result.rank);

    for (size_t i = 0; i < rrqr_result.rank; ++i) {
      z[i] = rhs[rrqr_result.permutation[i]];
    }

    z = solve(rrqr_result.R, std::move(z));

    for (size_t i = rrqr_result.rank; i < equalities_count; ++i) {
      Field expected_rhs = 0;

      for (size_t j = 0; j < rrqr_result.rank; ++j) {
        expected_rhs += rrqr_result.R[i, j] * z[j];
      }

      if (abs(rhs[rrqr_result.permutation[i]] - expected_rhs) >
          feasibility_tolerance_) {
        problem.proven_infeasible = true;
        return problem;
      }
    }

    // find linearly dependent rows using matrix after row reduction
    size_t new_rows_count = 0;

    for (size_t i = 0; i < n; ++i) {
      if (!problem.rhs_bounds[i].is_fixed()) {
        rows_mapping[i] = new_rows_count++;
        continue;
      }

      if (inverse_permutation[rows_mapping[i]] < rrqr_result.rank) {
        rows_mapping[i] = new_rows_count++;
        continue;
      }

      rows_mapping[i] = n;
    }

    // construct new problem
    for (size_t col = 0; col < d; ++col) {
      for (auto& [row, value] : problem.matrix.get_column(col)) {
        row = rows_mapping[row];
      }
    }
    problem.matrix.resize(new_rows_count, d);

    for (size_t row = 0; row < n; ++row) {
      if (rows_mapping[row] != n) {
        problem.rhs_bounds[rows_mapping[row]] = problem.rhs_bounds[row];
        problem.row_names[rows_mapping[row]] = problem.row_names[row];
      }
    }

    problem.rhs_bounds.resize(new_rows_count);
    problem.row_names.resize(new_rows_count);

    return problem;
  }

  Vector<Field> inverse(Vector<Field> solution) const override {
    return solution;
  }
};

}  // namespace presolve
