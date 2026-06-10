#pragma once

#include <cmath>
#include <ranges>
#include <unordered_set>
#include <vector>

#include "EtaFile.h"
#include "SCC.h"
#include "linalg/lu/FullPivotingLU.h"
#include "linalg/lu/Solve.h"
#include "utils/Accumulators.h"

namespace linalg {

struct LUPAConfig {
  size_t purge_after_iterations{100};
  size_t refactorize_after_iterations{500};
};

// LUP-Accelerated (LUPA)
// For a given matrix this class answers queries of form:
// get LUP-decomposition of submatrix of A formed by given columns
template <typename Field>
class LUPA {
  const CSCMatrix<Field>& A_;
  std::vector<size_t> columns_;

  FullPivotingLU<Field> factorizer_;

  // decomposition
  EtaFile<Field> us_;  // stores U1^-1 ... U(n+t)^-1
  EtaFile<Field> ls_;  // stores L1^-1 ... Ln^-1 R1^-1 ... Rt^-1
  Permutation P_;
  Permutation Q_;

  // Forrest-Tomlin update helpers
  size_t changes_since_refactorization_{0};
  size_t changes_since_purge_{0};

  const LUPAConfig config_;

  void purge() {
    ls_.purge();
    us_.purge();

    changes_since_purge_ = 0;
  }

  enum class UpdateResult { SUCCESS, NEED_REFACTORIZATION };

  UpdateResult forrest_tomlin_update(size_t current_column, size_t new_column) {
    auto [n, d] = A_.shape();

    Matrix<Field> column(n, 1, 0);
    for (const auto [row, value] : A_.get_column(new_column)) {
      column[row, 0] = value;
    }

    column = P_.apply(std::move(column));
    for (auto entry : ls_) {
      column = entry.apply(std::move(column));
    }

    // Permutation Q_ changes columns order. Here we find a column that becomes
    // current_column after applying Q_
    current_column = Q_.post_apply(current_column);

    Matrix<Field> r(n, 1, 0);

    auto itr = us_.begin();
    for (; itr != us_.end(); ++itr) {
      if ((*itr).pivot_index() == current_column) {
        for (auto [row, value] : (*itr).pivot_entries()) {
          if (row == current_column) {
            break;
          }
        }

        itr = us_.erase(itr);
        break;
      }
    }

    // TODO: it was noticed that r is often empty in setcover problem!
    // Check this for other problems

    for (; itr != us_.end(); ++itr) {
      Field diagonal = 0;
      Field main_value = 0;

      for (auto& [row, value] : (*itr).pivot_entries()) {
        if (row == current_column) {
          main_value = value;
          value = 0;
        } else if (row == (*itr).pivot_index()) {
          diagonal = value;
        }
      }

      if (!FieldTraits<Field>::is_nonzero(diagonal)) {
        return UpdateResult::NEED_REFACTORIZATION;
      }

      r[(*itr).pivot_index(), 0] = main_value / diagonal;
      r = (*itr).apply_transposed(std::move(r));
    }

    r[current_column, 0] = 1;

    Maximum<Field> r_max;
    for (size_t i = 0; i < n; ++i) {
      r_max.record(r[i, 0]);
    }

    if (*r_max > 10) {
      // std::println("refactorization: {}", *r_max.max());
      return UpdateResult::NEED_REFACTORIZATION;
    }

    // logging::log_value(*r_max.max(), "r_max.txt");

    ls_.push_back(current_column, r, EtaMatrixType::ROW);

    // add new eta matrix to U1 ... Un
    column = (*std::prev(ls_.cend())).apply(std::move(column));

    if (!FieldTraits<Field>::is_nonzero(column[current_column, 0])) {
      return UpdateResult::NEED_REFACTORIZATION;
    }

    const Field diagonal = column[current_column, 0];

    for (size_t i = 0; i < n; ++i) {
      column[i, 0] =
          i != current_column ? -column[i, 0] / diagonal : Field(1) / diagonal;
    }

    us_.push_back(current_column, column, EtaMatrixType::COLUMN);

    return UpdateResult::SUCCESS;
  }

 public:
  explicit LUPA(const CSCMatrix<Field>& A, LUPAConfig config = {})
      : A_(A),
        factorizer_(A.shape().first),
        us_(A.rows()),
        ls_(A.rows()),
        P_(Permutation::id(A.rows())),
        Q_(Permutation::id(A.rows())),
        config_(config) {}

  void set_columns(const std::vector<size_t>& columns) {
    assert(columns.size() == A_.shape().first);

    columns_ = columns;
    refactorize();
  }

  void change_column(size_t current_column, size_t new_column) {
    columns_[current_column] = new_column;
    ++changes_since_refactorization_;
    ++changes_since_purge_;

    if (changes_since_refactorization_ > config_.refactorize_after_iterations) {
      refactorize();
      return;
    }

    auto update_result = forrest_tomlin_update(current_column, new_column);

    if (update_result == UpdateResult::NEED_REFACTORIZATION) {
      refactorize();
      return;
    }

    if (changes_since_purge_ > config_.purge_after_iterations) {
      purge();
    }
  }

  void refactorize() {
    factorizer_.get(A_, columns_, P_, Q_, ls_, us_);

    changes_since_refactorization_ = 0;
    changes_since_purge_ = 0;
  }

  // solves Ax = b
  Vector<Field> solve_linear(Vector<Field> b) const {
    return linalg::solve_linear(std::move(b), P_, Q_, ls_, us_);
  }

  Vector<Field> solve_linear_transposed(Vector<Field> b) const {
    return linalg::solve_linear_transposed(std::move(b), P_, Q_, ls_, us_);
  }

  Matrix<Field> get_row(size_t row_index) const {
    size_t n = columns_.size();

    Matrix<Field> e(n, 1, 0);
    e[row_index, 0] = 1;

    return solve_linear_transposed(e);
  }

  // Returns inverse of the current matrix, reconstructed from LU-decomposition.
  // Note: This method is for testing, it is not optimized in any way.
  Matrix<Field> get_inverse() const {
    auto result = Matrix<Field>::identity(A_.rows());

    result = P_.apply(std::move(result));

    for (auto entry : ls_) {
      result = entry.apply(std::move(result));
    }

    for (auto entry : us_ | std::views::reverse) {
      result = entry.apply(std::move(result));
    }

    result = Q_.apply(std::move(result));

    return result;
  }

  // Returns current matrix, reconstructed from LU-decomposition.
  // Note: This method is for testing, it is not optimized in any way
  Matrix<Field> get_matrix() const {
    auto result = Matrix<Field>::identity(A_.rows());

    result = Q_.apply_transposed(std::move(result));

    for (auto entry : us_) {
      result = entry.apply_inverse(std::move(result));
    }

    for (auto entry : ls_ | std::views::reverse) {
      result = entry.apply_inverse(std::move(result));
    }

    result = P_.apply_transposed(std::move(result));

    return result;
  }

  size_t size() const { return ls_.size() + us_.size(); }
};

}  // namespace linalg
