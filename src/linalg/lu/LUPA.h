#pragma once

#include <cmath>
#include <ranges>
#include <unordered_set>
#include <vector>

#include "EtaFile.h"
#include "linalg/lu/FullPivotingLU.h"
#include "linalg/lu/Solve.h"
#include "utils/Accumulators.h"
#include "utils/Logging.h"

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

  // buffers for Forrest-Tomlin update
  Vector<Field> column_;
  Vector<Field> r_;

  // Determinant of B^-1
  Field det_;

  LUPAConfig config_;

  void purge() {
    ls_.purge();
    us_.purge();

    changes_since_purge_ = 0;
  }

  enum class UpdateResult { SUCCESS, NEED_REFACTORIZATION };

  UpdateResult forrest_tomlin_update(size_t current_column, size_t new_column) {
    auto [n, d] = A_.shape();

    // TODO: this computation is duplicated in primal simplex.
    column_ = A_.get_column_as_matrix(new_column);

    column_ = P_.apply(std::move(column_));
    for (auto entry : ls_) {
      column_ = entry.apply(std::move(column_));
    }

    // Permutation Q_ changes columns order. Here we find a column that becomes
    // current_column after applying Q_
    current_column = Q_.post_apply(current_column);

    auto itr = us_.begin();
    for (; itr != us_.end(); ++itr) {
      if ((*itr).pivot_index() == current_column) {
        det_ /= (*itr).det();
        itr = us_.erase(itr);

        break;
      }
    }

    // TODO: it was noticed that r is often empty in setcover problem!
    // Check this for other problems
    r_.resize(n);
    for (size_t i = 0; i < n; ++i) {
      r_[i] = 0;
    }

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

      r_[(*itr).pivot_index()] = main_value / diagonal;
      r_ = (*itr).apply_transposed(std::move(r_));
    }

    r_[current_column] = 1;

    Maximum<Field> r_max;
    for (size_t i = 0; i < n; ++i) {
      r_max.record(r_[i]);
    }

    // if (*r_max > 10) {
    //   std::println("refactorization: {}", *r_max);
    //   return UpdateResult::NEED_REFACTORIZATION;
    // }

    ls_.push_back(current_column, r_, EtaMatrixType::ROW);

    // add new eta matrix to U1 ... Un
    column_ = (*std::prev(ls_.cend())).apply(std::move(column_));

    const Field diagonal = column_[current_column];

    if (!FieldTraits<Field>::is_nonzero(diagonal)) {
      return UpdateResult::NEED_REFACTORIZATION;
    }

    det_ /= diagonal;

    for (size_t i = 0; i < n; ++i) {
      column_[i] =
          i != current_column ? -column_[i] / diagonal : Field(1) / diagonal;
    }

    us_.push_back(current_column, column_, EtaMatrixType::COLUMN);

    return UpdateResult::SUCCESS;
  }

  void guard_columns_set() const {
    if (columns_.size() != A_.rows()) {
      throw std::logic_error("Columns are not set.");
    }
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
    if (columns.size() != A_.rows()) {
      throw std::invalid_argument(
          std::format("Number of selected columns doesn't match the number of "
                      "rows: {} != {}",
                      columns.size(), A_.rows()));
    }

    columns_ = columns;
    refactorize();
  }

  // pivot_element is the \alpha_{pp}, where
  void change_column(size_t current_column, size_t new_column,
                     std::optional<Field> pivot_element = std::nullopt) {
    guard_columns_set();

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
    guard_columns_set();

    factorizer_.get(A_, columns_, P_, Q_, ls_, us_);

    det_ = 1;

    det_ *= P_.is_even() ? 1 : -1;
    det_ *= Q_.is_even() ? 1 : -1;

    for (const auto entry : ls_) {
      det_ *= entry.det();
    }
    for (const auto entry : us_) {
      det_ *= entry.det();
    }

    changes_since_refactorization_ = 0;
    changes_since_purge_ = 0;
  }

  size_t get_changes_since_refactorization() const {
    return changes_since_refactorization_;
  }

  // solves Ax = b
  Vector<Field> solve_linear(Vector<Field> b) const {
    guard_columns_set();

    auto result = linalg::solve_linear(std::move(b), P_, Q_, ls_, us_);
    return result;
  }

  Vector<Field> solve_linear_transposed(Vector<Field> b) const {
    guard_columns_set();

    return linalg::solve_linear_transposed(std::move(b), P_, Q_, ls_, us_);
  }

  Vector<Field> get_row(size_t row_index) const {
    guard_columns_set();

    Vector<Field> e(columns_.size());
    e[row_index] = 1;

    return solve_linear_transposed(std::move(e));
  }

  // Returns inverse of the current matrix, reconstructed from LU-decomposition.
  // Note: This method is for testing, it is not optimized in any way.
  Matrix<Field> get_inverse() const {
    guard_columns_set();

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
    guard_columns_set();

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

  Matrix<Field> get_l() const {
    guard_columns_set();

    auto result = Matrix<Field>::identity(A_.rows());

    for (auto entry : ls_ | std::views::reverse) {
      result = entry.apply_inverse(std::move(result));
    }

    return result;
  }

  Matrix<Field> get_u() const {
    guard_columns_set();

    auto result = Matrix<Field>::identity(A_.rows());

    for (auto entry : us_) {
      result = entry.apply_inverse(std::move(result));
    }

    return result;
  }

  size_t size() const { return ls_.size() + us_.size(); }

  Field det() const {
    guard_columns_set();
    return det_;
  }
};

}  // namespace linalg
