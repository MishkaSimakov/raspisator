#pragma once

#include <ranges>
#include <vector>

#include "CSCMatrix.h"
#include "EtaFile.h"
#include "Permutation.h"
#include "linear/matrix/LU.h"
#include "utils/Accumulators.h"
#include "utils/Logging.h"

namespace linalg {

// solves Ax = b, where PA = LU
// L must be without ones on the main diagonal
template <typename Field>
Matrix<Field> solve_linear(const CSCMatrix<Field>& L, const CSCMatrix<Field>& U,
                           const Permutation& P, const Matrix<Field>& b) {
  auto [n, _] = L.shape();

  Matrix<Field> result = P.apply(b);

  // solve Ly = Pb
  for (size_t col = 0; col < n; ++col) {
    for (const auto& [row, value] : L.get_column(col)) {
      result[row, 0] -= result[col, 0] * value;
    }
  }

  // solve Ux = y
  for (size_t i = 0; i < n; ++i) {
    size_t col = n - i - 1;

    for (const auto& [row, value] : U.get_column(col)) {
      if (row == col) {
        result[col, 0] /= value;
        break;
      }
    }

    for (const auto& [row, value] : U.get_column(col)) {
      if (row != col) {
        result[row, 0] -= result[col, 0] * value;
      }
    }
  }

  return result;
}

// solves A^T x = b, where PA = LU
// uses memory of b as output
// L must be without ones on the main diagonal
// returns a point y, such that x[i] = y[P[i]]
template <typename Field>
void solve_transposed_linear_inplace(const CSCMatrix<Field>& L,
                                     const CSCMatrix<Field>& U,
                                     const Permutation& P, Matrix<Field>& b) {
  auto [n, _] = L.shape();

  // solve U^T y = b
  for (size_t col = 0; col < n; ++col) {
    Field diagonal;

    // std::cout << b[col, 0] << " -= ";

    for (const auto& [row, value] : U.get_column(col)) {
      if (row != col) {
        b[col, 0] -= value * b[row, 0];

        // std::cout << value << "*" << b[row, 0] << " + ";
      } else {
        diagonal = value;
      }
    }

    // std::cout << " | " << b[col, 0] << " " << diagonal << "\n";
    b[col, 0] /= diagonal;

    if (!FieldTraits<Field>::is_nonzero(b[col, 0])) {
      b[col, 0] = 0;
    }
  }

  // solve L^T x = y
  for (size_t i = 0; i < n; ++i) {
    size_t col = n - i - 1;

    for (const auto& [row, value] : L.get_column(col)) {
      b[col, 0] -= b[row, 0] * value;
    }

    if (!FieldTraits<Field>::is_nonzero(b[col, 0])) {
      b[col, 0] = 0;
    }
  }
}

// same as solve_transposed_linear_inplace, but returns solution in a newly
// allocated memory
template <typename Field>
Matrix<Field> solve_transposed_linear(const CSCMatrix<Field>& L,
                                      const CSCMatrix<Field>& U,
                                      const Permutation& P,
                                      const Matrix<Field>& b) {
  auto copy = b;
  solve_transposed_linear_inplace(L, U, P, copy);
  return copy;
}

// LUP-Accelerated (LUPA)
// For a given matrix this class answers queries of form:
// get LUP-decomposition of submatrix of A formed by given columns
template <typename Field>
class LUPA {
  const CSCMatrix<Field>& A_;
  std::vector<size_t> columns_;

  // decomposition
  EtaFile<Field> us_;  // stores U1^-1 ... U(n+t)^-1
  EtaFile<Field> ls_;  // stores L1^-1 ... Ln^-1 R1^-1 ... Rt^-1
  Permutation P_;
  Permutation Q_;

  // Forrest-Tomlin update helpers
  size_t changes_since_refactorization_{0};
  size_t changes_since_purge_{0};
  bool force_refactorization_{false};

  // LU-decomposition implementation buffers
  std::vector<Field> dense_;
  std::vector<size_t> nonzero_indices_;

  CSCMatrix<Field> L_;
  CSCMatrix<Field> U_;
  std::vector<size_t> P_impl_;
  std::vector<size_t> Q_impl_;

  std::vector<bool> visited_;
  std::vector<size_t> parent_;
  std::vector<size_t> child_;

  void dfs(const CSCMatrix<Field>& L, size_t start,
           const std::vector<size_t>& rows_permutation) {
    if (visited_[start]) {
      return;
    }

    auto [n, d] = L.shape();

    // current is a row index in A numeration
    size_t current = start;
    parent_[current] = n;
    visited_[current] = true;

    while (current != n) {
      if (rows_permutation[current] == n) {
        nonzero_indices_.push_back(current);
        current = parent_[current];
        continue;
      }

      const auto& children = L.get_column(rows_permutation[current]);

      // if we visited all children, then exit the current node
      if (child_[current] >= children.size()) {
        nonzero_indices_.push_back(current);
        child_[current] = 0;
        current = parent_[current];
        continue;
      }

      size_t next = children[child_[current]].first;
      ++child_[current];

      if (visited_[next]) {
        continue;
      }

      visited_[next] = true;
      parent_[next] = current;

      current = next;
    }
  }

  // LUP decomposition for given submatrix of A.
  // First it selects only given columns from A.
  // Returns L and U as sparse matrices and permutation P such that:
  // - P[i] = j means that i-th row in A is j-th row in PA
  // - LU = PA
  // - L and U does not store zeros
  // - L does not store ones on the main diagonal
  // Note: It is not checked whether A is non-singular, but A must be
  // non-singular for algorithm to work.
  void get_lup(const std::vector<size_t>& columns) {
    size_t n = A_.shape().first;

    assert(L_.shape().first == n);
    assert(U_.shape().first == n);
    assert(P_impl_.size() >= n);
    assert(columns.size() == n);

    L_.clear();
    U_.clear();
    std::ranges::fill_n(P_impl_.begin(), n, n);

    for (size_t j = 0; j < n; ++j) {
      for (const auto& [index, value] : A_.get_column(columns[j])) {
        dense_[index] = value;
        dfs(L_, index, P_impl_);
      }

      Field largest_value = 0;
      size_t largest_value_row = n;

      for (size_t row : std::views::reverse(nonzero_indices_)) {
        if (P_impl_[row] == n) {
          if (largest_value < FieldTraits<Field>::abs(dense_[row])) {
            largest_value = FieldTraits<Field>::abs(dense_[row]);
            largest_value_row = row;
          }
        } else {
          for (const auto& [index, value] : L_.get_column(P_impl_[row])) {
            dense_[index] -= dense_[row] * value;
          }
        }
      }

      // logging::log_value(largest_value, "lu_pivoting_value.txt");
      assert(largest_value_row != n);

      P_impl_[largest_value_row] = j;

      Field diagonal_element;
      // copy dense into appropriate sparse columns of U and L
      U_.add_column();
      L_.add_column();

      for (size_t row : nonzero_indices_) {
        // apply rows permutation
        if (P_impl_[row] == n) {
          L_.push_to_last_column(row, dense_[row]);
        } else {
          U_.push_to_last_column(row, dense_[row]);
        }

        if (P_impl_[row] == j) {
          diagonal_element = dense_[row];
        }

        dense_[row] = 0;
        visited_[row] = false;
      }

      for (Field& value : L_.get_column(j) | std::views::values) {
        value /= diagonal_element;
      }

      nonzero_indices_.clear();
    }

    P_ = Permutation::from_vector(P_impl_);
    L_ = P_.apply(std::move(L_));
    U_ = P_.apply(std::move(U_));

    // calculate eta files
    us_.clear();
    ls_.clear();

    Maximum<Field> max;

    for (size_t i = 0; i < n; ++i) {
      // upper
      std::vector<std::pair<size_t, Field>> column;
      Field diagonal;

      for (auto [row, value] : U_.get_column(i)) {
        if (row == i) {
          diagonal = value;
          column.emplace_back(row, 1);
        } else {
          column.emplace_back(row, -value);
        }
      }

      for (auto& [row, value] : column) {
        value /= diagonal;
        max.record(FieldTraits<Field>::abs(value));
      }

      us_.push_back(i, column);

      // lower
      column.clear();

      column.emplace_back(i, 1);
      for (auto [row, value] : L_.get_column(i)) {
        column.emplace_back(row, -value);
      }

      ls_.push_back(i, column);
    }

    // logging::log_value(*max.max(), "max_us_value.txt");
  }

  void refactorize() {
    get_lup(columns_);

    changes_since_refactorization_ = 0;
    changes_since_purge_ = 0;
    force_refactorization_ = false;
  }

  void purge() {
    ls_.purge();
    us_.purge();

    changes_since_purge_ = 0;
  }

  void forrest_tomlin_update(size_t current_column, size_t new_column) {
    auto [n, d] = A_.shape();

    Matrix<Field> column(n, 1, 0);
    for (const auto [row, value] : A_.get_column(new_column)) {
      column[row, 0] = value;
    }

    column = P_.apply(std::move(column));
    for (auto entry : ls_) {
      column = ls_.apply(std::move(column), entry);
    }

    Matrix<Field> r(n, 1, 0);

    auto itr = us_.begin();
    for (; itr != us_.end(); ++itr) {
      if ((*itr).index == current_column) {
        itr = us_.erase(itr);
        break;
      }
    }

    // TODO: it was noticed that r is often empty in setcover problem!
    // Check this for other problems

    for (; itr != us_.end(); ++itr) {
      Field diagonal = 0;
      Field main_value = 0;

      for (auto& [row, value] : (*itr).values) {
        if (row == current_column) {
          main_value = value;
          value = 0;
        } else if (row == (*itr).index) {
          diagonal = value;
        }
      }

      assert(FieldTraits<Field>::is_nonzero(diagonal));

      r[(*itr).index, 0] = main_value / diagonal;
      r = us_.apply_transposed(std::move(r), *itr);
    }

    r[current_column, 0] = 1;

    Maximum<Field> r_max;
    for (size_t i = 0; i < n; ++i) {
      r_max.record(r[i, 0]);
    }

    if (*r_max.max() > 20) {
      force_refactorization_ = true;
    }

    // logging::log_value(*r_max.max(), "r_max.txt");

    ls_.push_back(current_column, r, EtaType::ROW);

    // add new eta matrix to U1 ... Un
    column = ls_.apply(std::move(column), *(--ls_.cend()));
    assert(FieldTraits<Field>::is_nonzero(column[current_column, 0]));

    Field diagonal = column[current_column, 0];

    Maximum<Field> max;

    for (size_t i = 0; i < n; ++i) {
      column[i, 0] =
          i != current_column ? -column[i, 0] / diagonal : Field(1) / diagonal;

      max.record(FieldTraits<Field>::abs(column[i, 0]));
    }

    // logging::log_value(*max.max(), "max_additional_us_value.txt");

    us_.push_back(current_column, column, EtaType::COLUMN);
  }

 public:
  explicit LUPA(const CSCMatrix<Field>& A)
      : A_(A),
        P_(Permutation::id(A.shape().first)),
        Q_(Permutation::id(A.shape().first)),
        dense_(A.shape().first, 0),
        L_(A_.shape().first),
        U_(A_.shape().first),
        P_impl_(A.shape().first),
        Q_impl_(A.shape().first),
        visited_(A.shape().first, false),
        parent_(A.shape().first, 0),
        child_(A.shape().first, 0) {
    nonzero_indices_.reserve(A.shape().first);
  }

  void set_columns(const std::vector<size_t>& columns) {
    columns_ = columns;
    get_lup(columns_);
  }

  void change_column(size_t current_column, size_t new_column) {
    columns_[current_column] = new_column;
    ++changes_since_refactorization_;
    ++changes_since_purge_;

    if (changes_since_refactorization_ > 500 || force_refactorization_) {
      refactorize();
    } else {
      forrest_tomlin_update(current_column, new_column);

      if (changes_since_purge_ > 100) {
        purge();
      }
    }
  }

  // solves Ax = b
  Matrix<Field> solve_linear(Matrix<Field> b) const {
    Matrix<Field> result = P_.apply(std::move(b));

    for (auto entry : ls_) {
      result = ls_.apply(std::move(result), entry);
    }
    for (auto entry : us_ | std::views::reverse) {
      result = us_.apply(std::move(result), entry);
    }

    return result;
  }

  Matrix<Field> solve_linear_transposed(Matrix<Field> b) const {
    for (auto entry : us_) {
      b = us_.apply_transposed(std::move(b), entry);
    }
    for (auto entry : ls_ | std::views::reverse) {
      b = ls_.apply_transposed(std::move(b), entry);
    }

    return P_.apply_transposed(std::move(b));
  }

  Matrix<Field> get_row(size_t row_index) const {
    size_t n = columns_.size();

    Matrix<Field> e(n, 1, 0);
    e[row_index, 0] = 1;

    return solve_linear_transposed(e);
  }

  auto get_lup() const { return std::tie(L_, U_, P_); }

  // This method is for testing, it is not optimized in any way
  Matrix<Field> get_inverse() const {
    auto result = Matrix<Field>::unity(A_.shape().first);

    result = P_.apply(std::move(result));

    for (auto entry : ls_) {
      result = ls_.apply(std::move(result), entry);
    }

    for (auto entry : us_ | std::views::reverse) {
      result = us_.apply(std::move(result), entry);
    }

    return result;
  }

  size_t size() const { return ls_.size() + us_.size(); }
};

template <typename Field>
std::tuple<CSCMatrix<Field>, CSCMatrix<Field>, Permutation> sparse_lup(
    const CSCMatrix<Field>& A, const std::vector<size_t>& columns) {
  auto lupa = LUPA<Field>(A);
  lupa.set_columns(columns);

  return lupa.get_lup();
}

}  // namespace linalg
