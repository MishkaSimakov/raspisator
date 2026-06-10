#pragma once

#include <cassert>
#include <vector>

#include "linalg/CSCMatrix.h"
#include "linalg/Matrix.h"
#include "linalg/Permutation.h"
#include "linalg/lu/EtaFile.h"
#include "linalg/lu/SingularityError.h"
#include "utils/Accumulators.h"

namespace linalg {

// This class performs LU decomposition with full pivoting for a sparse matrix.
// Suppose it received a matrix A. It returns P, Q, ls and us such that:
// 1. P is row permutation of A
// 2. Q is column permutation of A
// 3. us is eta-file containing L1^-1 ... Ln^-1
// 4. ls is eta-file containing U1^-1 ... Un^-1, where Li^-1 and Ui^-1 are
// column eta-matrices with pivot column i and L = L1 ... Ln, U = Un ... U1
// 5. PAQ = LU
template <typename Field>
class FullPivotingLU {
  size_t size_;

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

 public:
  explicit FullPivotingLU(size_t size)
      : size_(size),
        dense_(size, 0),
        L_(CSCMatrix<Field>::zeros(size)),
        U_(CSCMatrix<Field>::zeros(size)),
        P_impl_(size),
        Q_impl_(size),
        visited_(size, false),
        parent_(size, 0),
        child_(size, 0) {
    nonzero_indices_.reserve(size);
  }

  std::tuple<Permutation, Permutation, EtaFile<Field>, EtaFile<Field>> get(
      const CSCMatrix<Field>& A, const std::vector<size_t>& columns) {
    auto P = Permutation::id(size_);
    auto Q = Permutation::id(size_);

    EtaFile<Field> ls(size_);
    EtaFile<Field> us(size_);

    get(A, columns, P, Q, ls, us);
    return {std::move(P), std::move(Q), std::move(ls), std::move(us)};
  }

  void get(const CSCMatrix<Field>& A, const std::vector<size_t>& columns,
           Permutation& P, Permutation& Q, EtaFile<Field>& ls,
           EtaFile<Field>& us) {
    using std::abs;

    size_t n = size_;

    assert(A.shape().first == n && columns.size() == n);
    assert(L_.shape().first == n);
    assert(U_.shape().first == n);
    assert(P_impl_.size() >= n);
    assert(Q_impl_.size() >= n);
    assert(columns.size() == n);

    std::vector<size_t> rows_nonzeros(n, 0);
    for (size_t i = 0; i < n; ++i) {
      for (auto [row, value] : A.get_column(columns[i])) {
        ++rows_nonzeros[row];
      }
    }

    L_.clear();
    U_.clear();
    std::ranges::fill_n(P_impl_.begin(), n, n);
    std::ranges::fill_n(Q_impl_.begin(), n, n);

    for (size_t j = 0; j < n; ++j) {
      // choose pivot column with the least amount of elements
      ArgMinimum<size_t, std::less<>> min_nz_column;

      for (size_t i = 0; i < n; ++i) {
        if (Q_impl_[i] == n) {
          min_nz_column.record(i, A.get_column(columns[i]).size());
        }
      }

      assert(min_nz_column.has_value());
      const size_t pivot_column = min_nz_column->index;

      for (const auto& [index, value] : A.get_column(columns[pivot_column])) {
        dense_[index] = value;
        dfs(L_, index, P_impl_);
      }

      ArgMaximum<Field> max_value;

      for (const size_t row : std::views::reverse(nonzero_indices_)) {
        if (P_impl_[row] == n) {
          max_value.record(row, abs(dense_[row]));
        } else {
          for (const auto& [index, value] : L_.get_column(P_impl_[row])) {
            dense_[index] -= dense_[row] * value;
          }
        }
      }

      if (!max_value.has_value() ||
          !FieldTraits<Field>::is_nonzero(max_value->max)) {
        throw SingularityError();
      }

      // choose pivoting row
      Field threshold = 0.75;
      ArgMinimum<size_t, std::less<>> min_nz_row;

      for (size_t row : std::views::reverse(nonzero_indices_)) {
        if (P_impl_[row] == n &&
            abs(dense_[row]) > threshold * max_value->max) {
          min_nz_row.record(row, rows_nonzeros[row]);
        }
      }

      if (!min_nz_row.has_value()) {
        throw SingularityError();
      }

      const size_t pivot_row = min_nz_row->index;

      P_impl_[pivot_row] = j;
      Q_impl_[pivot_column] = j;

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

      for (auto& [_, value] : L_.get_column(j)) {
        value /= diagonal_element;
      }

      nonzero_indices_.clear();
    }

    P = Permutation::from_vector(P_impl_);

    std::vector<size_t> Q_transposed(n);
    for (size_t i = 0; i < n; ++i) {
      Q_transposed[Q_impl_[i]] = i;
    }
    Q = Permutation::from_vector(std::move(Q_transposed));
    L_ = P.apply(std::move(L_));
    U_ = P.apply(std::move(U_));

    // calculate eta files
    us.clear();
    ls.clear();

    Maximum<Field> max_u;
    Maximum<Field> max_l;
    size_t nonzeros = 0;

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

        max_u.record(abs(value));
        ++nonzeros;
      }

      us.push_back(i, column);

      // lower
      column.clear();

      column.emplace_back(i, 1);
      for (auto [row, value] : L_.get_column(i)) {
        column.emplace_back(row, -value);

        max_l.record(abs(value));
        ++nonzeros;
      }

      ls.push_back(i, column);
    }

    // logging::log_value(*max_u.max(), "max_u_value.txt");
    // logging::log_value(*max_l.max(), "max_l_value.txt");
    // logging::log_value(nonzeros, "lu_nonzeros_count.txt");
  }
};

}  // namespace linalg
