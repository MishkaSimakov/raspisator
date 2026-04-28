#pragma once

#include <cmath>
#include <ranges>
#include <unordered_set>
#include <vector>

#include "CSCMatrix.h"
#include "EtaFile.h"
#include "Permutation.h"
#include "SCC.h"
#include "linear/matrix/LU.h"
#include "linear/matrix/Norms.h"
#include "utils/Accumulators.h"
#include "utils/Logging.h"

namespace linalg {

struct SingularityError final : std::runtime_error {
  SingularityError()
      : std::runtime_error("Matrix is possibly singular in LU decomposition") {}
};

template <typename Field>
Field scale_factor(const Matrix<Field>&) {
  return 1;
}

template <>
inline double scale_factor(const Matrix<double>& b) {
  auto [n, d] = b.shape();

  Maximum<double> max;
  Minimum<double> min;

  for (size_t i = 0; i < d; ++i) {
    for (size_t j = 0; j < n; ++j) {
      max.record(FieldTraits<double>::abs(b[j, i]));

      if (FieldTraits<double>::is_nonzero(b[j, i])) {
        min.record(FieldTraits<double>::abs(b[j, i]));
      }
    }
  }

  if (!min.min() || !max.max()) {
    return 1;
  }

  double scale_factor = std::sqrt(*max.max() * *min.min());

  // round to power of 2
  int power = std::round(std::log2(scale_factor));

  return std::exp2(-power);
}

// solves Ax = b, where PAQ = LU
// ls is eta-file containing L1^-1 ... Ln^-1, where L = L1 ... Ln
// us is eta-file containing U1^-1 ... Un^-1, where U = Un ... U1
template <typename Field>
Matrix<Field> solve_linear(Matrix<Field> b, const Permutation& P,
                           const Permutation& Q, const EtaFile<Field>& ls,
                           const EtaFile<Field> us) {
  Matrix<Field> result = P.apply(std::move(b));

  auto sf = scale_factor(result);
  result /= sf;

  for (auto entry : ls) {
    result = ls.apply(std::move(result), entry);
  }

  for (auto entry : us | std::views::reverse) {
    result = us.apply(std::move(result), entry);
  }

  result = Q.apply(std::move(result));
  result *= sf;

  return result;
}

// solves A^T x = b, where PAQ = LU
// ls is eta-file containing L1^-1 ... Ln^-1, where L = L1 ... Ln
// us is eta-file containing U1^-1 ... Un^-1, where U = Un ... U1
template <typename Field>
Matrix<Field> solve_linear_transposed(Matrix<Field> b, const Permutation& P,
                                      const Permutation& Q,
                                      const EtaFile<Field>& ls,
                                      const EtaFile<Field>& us) {
  b = Q.apply_transposed(std::move(b));

  auto sf = scale_factor(b);
  b /= sf;

  for (auto entry : us) {
    b = us.apply_transposed(std::move(b), entry);
  }

  for (auto entry : ls | std::views::reverse) {
    b = ls.apply_transposed(std::move(b), entry);
  }

  b = P.apply_transposed(std::move(b));
  b *= sf;

  return b;
}

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
        L_(size),
        U_(size),
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

    EtaFile<Field> ls;
    EtaFile<Field> us;

    get(A, columns, P, Q, ls, us);
    return {std::move(P), std::move(Q), std::move(ls), std::move(us)};
  }

  void get(const CSCMatrix<Field>& A, const std::vector<size_t>& columns,
           Permutation& P, Permutation& Q, EtaFile<Field>& ls,
           EtaFile<Field>& us) {
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
      ArgMinimum<size_t, std::less<>> nz_count;

      for (size_t i = 0; i < n; ++i) {
        if (Q_impl_[i] == n) {
          nz_count.record(i, A.get_column(columns[i]).size());
        }
      }

      const size_t pivot_column = *nz_count.argmin();

      for (const auto& [index, value] : A.get_column(columns[pivot_column])) {
        dense_[index] = value;
        dfs(L_, index, P_impl_);
      }

      ArgMaximum<Field> max_value;

      for (const size_t row : std::views::reverse(nonzero_indices_)) {
        if (P_impl_[row] == n) {
          max_value.record(row, FieldTraits<Field>::abs(dense_[row]));
        } else {
          for (const auto& [index, value] : L_.get_column(P_impl_[row])) {
            dense_[index] -= dense_[row] * value;
          }
        }
      }

      if (!max_value.max().has_value() ||
          !FieldTraits<Field>::is_nonzero(*max_value.max())) {
        throw SingularityError();
      }

      // choose pivoting row
      Field threshold = 0.75;
      ArgMinimum<size_t, std::less<>> row_nz_count;

      for (size_t row : std::views::reverse(nonzero_indices_)) {
        if (P_impl_[row] == n && FieldTraits<Field>::abs(dense_[row]) >
                                     threshold * *max_value.max()) {
          row_nz_count.record(row, rows_nonzeros[row]);
        }
      }

      if (!row_nz_count.argmin().has_value()) {
        throw SingularityError();
      }

      size_t pivot_row = *row_nz_count.argmin();

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

        max_u.record(FieldTraits<Field>::abs(value));
        ++nonzeros;
      }

      us.push_back(i, column);

      // lower
      column.clear();

      column.emplace_back(i, 1);
      for (auto [row, value] : L_.get_column(i)) {
        column.emplace_back(row, -value);

        max_l.record(FieldTraits<Field>::abs(value));
        ++nonzeros;
      }

      ls.push_back(i, column);
    }

    // logging::log_value(*max_u.max(), "max_u_value.txt");
    // logging::log_value(*max_l.max(), "max_l_value.txt");
    // logging::log_value(nonzeros, "lu_nonzeros_count.txt");
  }
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
      column = ls_.apply(std::move(column), entry);
    }

    // Permutation Q_ changes columns order. Here we find a column that becomes
    // current_column after applying Q_
    current_column = Q_.post_apply(current_column);

    Matrix<Field> r(n, 1, 0);

    auto itr = us_.begin();
    for (; itr != us_.end(); ++itr) {
      if ((*itr).index == current_column) {
        for (auto [row, value] : (*itr).values) {
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

      for (auto& [row, value] : (*itr).values) {
        if (row == current_column) {
          main_value = value;
          value = 0;
        } else if (row == (*itr).index) {
          diagonal = value;
        }
      }

      if (!FieldTraits<Field>::is_nonzero(diagonal)) {
        throw SingularityError();
      }

      r[(*itr).index, 0] = main_value / diagonal;
      r = us_.apply_transposed(std::move(r), *itr);
    }

    r[current_column, 0] = 1;

    Maximum<Field> r_max;
    for (size_t i = 0; i < n; ++i) {
      r_max.record(r[i, 0]);
    }

    if (*r_max.max() > 10) {
      // std::println("refactorization: {}", *r_max.max());
      return UpdateResult::NEED_REFACTORIZATION;
    }

    // logging::log_value(*r_max.max(), "r_max.txt");

    ls_.push_back(current_column, r, EtaType::ROW);

    // add new eta matrix to U1 ... Un
    column = ls_.apply(std::move(column), *(--ls_.cend()));

    if (!FieldTraits<Field>::is_nonzero(column[current_column, 0])) {
      throw SingularityError();
    }

    Field diagonal = column[current_column, 0];

    for (size_t i = 0; i < n; ++i) {
      column[i, 0] =
          i != current_column ? -column[i, 0] / diagonal : Field(1) / diagonal;
    }

    us_.push_back(current_column, column, EtaType::COLUMN);

    return UpdateResult::SUCCESS;
  }

 public:
  explicit LUPA(const CSCMatrix<Field>& A)
      : A_(A),
        factorizer_(A.shape().first),
        P_(Permutation::id(A.shape().first)),
        Q_(Permutation::id(A.shape().first)) {}

  void set_columns(const std::vector<size_t>& columns) {
    assert(columns.size() == A_.shape().first);

    columns_ = columns;
    refactorize();
  }

  void change_column(size_t current_column, size_t new_column) {
    columns_[current_column] = new_column;
    ++changes_since_refactorization_;
    ++changes_since_purge_;

    if (changes_since_refactorization_ > 500) {
      refactorize();
      return;
    }

    auto update_result = forrest_tomlin_update(current_column, new_column);

    if (update_result == UpdateResult::NEED_REFACTORIZATION) {
      refactorize();
      return;
    }

    if (changes_since_purge_ > 100) {
      purge();
    }
  }

  void refactorize() {
    factorizer_.get(A_, columns_, P_, Q_, ls_, us_);

    changes_since_refactorization_ = 0;
    changes_since_purge_ = 0;
  }

  // solves Ax = b
  Matrix<Field> solve_linear(Matrix<Field> b) const {
    return linalg::solve_linear(std::move(b), P_, Q_, ls_, us_);
  }

  Matrix<Field> solve_linear_transposed(Matrix<Field> b) const {
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
    auto result = Matrix<Field>::unity(A_.shape().first);

    result = P_.apply(std::move(result));

    for (auto entry : ls_) {
      result = ls_.apply(std::move(result), entry);
    }

    for (auto entry : us_ | std::views::reverse) {
      result = us_.apply(std::move(result), entry);
    }

    result = Q_.apply(std::move(result));

    return result;
  }

  // Returns current matrix, reconstructed from LU-decomposition.
  // Note: This method is for testing, it is not optimized in any way
  Matrix<Field> get_matrix() const {
    auto result = Matrix<Field>::unity(A_.shape().first);

    result = Q_.apply_transposed(std::move(result));

    for (auto entry : us_) {
      result = us_.apply_inverse(std::move(result), entry);
    }

    for (auto entry : ls_ | std::views::reverse) {
      result = ls_.apply_inverse(std::move(result), entry);
    }

    result = P_.apply_transposed(std::move(result));

    return result;
  }

  size_t size() const { return ls_.size() + us_.size(); }
};

}  // namespace linalg
