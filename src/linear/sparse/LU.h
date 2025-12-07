#pragma once

#include <vector>

#include "CSCMatrix.h"

namespace linalg {

template <typename Field>
void apply_permutation_inplace(CSCMatrix<Field>& matrix,
                               const std::vector<size_t>& permutation) {
  auto [n, d] = matrix.shape();

  for (size_t col = 0; col < d; ++col) {
    for (size_t& row : matrix.get_column(col) | std::views::keys) {
      row = permutation[row];
    }
  }
}

template <typename Field>
CSCMatrix<Field> apply_permutation(const CSCMatrix<Field>& matrix,
                                   const std::vector<size_t>& permutation) {
  CSCMatrix<Field> result = matrix;
  apply_permutation_inplace(result, permutation);

  return result;
}

// LUP-Accelerated (LUPA)
// For given matrix A this class answers queries of form:
// get LUP-decomposition of submatrix of A formed by given columns
template <typename Field>
class LUPA {
  const CSCMatrix<Field>& A_;

  std::vector<Field> dense;
  std::vector<size_t> nonzero_indices;

  // dfs vectors
  std::vector<bool> visited;
  std::vector<size_t> parent;
  std::vector<size_t> child;

  void dfs(const CSCMatrix<Field>& L, size_t start,
           const std::vector<size_t>& rows_permutation) {
    if (visited[start]) {
      return;
    }

    auto [n, d] = L.shape();

    // current is a row index in A numeration
    size_t current = start;
    parent[current] = n;
    visited[current] = true;

    while (current != n) {
      if (rows_permutation[current] == n) {
        nonzero_indices.push_back(current);
        current = parent[current];
        continue;
      }

      const auto& children = L.get_column(rows_permutation[current]);

      // if we visited all children, then exit the current node
      if (child[current] >= children.size()) {
        nonzero_indices.push_back(current);
        child[current] = 0;
        current = parent[current];
        continue;
      }

      size_t next = children[child[current]].first;
      ++child[current];

      if (visited[next]) {
        continue;
      }

      visited[next] = true;
      parent[next] = current;

      current = next;
    }
  }

 public:
  explicit LUPA(const CSCMatrix<Field>& A)
      : A_(A),
        dense(A.shape().first, 0),
        visited(A.shape().first, false),
        parent(A.shape().first, 0),
        child(A.shape().first, 0) {
    nonzero_indices.reserve(A.shape().first);
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
  void get_lup(const std::vector<size_t>& columns, CSCMatrix<Field>& L,
               CSCMatrix<Field>& U, std::vector<size_t>& P) {
    size_t n = A_.shape().first;

    assert(L.shape().first == n);
    assert(U.shape().first == n);
    assert(P.size() >= n);
    assert(columns.size() == n);

    L.clear();
    U.clear();
    std::ranges::fill_n(P.begin(), n, n);

    for (size_t j = 0; j < n; ++j) {
      for (const auto& [index, value] : A_.get_column(columns[j])) {
        dense[index] = value;
        dfs(L, index, P);
      }

      Field largest_value = 0;
      size_t largest_value_row = n;

      for (size_t row : std::views::reverse(nonzero_indices)) {
        if (P[row] == n) {
          if (largest_value < FieldTraits<Field>::abs(dense[row])) {
            largest_value = FieldTraits<Field>::abs(dense[row]);
            largest_value_row = row;
          }
        } else {
          for (const auto& [index, value] : L.get_column(P[row])) {
            dense[index] -= dense[row] * value;
          }
        }
      }

      assert(largest_value_row != n);

      P[largest_value_row] = j;

      Field diagonal_element;
      // copy dense into appropriate sparse columns of U and L
      U.add_column();
      L.add_column();

      for (size_t row : nonzero_indices) {
        // apply rows permutation
        if (P[row] == n) {
          L.push_to_last_column(row, dense[row]);
        } else {
          U.push_to_last_column(row, dense[row]);
        }

        if (P[row] == j) {
          diagonal_element = dense[row];
        }

        dense[row] = 0;
        visited[row] = false;
      }

      for (Field& value : L.get_column(j) | std::views::values) {
        value /= diagonal_element;
      }

      nonzero_indices.clear();
    }

    linalg::apply_permutation_inplace(L, P);
    linalg::apply_permutation_inplace(U, P);
  }
};

template <typename Field>
std::tuple<CSCMatrix<Field>, CSCMatrix<Field>, std::vector<size_t>> sparse_lup(
    const CSCMatrix<Field>& A, const std::vector<size_t>& columns) {
  size_t n = A.shape().first;

  CSCMatrix<Field> L(n);
  CSCMatrix<Field> U(n);
  std::vector<size_t> rows_permutation(n);

  LUPA<Field>(A).get_lup(columns, L, U, rows_permutation);

  return {std::move(L), std::move(U), std::move(rows_permutation)};
}


// solves Ax = b, where PA = LU
// L must be without ones on the main diagonal
template <typename Field>
Matrix<Field> solve_linear(const CSCMatrix<Field>& L, const CSCMatrix<Field>& U,
                           const std::vector<size_t>& P,
                           const Matrix<Field>& b) {
  auto [n, _] = L.shape();

  // construct Pb
  Matrix<Field> result(n, 1);
  for (size_t i = 0; i < n; ++i) {
    result[P[i], 0] = b[i, 0];
  }

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
                                     const std::vector<size_t>& P,
                                     Matrix<Field>& b) {
  auto [n, _] = L.shape();

  // solve U^T y = b
  for (size_t col = 0; col < n; ++col) {
    Field diagonal;

    for (const auto& [row, value] : U.get_column(col)) {
      if (row != col) {
        b[col, 0] -= value * b[row, 0];
      } else {
        diagonal = value;
      }
    }

    b[col, 0] /= diagonal;
  }

  // solve L^T x = y
  for (size_t i = 0; i < n; ++i) {
    size_t col = n - i - 1;

    for (const auto& [row, value] : L.get_column(col)) {
      b[col, 0] -= b[row, 0] * value;
    }
  }
}

// same as solve_transposed_linear_inplace, but returns solution in a newly
// allocated memory
template <typename Field>
Matrix<Field> solve_transposed_linear(const CSCMatrix<Field>& L,
                                      const CSCMatrix<Field>& U,
                                      const std::vector<size_t>& P,
                                      const Matrix<Field>& b) {
  auto copy = b;
  solve_transposed_linear_inplace(L, U, P, copy);
  return copy;
}

}  // namespace linalg
