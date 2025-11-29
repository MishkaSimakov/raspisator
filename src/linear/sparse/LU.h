#pragma once

#include <iostream>
#include <numeric>
#include <set>

#include "CSCMatrix.h"

namespace linalg {

template <typename Field>
CSCMatrix<Field> apply_permutation(const CSCMatrix<Field>& matrix,
                                   const std::vector<size_t>& permutation) {
  auto [n, d] = matrix.shape();

  CSCMatrix<Field> result(n);

  for (size_t col = 0; col < d; ++col) {
    result.add_column();

    for (const auto& [row, value] : matrix.get_column(col)) {
      result.push_to_last_column(permutation[row], value);
    }
  }

  return result;
}

template <typename Field>
CSCMatrix<Field> apply_partial_permutation(
    const CSCMatrix<Field>& matrix, const std::vector<size_t>& permutation) {
  auto [n, d] = matrix.shape();

  size_t valid_count = 0;
  for (size_t i : permutation) {
    if (i < d) {
      ++valid_count;
    }
  }

  CSCMatrix<Field> result(valid_count);

  for (size_t col = 0; col < d; ++col) {
    result.add_column();

    for (const auto& [row, value] : matrix.get_column(col)) {
      if (permutation[row] < d) {
        result.push_to_last_column(permutation[row], value);
      }
    }
  }

  return result;
}

template <typename Field>
void dfs(const CSCMatrix<Field>& L, size_t start,
         const std::vector<size_t>& rows_permutation,
         std::vector<size_t>& nonzero_indices, std::vector<bool>& visited,
         std::vector<size_t>& parent, std::vector<size_t>& child) {
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

// LUP decomposition for sparse matrices.
// Returns L and U as sparse matrices and permutation P such that:
// - P[i] = j means that i-th row in A is j-th row in PA
// - LU = PA
// - L and U does not store zeros
template <typename Field>
std::tuple<CSCMatrix<Field>, CSCMatrix<Field>, std::vector<size_t>> sparse_lu(
    const CSCMatrix<Field>& A) {
  size_t n = A.shape().first;
  if (n != A.shape().second) {
    throw std::invalid_argument("Matrix A must be square.");
  }

  CSCMatrix<Field> L(n);
  CSCMatrix<Field> U(n);

  std::vector<size_t> rows_permutation(n, n);

  std::vector<Field> dense(n, 0);
  std::vector<size_t> nonzero_indices;

  // dfs vectors
  std::vector<bool> visited(n, false);
  std::vector<size_t> parent(n, 0);
  std::vector<size_t> child(n, 0);

  for (size_t j = 0; j < n; ++j) {
    // dense vector representation of computed column

    // TODO: non-recursive
    for (const auto& [index, value] : A.get_column(j)) {
      dense[index] = value;
      dfs(L, index, rows_permutation, nonzero_indices, visited, parent, child);
    }

    assert(std::set(nonzero_indices.begin(), nonzero_indices.end()).size() ==
           nonzero_indices.size());

    // TODO: it seems that this two loops can be merged
    // compute u_j
    Field largest_value;
    size_t largest_value_row;

    for (size_t row : std::views::reverse(nonzero_indices)) {
      if (rows_permutation[row] == n) {
        if (largest_value < FieldTraits<Field>::abs(dense[row])) {
          largest_value = FieldTraits<Field>::abs(dense[row]);
          largest_value_row = row;
        }
      } else {
        for (const auto& [index, value] : L.get_column(rows_permutation[row])) {
          dense[index] -= dense[row] * value;
        }
      }
    }

    rows_permutation[largest_value_row] = j;

    Field diagonal_element;
    // copy dense into appropriate sparse columns of U and L
    U.add_column();
    L.add_column();

    for (size_t row : nonzero_indices) {
      // apply rows permutation
      if (rows_permutation[row] == n) {
        L.push_to_last_column(row, dense[row]);
      } else {
        U.push_to_last_column(row, dense[row]);
      }

      if (rows_permutation[row] == j) {
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

  // apply P^-T to L and U
  // L -> P^-T L, U -> P^-T U
  // since P is a permutation matrix, P^-T = P

  return {linalg::apply_permutation(L, rows_permutation),
          linalg::apply_permutation(U, rows_permutation),
          std::move(rows_permutation)};
}

}  // namespace linalg
