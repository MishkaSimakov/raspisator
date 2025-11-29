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
void dfs(const CSCMatrix<Field>& L, size_t start,
         const std::vector<size_t>& rows_permutation,
         std::vector<size_t>& nonzero_indices, std::vector<bool>& visited,
         std::vector<size_t>& parent, std::vector<size_t>& child) {
  if (visited[start]) {
    return;
  }

  auto [n, d] = L.shape();

  size_t current = start;
  parent[current] = n;
  visited[current] = true;

  while (current != n) {
    if (current >= d) {
      nonzero_indices.push_back(current);
      current = parent[current];
      continue;
    }

    const auto& children = L.get_column(current);

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
    CSCMatrix<Field> A) {
  size_t n = A.shape().first;
  if (n != A.shape().second) {
    throw std::invalid_argument("Matrix A must be square.");
  }

  CSCMatrix<Field> L(n);
  CSCMatrix<Field> U(n);

  std::vector<size_t> rows_permutation(n);
  std::iota(rows_permutation.begin(), rows_permutation.end(), 0);

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

    assert(std::set(nonzero_indices.begin(), nonzero_indices.end()).size() == nonzero_indices.size());

    // TODO: it seems that this two loops can be merged
    // compute u_j
    for (size_t row : std::views::reverse(nonzero_indices)) {
      if (row >= j) {
        continue;
      }

      for (const auto& [index, value] : L.get_column(row)) {
        dense[index] -= dense[row] * value;
      }
    }

    size_t largest_element_row = j;
    for (size_t row : std::views::reverse(nonzero_indices)) {
      if (row >= j) {
        if (FieldTraits<Field>::abs(dense[largest_element_row]) <
            FieldTraits<Field>::abs(dense[row])) {
          largest_element_row = row;
        }
      }
    }

    A.swap_rows(largest_element_row, j);
    std::swap(rows_permutation[j], rows_permutation[largest_element_row]);
    L.swap_rows(largest_element_row, j);

    Field diagonal_element;
    // copy dense into appropriate sparse columns of U and L
    U.add_column();
    L.add_column();

    for (size_t row : nonzero_indices) {
      // apply rows permutation
      size_t permuted_row = row;
      if (row == largest_element_row) {
        permuted_row = j;
      } else if (row == j) {
        permuted_row = largest_element_row;
      }

      if (permuted_row > j) {
        L.push_to_last_column(permuted_row, dense[row]);
      } else {
        U.push_to_last_column(permuted_row, dense[row]);
      }

      if (permuted_row == j) {
        diagonal_element = dense[row];
      }

      dense[row] = 0;
      visited[row] = false;
    }

    for (Field& value : L.get_column(j) | std::views::values) {
      value /= diagonal_element;
    }

    nonzero_indices.clear();

    // std::println("----- j = {} -----", j);
    // std::cout << A << std::endl;
  }

  std::vector<size_t> inverse(rows_permutation.size());
  for (size_t i = 0; i < rows_permutation.size(); ++i) {
    inverse[rows_permutation[i]] = i;
  }

  return {std::move(L), std::move(U), std::move(inverse)};
}

}  // namespace linalg
