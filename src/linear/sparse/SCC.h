#pragma once

#include <cassert>

#include "CSCMatrix.h"

// Block triangular form is very useful for LU decomposition. It is closely
// related with strongly connected components (SCC) of a sparsity graph S(A).

namespace linalg {

template <typename Field>
class TarjanSCC {
  size_t size_;

  size_t index_;
  std::vector<bool> visited_;
  std::vector<bool> in_stack_;
  std::vector<size_t> low_;
  std::vector<size_t> id_;

  std::vector<size_t> stack_;

  void dfs(size_t node, const CSCMatrix<Field>& matrix,
           const std::vector<size_t>& columns, std::vector<size_t>& roots) {
    visited_[node] = true;
    id_[node] = index_;
    low_[node] = index_;

    ++index_;

    stack_.push_back(node);
    in_stack_[node] = true;

    for (auto [row, _] : matrix.get_column(columns[node])) {
      if (!visited_[row]) {
        dfs(row, matrix, columns, roots);
        low_[node] = std::min(low_[node], low_[row]);
      } else if (in_stack_[node]) {
        low_[node] = std::min(low_[node], id_[row]);
      }
    }

    if (low_[node] == id_[node]) {
      while (true) {
        size_t u = stack_.back();

        in_stack_[u] = false;
        stack_.pop_back();

        roots[u] = node;

        if (u == node) {
          break;
        }
      }
    }
  }

 public:
  explicit TarjanSCC(size_t n)
      : size_(n),
        index_(0),
        visited_(n, false),
        in_stack_(n, false),
        low_(n),
        id_(n) {}

  std::vector<size_t> get(const CSCMatrix<Field>& matrix,
                          const std::vector<size_t>& columns) {
    std::vector<size_t> roots(size_);

    get(matrix, columns, roots);

    return roots;
  }

  void get(const CSCMatrix<Field>& matrix, const std::vector<size_t>& columns,
           std::vector<size_t>& roots) {
    assert(matrix.shape().first == size_ && columns.size() == size_);
    assert(roots.size() == size_);

    for (size_t i = 0; i < size_; ++i) {
      if (!visited_[i]) {
        dfs(i, matrix, columns, roots);
      }
    }

    // reset state
    std::fill_n(visited_.begin(), size_, false);
    index_ = 0;
  }
};

}  // namespace linalg
