#pragma once

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "linalg/Linalg.h"
#include "linear/model/Bound.h"

namespace problem {

// Core LP problem representation. It is shared by LP and StandardLP.
// Note: Objective is to MAXIMIZE cost.
template <typename Field>
struct CoreLP {
  CSCMatrix<Field> matrix;
  Vector<Field> cost;

  Field cost_offset{0};

  std::vector<Bound<Field>> var_bounds;

  std::string name;
  std::string cost_name;

  std::vector<std::string> var_names;
  std::vector<std::string> row_names;

  CoreLP() = default;

  // Construct empty but valid problem with @rows constraints and @cols
  // variables. All names are empty.
  CoreLP(size_t rows, size_t cols)
      : matrix(CSCMatrix<Field>::zeros(rows, cols)),
        cost(cols),
        var_bounds(cols),
        var_names(cols),
        row_names(rows) {}

  // for debugging purposes, throws if the problem is incorrect
  void validate() const {
    const auto [n, d] = matrix.shape();

    if (cost.size() != d) {
      throw std::runtime_error("Wrong cost vector size.");
    }

    if (var_bounds.size() != d) {
      throw std::runtime_error("Wrong variable bounds vector size.");
    }

    if (var_names.size() != d) {
      throw std::runtime_error("Wrong variable names vector size.");
    }

    if (row_names.size() != n) {
      throw std::runtime_error("Wrong row names vector size.");
    }
  }

  std::string row_name(size_t index) const {
    return row_names[index].empty() ? "r" + std::to_string(index)
                                    : row_names[index];
  }

  std::string var_name(size_t index) const {
    return var_names[index].empty() ? "x" + std::to_string(index)
                                    : var_names[index];
  }
};

}  // namespace problem
