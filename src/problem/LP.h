#pragma once

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "linear/model/Bound.h"
#include "linear/sparse/CSCMatrix.h"

namespace problem {

template <typename Field>
struct LP {
  CSCMatrix<Field> matrix;
  std::vector<Field> cost;

  Field cost_offset{0};

  std::vector<Bound<Field>> var_bounds;
  std::vector<Bound<Field>> rhs_bounds;

  std::vector<Bound<Field>> implied_var_bounds;

  std::string name;
  std::string cost_name;

  std::vector<std::string> var_names;
  std::vector<std::string> row_names;

  bool proven_infeasible{false};
  bool proven_unbounded{false};

  // for debugging purposes, throws if problem is not correct
  void validate() const {
    const auto [n, d] = matrix.shape();

    if (cost.size() != d) {
      throw std::runtime_error("Wrong cost vector size.");
    }

    if (var_bounds.size() != d) {
      throw std::runtime_error("Wrong variable bounds vector size.");
    }

    if (rhs_bounds.size() != n) {
      throw std::runtime_error("Wrong RHS bounds vector size.");
    }

    if (implied_var_bounds.size() != d) {
      throw std::runtime_error("Wrong implied variable bounds vector size.");
    }

    if (var_names.size() != d) {
      throw std::runtime_error("Wrong variable names vector size.");
    }

    if (row_names.size() != n) {
      throw std::runtime_error("Wrong row names vector size.");
    }
  }
};

namespace detail {

template <typename Field>
void print_sparse_row(std::ostream& os, const SparseVector<Field>& row,
                      const std::vector<std::string>& names) {
  for (size_t i = 0; i < row.size(); ++i) {
    if (i != 0 && row[i].second >= 0) {
      os << " + ";
    } else {
      os << " ";
    }

    auto var_name = names[row[i].first].empty()
                        ? "x" + std::to_string(row[i].first)
                        : names[row[i].first];

    os << row[i].second << " " << var_name;
  }
}

}  // namespace detail

template <typename Field>
std::ostream& operator<<(std::ostream& os, const LP<Field>& problem) {
  std::println(os, "Problem: {}", problem.name);
  std::println(os, "  Status: proven_infeasible = {}, proven_unbounded = {}",
               problem.proven_infeasible, problem.proven_unbounded);

  SparseVector<Field> sparse_cost;
  for (size_t i = 0; i < problem.cost.size(); ++i) {
    if (FieldTraits<Field>::is_nonzero(problem.cost[i])) {
      sparse_cost.emplace_back(i, problem.cost[i]);
    }
  }

  os << "  min " << problem.cost_name << " = ";
  detail::print_sparse_row(os, sparse_cost, problem.var_names);
  os << " + " << problem.cost_offset << "\n";

  const auto transposed = problem.matrix.get_transposed();
  for (size_t i = 0; i < transposed.size(); ++i) {
    auto row_name = problem.row_names[i].empty() ? "r" + std::to_string(i)
                                                 : problem.row_names[i];

    os << problem.row_names[i] << ":\t";

    detail::print_sparse_row(os, transposed[i], problem.var_names);

    os << " \\in " << problem.rhs_bounds[i] << "\n";
  }

  return os;
}

}  // namespace problem
