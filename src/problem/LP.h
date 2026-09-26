#pragma once

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CoreLP.h"
#include "detail/ExpressionPrinter.h"
#include "linalg/Linalg.h"
#include "problem/Bound.h"

namespace problem {

// LP problem representation suitable for solving.
template <typename Field>
struct StandardLP;

// LP problem representation suitable for presolve.
template <typename Field>
struct LP : CoreLP<Field> {
  std::vector<Bound<Field>> rhs_bounds;

  std::vector<Bound<Field>> implied_var_bounds;

  bool proven_infeasible{false};
  bool proven_unbounded{false};

  LP() = default;

  // Note: read comment in similar CoreLP constructor
  LP(size_t rows, size_t cols)
      : CoreLP<Field>(rows, cols), rhs_bounds(rows), implied_var_bounds(cols) {}

  explicit LP(const StandardLP<Field>&);

  // for debugging purposes, throws if problem is not correct
  void validate() const {
    CoreLP<Field>::validate();

    const auto [n, d] = this->matrix.shape();

    if (rhs_bounds.size() != n) {
      throw std::runtime_error("Wrong RHS bounds vector size.");
    }

    if (implied_var_bounds.size() != d) {
      throw std::runtime_error("Wrong implied variable bounds vector size.");
    }
  }
};

template <typename Field>
std::ostream& operator<<(std::ostream& os, const LP<Field>& problem) {
  using std::abs;

  std::println(os, "Problem: {}", problem.name);
  std::println(os, "  Status: proven_infeasible = {}, proven_unbounded = {}",
               problem.proven_infeasible, problem.proven_unbounded);

  // print cost
  {
    detail::ExpressionPrinter printer{os};

    printer.name("max " + problem.safe_cost_name());

    for (size_t i = 0; i < problem.cost.size(); ++i) {
      if (abs(problem.cost[i]) > FieldTraits<Field>::tolerance) {
        printer.print(problem.cost[i], problem.var_name(i));
      }
    }

    if (abs(problem.cost_offset) > FieldTraits<Field>::tolerance) {
      printer.print(problem.cost_offset);
    }
    os << "\n";
  }

  // print constraints
  for (size_t row = 0; row < problem.matrix.rows(); ++row) {
    detail::ExpressionPrinter printer{os};

    printer.name(problem.row_name(row));

    for (size_t col = 0; col < problem.matrix.cols(); ++col) {
      if (auto value = problem.matrix.at(row, col)) {
        printer.print(value, problem.var_name(col));
      }
    }

    os << " in " << problem.rhs_bounds[row] << "\n";
  }

  // print variables
  for (size_t i = 0; i < problem.matrix.cols(); ++i) {
    os << problem.var_name(i) << " in " << problem.var_bounds[i];

    if (problem.implied_var_bounds[i] != problem.var_bounds[i]) {
      os << "(implied bound: " << problem.implied_var_bounds[i] << ")";
    }

    os << "\n";
  }

  return os;
}

}  // namespace problem
