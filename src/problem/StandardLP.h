#pragma once

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CoreLP.h"
#include "LP.h"
#include "detail/ExpressionPrinter.h"
#include "linalg/Linalg.h"

namespace problem {

// LP problem representation suitable for solving.
template <typename Field>
struct StandardLP : CoreLP<Field> {
  Vector<Field> rhs;

  StandardLP() = default;

  // Note: read comment in similar CoreLP constructor
  StandardLP(size_t rows, size_t cols) : CoreLP<Field>(rows, cols), rhs(rows) {}

  explicit StandardLP(const LP<Field>& other)
      : CoreLP<Field>(other), rhs(other.rhs_bounds.size()) {
    for (size_t i = 0; i < other.rhs_bounds.size(); ++i) {
      if (!other.rhs_bounds[i].is_fixed()) {
        throw std::invalid_argument(
            "Can not convert LP to StandardLP because some constraints have "
            "non-zero range.");
      }

      rhs[i] = *other.rhs_bounds[i].lower;
    }
  }

  // for debugging purposes, throws if problem is not correct
  void validate() const {
    CoreLP<Field>::validate();

    const auto [n, d] = this->matrix.shape();

    if (rhs.size() != n) {
      throw std::runtime_error("Wrong RHS vector size.");
    }
  }
};

template <typename Field>
LP<Field>::LP(const StandardLP<Field>& other) : CoreLP<Field>(other) {
  implied_var_bounds = other.var_bounds;

  rhs_bounds.resize(other.rhs.size());
  for (size_t i = 0; i < other.rhs.size(); ++i) {
    rhs_bounds[i] = Bound<Field>{other.rhs[i], other.rhs[i]};
  }
}

template <typename Field>
std::ostream& operator<<(std::ostream& os, const StandardLP<Field>& problem) {
  using std::abs;

  std::println(os, "Problem: {}", problem.name);

  // print cost
  {
    detail::ExpressionPrinter printer{os};

    printer.name(problem.cost_name);

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

    os << " = " << problem.rhs[row] << "\n";
  }

  // print variables
  for (size_t i = 0; i < problem.matrix.cols(); ++i) {
    os << problem.var_name(i) << " in " << problem.var_bounds[i] << "\n";
  }

  return os;
}

}  // namespace problem
