#pragma once

#include "MILP.h"
#include "StandardLP.h"

namespace problem {

template <typename Field>
struct StandardMILP : StandardLP<Field> {
  std::vector<bool> is_integer;

  StandardMILP() = default;

  // Note: read comment in similar CoreLP constructor
  StandardMILP(size_t rows, size_t cols)
      : StandardLP<Field>(rows, cols), is_integer(cols) {}

  explicit StandardMILP(const MILP<Field>& other)
      : StandardLP<Field>(other), is_integer(other.is_integer) {}

  // for debugging purposes, throws if problem is not correct
  void validate() const {
    StandardLP<Field>::validate();

    const auto [n, d] = this->matrix.shape();

    if (is_integer.size() != d) {
      throw std::runtime_error("Wrong integrality vector size.");
    }
  }
};

template <typename Field>
MILP<Field>::MILP(const StandardMILP<Field>& other) : LP<Field>(other) {
  is_integer = other.is_integer;
  implied_is_integer = other.is_integer;
}

template <typename Field>
std::ostream& operator<<(std::ostream& os, const StandardMILP<Field>& problem) {
  using std::abs;

  std::println(os, "Problem: {}", problem.name);

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

    os << " = " << problem.rhs[row] << "\n";
  }

  // print variables
  for (size_t i = 0; i < problem.matrix.cols(); ++i) {
    os << problem.var_name(i) << " in " << problem.var_bounds[i];

    if (problem.implied_var_bounds[i] != problem.var_bounds[i]) {
      os << "(implied bound: " << problem.implied_var_bounds[i] << ")";
    }

    if (problem.is_integer[i]) {
      os << " and integer";
    }

    os << "\n";
  }

  return os;
}

}  // namespace problem
