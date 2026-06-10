#pragma once

#include "LP.h"

namespace problem {

template <typename Field>
struct MILP : LP<Field> {
  std::vector<bool> is_integer;
  std::vector<bool> implied_is_integer;

  // for debugging purposes, throws if problem is not correct
  void validate() const {
    LP<Field>::validate();

    const auto [n, d] = this->matrix.shape();

    if (is_integer.size() != d) {
      throw std::runtime_error("Wrong integrality vector size.");
    }

    if (implied_is_integer.size() != d) {
      throw std::runtime_error("Wrong implied integrality vector size.");
    }
  }
};

template <typename Field>
std::ostream& operator<<(std::ostream& os, const MILP<Field>& problem) {
  using std::abs;

  std::println(os, "Problem: {}", problem.name);
  std::println(os, "  Status: proven_infeasible = {}, proven_unbounded = {}",
               problem.proven_infeasible, problem.proven_unbounded);

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

    os << " in " << problem.rhs_bounds[row] << "\n";
  }

  // print variables
  for (size_t i = 0; i < problem.matrix.cols(); ++i) {
    os << problem.var_name(i) << " in " << problem.var_bounds[i];

    if (problem.implied_var_bounds[i] != problem.var_bounds[i]) {
      os << "(implied bound: " << problem.implied_var_bounds[i] << ")";
    }

    if (problem.is_integer[i]) {
      os << " and integer";
    } else if (problem.implied_is_integer[i]) {
      os << "and implied integer";
    }

    os << "\n";
  }

  return os;
}

}  // namespace problem
