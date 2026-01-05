#pragma once

#include <vector>

#include "MILPProblem.h"
#include "linear/matrix/Matrix.h"

template <typename Field>
struct MILPProblemAsMatrices {
  Matrix<Field> A;
  Matrix<Field> b;
  Matrix<Field> c;

  std::vector<Field> lower;
  std::vector<Field> upper;

  std::vector<VariableType> variables;
};

template <typename Field>
MILPProblemAsMatrices<Field> to_matrices(const MILPProblem<Field>& problem) {
  size_t n = problem.constraints.size();
  size_t d = problem.variables.size();

  auto enumeration = problem.enumerate_variables();

  Matrix<Field> A(n, d, 0);
  Matrix<Field> b(n, 1, 0);
  Matrix<Field> c(1, d, 0);

  for (const auto& [var, coef] : problem.objective.get_variables()) {
    c[0, enumeration.at(var)] = coef;
  }

  for (size_t i = 0; i < n; ++i) {
    const Constraint<Field>& constraint = problem.constraints[i];

    if (constraint.type != ConstraintType::EQUAL_ZERO) {
      throw std::runtime_error(
          "Non-equality constraints can not be transformed into matrices.");
    }

    for (const auto& [var, coef] : constraint.expr.get_variables()) {
      A[i, enumeration.at(var)] = coef;
    }

    b[i, 0] = -constraint.expr.get_shift();
  }

  std::vector<VariableType> types(d);
  std::vector<Field> lower(d);
  std::vector<Field> upper(d);

  for (const auto& info : problem.variables) {
    types[enumeration.at(info.name)] = info.type;
    lower[enumeration.at(info.name)] = info.lower_bound;
    upper[enumeration.at(info.name)] = info.upper_bound;
  }

  return {A, b, c, lower, upper, types};
}
