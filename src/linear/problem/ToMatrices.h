#pragma once

#include <vector>

#include "MILPProblem.h"
#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "linear/model/LP.h"

using linalg::Matrix, linalg::Vector;

template <typename Field>
struct MILPProblemAsMatrices {
  Matrix<Field> A;
  Vector<Field> b;
  Vector<Field> c;

  Bounds<Field> bounds;
  std::vector<VariableType> variables;
};

template <typename Field>
MILPProblemAsMatrices<Field> to_matrices(const MILPProblem<Field>& problem) {
  size_t n = problem.constraints.size();
  size_t d = problem.variables.size();

  auto enumeration = problem.enumerate_variables();

  Matrix<Field> A(n, d);
  Vector<Field> b(n);
  Vector<Field> c(d);

  for (const auto& [var, coef] : problem.objective.get_variables()) {
    c[enumeration.at(var)] = coef;
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

    b[i] = -constraint.expr.get_shift();
  }

  std::vector<VariableType> types(d);
  Bounds<Field> bounds(d);

  for (const auto& info : problem.variables) {
    types[enumeration.at(info.name)] = info.type;
    bounds[enumeration.at(info.name)] = info.bound;
  }

  return {A, b, c, bounds, types};
}
