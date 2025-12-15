#pragma once

#include <algorithm>

#include "BaseOptimizer.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/RowBasis.h"
#include "linear/problem/MILPProblem.h"

template <typename Field>
class RemoveLinearlyDependentConstraints final : public BaseOptimizer<Field> {
  size_t equalities_count(const MILPProblem<Field>& problem) {
    size_t count = 0;

    for (const auto& constraint : problem.constraints) {
      if (constraint.type == ConstraintType::EQUAL) {
        ++count;
      }
    }

    return count;
  }

  std::pair<Matrix<Field>, Matrix<Field>> as_matrix(
      const MILPProblem<Field>& problem) {
    size_t n = equalities_count(problem);
    size_t d = problem.variables.size();

    auto enumeration = problem.enumerate_variables();

    Matrix<Field> A(n, d, 0);
    Matrix<Field> b(n, 1, 0);

    size_t i = 0;
    for (const auto& constraint : problem.constraints) {
      if (constraint.type != ConstraintType::EQUAL) {
        continue;
      }

      for (const auto& [var, coef] : constraint.lhs.get_variables()) {
        A[i, enumeration.at(var)] = coef;
      }

      b[i, 0] = constraint.rhs;
      ++i;
    }

    return std::pair{A, b};
  }

 public:
  RemoveLinearlyDependentConstraints() = default;

  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    auto [A, b] = as_matrix(problem);

    auto row_basis = linalg::get_row_basis(A);

    // check feasibility
    auto row_basis_with_b = linalg::get_row_basis(linalg::hstack(A, b));
    if (row_basis.size() < row_basis_with_b.size()) {
      throw std::runtime_error("Provided system is trivially infeasible.");
    }

    // preserve only linearly independent rows
    size_t equality_index = 0;
    for (auto itr = problem.constraints.begin();
         itr != problem.constraints.end();) {
      if (itr->type != ConstraintType::EQUAL) {
        ++itr;
        continue;
      }

      if (std::ranges::find(row_basis, equality_index) != row_basis.end()) {
        ++itr;
      } else {
        itr = problem.constraints.erase(itr);
      }

      ++equality_index;
    }

    return problem;
  }

  Matrix<Field> inverse(const Matrix<Field>& point) override { return point; }
};
