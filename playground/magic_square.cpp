#include <iostream>

#include "../src/linear/matrix/RowBasis.h"
#include "linear/BigInteger.h"
#include "linear/BranchAndBound.h"
#include "linear/MPS.h"
#include "linear/SimplexMethod.h"

using Field = Rational;

int main() {
  // MILP formulation is from this paper:
  // https://arxiv.org/pdf/2504.20017v2

  size_t n = 3;
  size_t a_min = 1;

  size_t a_max = a_min + n * n - 1;
  size_t variables_count = n * n * n * n;
  size_t magic_constant = n * (a_min + a_max) / 2;

  auto get_var_index = [n](size_t i, size_t j, size_t k) {
    return k * n * n + j * n + i;
  };

  Matrix<Field> A(0, variables_count);
  Matrix<Field> b(0, 1);
  Matrix<Field> c(1, variables_count, 0);

  // 2.1 a
  for (size_t k = 0; k < n * n; ++k) {
    Matrix<Field> constraint(1, variables_count);

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        constraint[0, get_var_index(i, j, k)] = 1;
      }
    }

    A = linalg::vstack(A, constraint);
    b = linalg::vstack(b, Matrix<Field>::item(1));
  }

  // 2.1 b
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      Matrix<Field> constraint(1, variables_count);

      for (size_t k = 0; k < n * n; ++k) {
        constraint[0, get_var_index(i, j, k)] = 1;
      }

      A = linalg::vstack(A, constraint);
      b = linalg::vstack(b, Matrix<Field>::item(1));
    }
  }

  // 2.1 c
  for (size_t i = 0; i < n; ++i) {
    Matrix<Field> constraint(1, variables_count);

    for (size_t j = 0; j < n; ++j) {
      for (size_t k = 0; k < n * n; ++k) {
        constraint[0, get_var_index(i, j, k)] = k + a_min;
      }
    }

    A = linalg::vstack(A, constraint);
    b = linalg::vstack(b, Matrix<Field>::item(magic_constant));
  }

  // 2.1 d
  for (size_t j = 0; j < n; ++j) {
    Matrix<Field> constraint(1, variables_count);

    for (size_t i = 0; i < n; ++i) {
      for (size_t k = 0; k < n * n; ++k) {
        constraint[0, get_var_index(i, j, k)] = k + a_min;
      }
    }

    A = linalg::vstack(A, constraint);
    b = linalg::vstack(b, Matrix<Field>::item(magic_constant));
  }

  // 2.1 e
  {
    Matrix<Field> constraint(1, variables_count);

    for (size_t i = 0; i < n; ++i) {
      for (size_t k = 0; k < n * n; ++k) {
        constraint[0, get_var_index(i, i, k)] = k + a_min;
      }
    }

    A = linalg::vstack(A, constraint);
    b = linalg::vstack(b, Matrix<Field>::item(magic_constant));
  }

  // 2.1 d
  {
    Matrix<Field> constraint(1, variables_count);

    for (size_t i = 0; i < n; ++i) {
      for (size_t k = 0; k < n * n; ++k) {
        constraint[0, get_var_index(n - i - 1, i, k)] = k + a_min;
      }
    }

    A = linalg::vstack(A, constraint);
    b = linalg::vstack(b, Matrix<Field>::item(magic_constant));
  }

  // all variables are integer
  std::vector<VariableType> integer(variables_count);
  for (size_t i = 0; i < variables_count; ++i) {
    integer[i] = VariableType::INTEGER;
  }

  // resulting constraints are not linearly independent
  // perform row elimination
  auto row_basis = linalg::get_row_basis(A);

  Matrix<Field> reduced_A(row_basis.size(), variables_count);
  Matrix<Field> reduced_b(row_basis.size(), 1);

  for (size_t i = 0; i < row_basis.size(); ++i) {
    for (size_t col = 0; col < variables_count; ++col) {
      reduced_A[i, col] = A[row_basis[i], col];
    }

    reduced_b[i, 0] = b[row_basis[i], 0];
  }

  std::cout << reduced_A << std::endl;
  std::cout << reduced_b << std::endl;
  std::cout << c << std::endl;

  // solve the problem
  MILPProblem<Field> problem(reduced_A, reduced_b, c, integer);

  BranchAndBound<Field, SimplexMethod<Field>> solver(problem);
  auto solution = solver.solve();

  std::cout << GraphvizBuilder<Field>().build(solver.get_root()) << std::endl;

  if (std::holds_alternative<FiniteMILPSolution<Field>>(solution)) {
    std::cout << "solution: "
              << std::get<FiniteMILPSolution<Field>>(solution).point
              << std::endl;
  }

  return 0;
}
