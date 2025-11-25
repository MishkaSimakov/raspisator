#include <iostream>

#include "../src/linear/bb/BranchAndBound.h"
#include "../src/linear/matrix/RowBasis.h"
#include "linear/BigInteger.h"
#include "linear/MPS.h"
#include "linear/SimplexMethod.h"
#include "linear/builder/ProblemBuilder.h"

using Field = double;

int main() {
  // MILP formulation is from this paper:
  // https://arxiv.org/pdf/2504.20017v2

  size_t n = 3;
  size_t a_min = 1;

  size_t a_max = a_min + n * n - 1;
  size_t magic_constant = n * (a_min + a_max) / 2;

  auto builder = ProblemBuilder<Field>();
  std::vector<std::vector<std::vector<Variable<Field>>>> vars(n);

  for (size_t i = 0; i < n; ++i) {
    vars[i].resize(n);
    for (size_t j = 0; j < n; ++j) {
      for (size_t k = 0; k < n * n; ++k) {
        vars[i][j].push_back(builder.new_variable(
            fmt::format("x({}, {}, {})", i, j, k), VariableType::INTEGER));
      }
    }
  }

  // 2.1 a
  for (size_t k = 0; k < n * n; ++k) {
    Expression<Field> k_count;

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        k_count += vars[i][j][k];
      }
    }

    builder.add_constraint(k_count == Expression<Field>{1});
  }

  // 2.1 b
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      Expression<Field> cell_count;

      for (size_t k = 0; k < n * n; ++k) {
        cell_count += vars[i][j][k];
      }

      builder.add_constraint(cell_count == Expression<Field>{1});
    }
  }

  // 2.1 c
  for (size_t i = 0; i < n; ++i) {
    Expression<Field> row_sum;

    for (size_t j = 0; j < n; ++j) {
      for (size_t k = 0; k < n * n; ++k) {
        row_sum += vars[i][j][k] * (k + a_min);
      }
    }

    builder.add_constraint(row_sum == Expression<Field>(magic_constant));
  }

  // 2.1 d
  for (size_t j = 0; j < n; ++j) {
    Expression<Field> col_sum;

    for (size_t i = 0; i < n; ++i) {
      for (size_t k = 0; k < n * n; ++k) {
        col_sum += vars[i][j][k] * (k + a_min);
      }
    }

    builder.add_constraint(col_sum == Expression<Field>(magic_constant));
  }

  // 2.1 e
  {
    Expression<Field> diagonal_sum;

    for (size_t i = 0; i < n; ++i) {
      for (size_t k = 0; k < n * n; ++k) {
        diagonal_sum += vars[i][i][k] * (k + a_min);
      }
    }

    builder.add_constraint(diagonal_sum == Expression<Field>(magic_constant));
  }

  // 2.1 d
  {
    Expression<Field> diagonal_sum;

    for (size_t i = 0; i < n; ++i) {
      for (size_t k = 0; k < n * n; ++k) {
        diagonal_sum += vars[n - i - 1][i][k] * (k + a_min);
      }
    }

    builder.add_constraint(diagonal_sum == Expression<Field>(magic_constant));
  }

  // solve the problem
  auto problem = builder.get_problem();

  BranchAndBound<Field, SimplexMethod<Field>> solver(problem);
  auto solution = solver.solve();

  std::cout << GraphvizBuilder<Field>().build(solver.get_root()) << std::endl;

  if (std::holds_alternative<FiniteMILPSolution<Field>>(solution)) {
    auto point = std::get<FiniteMILPSolution<Field>>(solution).point;
    std::cout << "point: " << linalg::transposed(point) << std::endl;

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        for (size_t k = 0; k < n * n; ++k) {
          if (FieldTraits<Field>::is_strictly_positive(
                  builder.extract_variable(point, vars[i][j][k]))) {
            std::cout << k + a_min << " ";
            break;
          }
        }
      }

      std::cout << "\n";
    }
  }

  return 0;
}
