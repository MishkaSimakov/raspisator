#pragma once
#include <variant>

#include "Matrix.h"
#include "linear/Solution.h"
#include "utils/Variant.h"

// basic feasible solution
template <typename Field>
struct BFS {
  Matrix<Field> point;
  std::vector<size_t> basic_variables;

  static BFS construct_nondegenerate(Matrix<Field> point) {
    if (point.get_width() != 1) {
      throw std::invalid_argument("shape of the point must be (d, 1).");
    }

    std::vector<size_t> basic_variables;
    for (size_t i = 0; i < point.get_height(); ++i) {
      if (point[i, 0] != 0) {
        basic_variables.push_back(i);
      }
    }

    return BFS(point, basic_variables);
  }
};

// Solves cx -> max, Ax = b, x >= 0
// A is (n, d) matrix, b is (n, 1) matrix, c is (1, d) matrix
// it is assumed that n < d

// Important Note: this method is NOT numerically stable!
// Therefore, the result may be incorrect when using Field = double or float.
// Use Rational instead.
template <typename Field>
class SimplexMethod {
  Matrix<Field> A_;
  Matrix<Field> b_;
  Matrix<Field> c_;

 public:
  SimplexMethod(Matrix<Field> A, Matrix<Field> b, Matrix<Field> c)
      : A_(std::move(A)), b_(std::move(b)), c_(std::move(c)) {
    // check sizes
    auto [n, d] = A_.shape();

    if (n >= d) {
      throw std::invalid_argument(
          "Solve a system of linear equations instead.");
    }

    if (b_.shape() != std::pair{n, 1}) {
      throw std::invalid_argument("Matrix b has wrong dimensions.");
    }

    if (c_.shape() != std::pair{1, d}) {
      throw std::invalid_argument("Matrix c has wrong dimensions.");
    }
  }

  // bfs is (d, 1) matrix
  // returns tableau and a vector holding indices of basic variables
  Matrix<Field> initialize_tableau(BFS<Field> bfs) {
    auto [n, d] = A_.shape();

    auto tableau = Matrix<Field>(n + 1, d + 1);

    // build A_b
    auto A_b = Matrix<Field>(n, n);
    auto c_b = Matrix<Field>(1, n);

    for (size_t i = 0; i < bfs.basic_variables.size(); ++i) {
      // copy i-th column into A_b
      c_b[0, i] = c_[0, bfs.basic_variables[i]];

      for (size_t j = 0; j < n; ++j) {
        A_b[j, i] = A_[j, bfs.basic_variables[i]];
      }
    }

    //
    auto inv = A_b.inverse();

    auto X = inv * A_;
    auto deltas = c_b * X - c_;
    auto b_x = inv * b_;
    auto value = c_b * b_x;

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < d; ++j) {
        tableau[i, j + 1] = X[i, j];
      }
    }

    for (size_t i = 0; i < n; ++i) {
      tableau[i, 0] = b_x[i, 0];
    }

    for (size_t j = 0; j < d; ++j) {
      tableau[n, j + 1] = deltas[0, j];
    }

    tableau[n, 0] = value[0, 0];

    return tableau;
  }

  static std::pair<size_t, Field> find_entering_variable(
      const Matrix<Field>& tableau) {
    Field min_delta_value = 0;
    size_t min_delta_index = 0;

    // TODO:
    // it is said in this article that in case of tie
    // decision variables must have priority over slack variables
    // https://www.nascollege.org/econtent/ecotent-10-4-20/DR%20K%20K%20KANSAL/L%2010%20M%20COM%2020-4-E.pdf
    for (size_t j = 0; j < tableau.get_width() - 1; ++j) {
      Field delta = tableau[tableau.get_height() - 1, j + 1];

      if (delta <= min_delta_value) {
        min_delta_value = delta;
        min_delta_index = j;
      }
    }

    return std::pair{min_delta_index, min_delta_value};
  }

  static void pivot(Matrix<Field>& tableau, std::vector<size_t>& basic_vars,
                    size_t entering_var, size_t leaving_var) {
    basic_vars[leaving_var] = entering_var;
    tableau.gaussian_elimination(leaving_var, entering_var + 1);
  }

  static std::optional<size_t> find_leaving_variable(
      const Matrix<Field>& tableau, size_t entering_var) {
    std::optional<Field> min_t_value = std::nullopt;
    size_t min_t_index = 0;

    for (size_t i = 0; i < tableau.get_height() - 1; ++i) {
      if (tableau[i, entering_var + 1] <= 0) {
        continue;
      }

      Field t = tableau[i, 0] / tableau[i, entering_var + 1];

      if (!min_t_value || t < *min_t_value) {
        min_t_value = t;
        min_t_index = i;
      }
    }

    return min_t_value ? std::optional{min_t_index} : std::nullopt;
  }

  // finds maximum starting from bfs (basic feasible solution)
  std::variant<FiniteLPSolution<Field>, InfiniteSolution> solve_from(
      BFS<Field> bfs) {
    auto [n, d] = A_.shape();

    if (bfs.point.shape() != std::pair{d, 1}) {
      throw std::invalid_argument("bfs has wrong shape.");
    }

    auto tableau = initialize_tableau(bfs);

    while (true) {
      auto [entering_var, entering_var_delta] = find_entering_variable(tableau);

      if (entering_var_delta >= 0) {
        // solution is found
        Matrix<Field> point(d, 1, 0);

        for (size_t i = 0; i < n; ++i) {
          point[bfs.basic_variables[i], 0] = tableau[i, 0];
        }

        Field value = tableau[n, 0];

        return FiniteLPSolution{std::move(point),
                                std::move(bfs.basic_variables),
                                std::move(value), std::move(tableau)};
      }

      auto leaving_var = find_leaving_variable(tableau, entering_var);

      if (!leaving_var) {
        return InfiniteSolution{};
      }

      pivot(tableau, bfs.basic_variables, entering_var, *leaving_var);
    }
  }

  // algorithm is taken from
  // https://people.orie.cornell.edu/dpw/orie6300/Lectures/lec12.pdf
  std::optional<BFS<Field>> find_bfs() {
    auto [n, d] = A_.shape();

    // ensure that b >= 0
    for (size_t i = 0; i < n; ++i) {
      if (b_[i, 0] < 0) {
        // multiply whole equation by -1
        b_[i, 0] = -b_[i, 0];

        for (size_t j = 0; j < d; ++j) {
          A_[i, j] = -A_[i, j];
        }
      }
    }

    // add artificial variables
    auto A_new = A_.get_extended(n, d + n, 0);
    auto c_new = Matrix<Field>(1, d + n, 0);
    auto bfs_new = Matrix<Field>(d + n, 1, 0);

    for (size_t i = 0; i < n; ++i) {
      A_new[i, i + d] = 1;
      c_new[0, i + d] = -1;
    }

    for (size_t i = 0; i < n; ++i) {
      bfs_new[i + d, 0] = b_[i, 0];
    }

    std::vector<size_t> basic_vars(n);
    for (size_t i = 0; i < n; ++i) {
      basic_vars[i] = d + i;
    }

    auto solver = SimplexMethod(std::move(A_new), b_, std::move(c_new));

    // InfiniteSolution is impossible here
    auto solution = std::get<FiniteLPSolution<Field>>(
        solver.solve_from(BFS{bfs_new, basic_vars}));

    // Case 1
    if (solution.value < 0) {
      return std::nullopt;
    }

    // Case 2: try to eliminate artificial variables from basic variables (if
    // there are any) using pivot operation

    for (size_t i = 0; i < solution.basic_variables.size(); ++i) {
      if (solution.basic_variables[i] < d) {
        continue;
      }

      // artificial basic variable -> try to pivot
      bool pivoted = false;

      for (size_t j = 0; j < d; ++j) {
        if (solution.tableau[i, j + 1] != 0) {
          SimplexMethod::pivot(solution.tableau, solution.basic_variables, j,
                               i);
          pivoted = true;
        }
      }

      if (!pivoted) {
        // rows of A are linearly dependent
        throw std::runtime_error(
            "Linear dependent rows in A. This case is not implemented yet.");
      }
    }

    auto bfs = Matrix<Field>(d, 1);
    for (size_t i = 0; i < n; ++i) {
      bfs[i, 0] = solution.point[i, 0];
    }

    return BFS<Field>{bfs, std::move(solution.basic_variables)};
  }

  // same as solve_from, but automatically finds bfs
  using SimplexSolution = std::variant<FiniteLPSolution<Field>,
                                       InfiniteSolution, NoFeasibleElements>;

  SimplexSolution solve() {
    auto bfs = find_bfs();

    if (!bfs.has_value()) {
      return NoFeasibleElements{};
    }

    return variant_cast<SimplexSolution>(solve_from(*bfs));
  }
};
