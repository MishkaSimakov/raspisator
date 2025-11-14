#pragma once
#include <variant>

#include "Matrix.h"

template <typename Field>
struct FiniteSolution {
  Matrix<Field> point;
  Field value;
};

struct InfiniteSolution {};

// Solves cx -> max, Ax = b, x >= 0
// A is (n, d) matrix, b is (n, 1) matrix, c is (1, d) matrix
// it is assumed that n < d
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
  std::pair<Matrix<Field>, std::vector<size_t>> initialize_tableau(
      const Matrix<Field>& bfs) {
    auto [n, d] = A_.shape();

    auto tableau = Matrix<Field>(n + 1, d + 1);

    // build A_b
    auto A_b = Matrix<Field>(n, n);
    auto c_b = Matrix<Field>(1, n);

    size_t basic_vars_counter = 0;
    std::vector<size_t> basic_vars_indices(n);

    for (size_t i = 0; i < d; ++i) {
      if (bfs[i, 0] == 0) {  // TODO: Field?
        continue;
      }

      // copy i-th column into A_b
      basic_vars_indices[basic_vars_counter] = i;
      c_b[0, basic_vars_counter] = c_[0, i];

      for (size_t j = 0; j < n; ++j) {
        A_b[j, basic_vars_counter] = A_[j, i];
      }

      ++basic_vars_counter;
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

    return {tableau, basic_vars_indices};
  }

  std::pair<size_t, Field> find_entering_variable(
      const Matrix<Field>& tableau) const {
    Field min_delta_value = 0;
    size_t min_delta_index = 0;

    auto [d, n] = A_.shape();

    for (size_t j = 0; j < n; ++j) {
      Field delta = tableau[d, j + 1];

      if (delta <= min_delta_value) {
        min_delta_value = delta;
        min_delta_index = j;
      }
    }

    return std::pair{min_delta_index, min_delta_value};
  }

  void transform_tableau(Matrix<Field>& tableau,
                         std::vector<size_t>& basic_vars, size_t entering_var,
                         size_t leaving_var) {
    basic_vars[leaving_var] = entering_var;
    tableau.gaussian_elimination(leaving_var, entering_var + 1);
  }

  std::optional<size_t> find_leaving_variable(const Matrix<Field>& tableau,
                                              size_t entering_var) const {
    auto [d, n] = A_.shape();

    std::optional<Field> min_t_value = std::nullopt;
    size_t min_t_index = 0;

    for (size_t i = 0; i < d; ++i) {
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
  // returns solution ( (d, 1) matrix) and maximum value
  std::variant<FiniteSolution<Field>, InfiniteSolution> solve_from(
      const Matrix<Field>& bfs) {
    auto [n, d] = A_.shape();

    if (bfs.shape() != std::pair{d, 1}) {
      throw std::invalid_argument("bfs has wrong shape.");
    }

    auto [tableau, basic_vars] = initialize_tableau(bfs);

    while (true) {
      auto [entering_var, entering_var_delta] = find_entering_variable(tableau);

      if (entering_var_delta >= 0) {
        // solution is found
        Matrix<Field> point(d, 1);
        for (size_t i = 0; i < n; ++i) {
          point[basic_vars[i], 0] = tableau[i, 0];
        }

        Field value = tableau[n, 0];
        return FiniteSolution<Field>{point, value};
      }

      auto leaving_var = find_leaving_variable(tableau, entering_var);

      if (!leaving_var) {
        return InfiniteSolution{};
      }

      transform_tableau(tableau, basic_vars, entering_var, *leaving_var);
    }
  }
};
