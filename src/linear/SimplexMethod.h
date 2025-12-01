#pragma once

#include <algorithm>
#include <fstream>
#include <ranges>
#include <variant>

#include "linear/model/LP.h"
#include "matrix/LU.h"
#include "matrix/Matrix.h"
#include "matrix/RowBasis.h"
#include "sparse/LU.h"
#include "utils/Variant.h"

// Double, double toil and trouble;
// Fire burn and caldron bubble.
// - Macbeth
//
// There are a lot of problems with a double. This version of SimplexMethod
// is trying to address them. It is trying to be numerically stable at least.

// Solves cx -> max, Ax = b, x >= 0
// A is (n, d) matrix, b is (n, 1) matrix, c is (1, d) matrix
// it is assumed that n < d
template <typename Field>
class SimplexMethod {
  CSCMatrix<Field> A_;

  Matrix<Field> b_;
  Matrix<Field> c_;

 public:
  SimplexMethod(CSCMatrix<Field> A, Matrix<Field> b, Matrix<Field> c)
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

  Field estimate_rounding_errors(const Matrix<Field>& tableau,
                                 const std::vector<size_t>& basic_vars) const {
    auto [n, d] = A_.shape();

    auto x = Matrix<Field>(d, 1);
    for (size_t i = 0; i < n; ++i) {
      x[basic_vars[i], 0] = tableau[i, 0];
    }

    auto error = linalg::to_dense(A_) * x - b_;

    Field max = 0;
    for (size_t i = 0; i < n; ++i) {
      max = std::max(max, FieldTraits<Field>::abs(error[i, 0]));
    }
    return max;
  }

  // TODO:
  // it is said in this article that in case of tie
  // decision variables must have priority over slack variables
  // https://www.nascollege.org/econtent/ecotent-10-4-20/DR%20K%20K%20KANSAL/L%2010%20M%20COM%2020-4-E.pdf
  std::pair<size_t, Field> find_entering_variable(
      const Matrix<Field>& dual_point,
      const std::vector<size_t>& basic_vars) const {
    Field min_change = 0;
    size_t min_change_index = 0;

    auto [n, d] = A_.shape();

    for (size_t i = 0; i < d; ++i) {
      // TODO: store basic vars as a set?
      if (std::ranges::find(basic_vars, i) != basic_vars.end()) {
        continue;
      }

      Field change = -c_[0, i];
      for (const auto& [row, value] : A_.get_column(i)) {
        change += dual_point[row, 0] * value;
      }

      if (change < min_change) {
        min_change_index = i;
        min_change = change;
      }
    }

    return {min_change_index, min_change};
  }

  static void pivot(std::vector<size_t>& basic_vars, size_t entering_var,
                    size_t leaving_var) {
    basic_vars[leaving_var] = entering_var;
  }

  std::optional<size_t> find_leaving_variable(
      const Matrix<Field>& point, const Matrix<Field>& entering_coords) const {
    std::optional<Field> min_t_value = std::nullopt;
    size_t min_t_index = 0;

    size_t n = A_.shape().first;

    for (size_t i = 0; i < n; ++i) {
      if (!FieldTraits<Field>::is_strictly_positive(entering_coords[i, 0])) {
        continue;
      }

      Field t = point[i, 0] / entering_coords[i, 0];

      if (!min_t_value || t < *min_t_value) {
        min_t_value = t;
        min_t_index = i;
      }
    }

    return min_t_value ? std::optional{min_t_index} : std::nullopt;
  }

  auto get_point(const Matrix<Field>& L, const Matrix<Field>& U,
                 const std::vector<size_t>& P) {
    auto Pb = linalg::apply_permutation(b_, P);

    // solve Ly = Pb
    auto y = linalg::solve_lower(L, Pb, std::true_type{});

    // solve Ux = y
    auto x = linalg::solve_upper(U, y, std::false_type{});

    return x;
  }

  auto get_entering_coordinates(const CSCMatrix<Field>& L,
                                const CSCMatrix<Field>& U,
                                const std::vector<size_t>& P,
                                size_t entering_var) {
    size_t n = A_.shape().first;
    Matrix<Field> b(n, 1);

    for (const auto& [row, value] : A_.get_column(entering_var)) {
      b[row, 0] = value;
    }

    return linalg::solve_linear(L, U, P, b);
  }

  auto get_dual_point(const CSCMatrix<Field>& L, const CSCMatrix<Field>& U,
                      const std::vector<size_t>& P,
                      const std::vector<size_t>& basic_vars) {
    size_t n = A_.shape().first;

    Matrix<Field> cb(n, 1);
    for (size_t i = 0; i < n; ++i) {
      cb[i, 0] = c_[0, basic_vars[i]];
    }

    return linalg::solve_transposed_linear(L, U, P, cb);
  }

  // finds maximum starting from bfs (basic feasible solution)
  std::variant<FiniteLPSolution<Field>, InfiniteSolution> solve_from(
      BFS<Field> bfs) {
    auto [n, d] = A_.shape();

    if (bfs.point.shape() != std::pair{d, 1}) {
      throw std::invalid_argument("bfs has wrong shape.");
    }
    if (bfs.basic_variables.size() != n) {
      throw std::invalid_argument("bfs has wrong size.");
    }

    while (true) {
      auto [L, U, P] = linalg::sparse_lup(A_, bfs.basic_variables);

      // obtain point associated with given basic variables by solving
      // Bu = b
      auto point = linalg::solve_linear(L, U, P, b_);

      // obtain a solution of a dual problem by solving
      // B^T u = c_b
      auto dual_point = get_dual_point(L, U, P, bfs.basic_variables);

      auto [entering_var, change] =
          find_entering_variable(dual_point, bfs.basic_variables);

      if (!FieldTraits<Field>::is_strictly_negative(change)) {
        // solution is found
        Matrix<Field> extended_point(d, 1, 0);

        for (size_t i = 0; i < n; ++i) {
          extended_point[bfs.basic_variables[i], 0] = point[i, 0];
        }

        Field value = (c_ * extended_point)[0, 0];

        return FiniteLPSolution{std::move(extended_point),
                                std::move(bfs.basic_variables),
                                std::move(value)};
      }

      // obtain coordinates of entering variable by solving
      // By = A_s, where s is the entering variable index
      auto entering_coords = get_entering_coordinates(L, U, P, entering_var);

      auto leaving_var = find_leaving_variable(point, entering_coords);

      if (!leaving_var) {
        return InfiniteSolution{};
      }

      pivot(bfs.basic_variables, entering_var, *leaving_var);
    }
  }

  // This function tries to construct bfs using old bfs as a starting point
  // restrictions on old_bfs:
  // 1. old_bfs forms the largest linearly independent set of columns
  // 2. old_bfs satisfies Ax = b
  // 3. variable with index negative_index is the only
  //    negative variable in old_bfs
  std::optional<BFS<Field>> reconstruct_bfs(BFS<Field> old_bfs,
                                            size_t negative_index) {
    // replacing one column from bfs with artificial column
    // TODO: choice of replaced column may be important, investigate this

    auto [n, d] = A_.shape();

    auto c_new = Matrix<Field>(1, d + 1, 0);
    c_new[0, d] = -1;

    auto bfs_new = old_bfs.point.get_extended(d + 1, 1, 0);
    bfs_new[d, 0] = old_bfs.point[negative_index, 0] * -1;
    bfs_new[negative_index, 0] = 0;

    size_t index_in_bfs =
        std::ranges::find(old_bfs.basic_variables, negative_index) -
        old_bfs.basic_variables.begin();

    if (index_in_bfs == old_bfs.basic_variables.size()) {
      throw std::invalid_argument(
          "Negative index is not among basic variables. This means that "
          "provided bfs is incorrect.");
    }

    old_bfs.basic_variables[index_in_bfs] = d;

    CSCMatrix<Field> A_new(A_);
    A_new.add_column();
    for (const auto& [row, value] : A_.get_column(negative_index)) {
      A_new.push_to_last_column(row, -value);
    }

    auto solver = SimplexMethod(A_new, b_, c_new);

    // InfiniteSolution is impossible here
    auto solution = std::get<FiniteLPSolution<Field>>(
        solver.solve_from(BFS{bfs_new, old_bfs.basic_variables}));

    // Case 1
    if (FieldTraits<Field>::is_strictly_negative(solution.value)) {
      return std::nullopt;
    }

    // Case 2: try to eliminate artificial variables from basic variables (if
    // there are any) using pivot operation

    std::ranges::sort(solution.basic_variables);

    // Case 2 (a): all basic variables are real
    if (solution.basic_variables.back() < d) {
      return BFS<Field>{solution.point[{0, d}, 0],
                        std::move(solution.basic_variables)};
    }

    // Case 2 (b): there are artificial variables amongst basic variables.
    // Trying to replace them with a real ones.
    std::vector<size_t> real_basic_vars;
    for (size_t basic_var : solution.basic_variables) {
      if (basic_var < d) {
        real_basic_vars.push_back(basic_var);
      }
    }

    auto row_basis = linalg::complete_row_basis(
        linalg::to_dense_transposed(A_), real_basic_vars);

    if (row_basis.size() < n) {
      // rows of A are linearly dependent
      throw std::runtime_error(
          "Linear dependent rows in A. This case is not implemented yet.");
    }

    return BFS<Field>{solution.point[{0, d}, 0], std::move(row_basis)};
  }

  // algorithm is taken from
  // https://people.orie.cornell.edu/dpw/orie6300/Lectures/lec12.pdf
  std::optional<BFS<Field>> find_bfs() {
    auto [n, d] = A_.shape();

    // ensure that b >= 0
    for (size_t col = 0; col < d; ++col) {
      for (std::tuple<size_t&, Field&> element : A_.get_column(col)) {
        auto& [row, value] = element;
        if (b_[row, 0] < 0) {
          value *= -1;
        }
      }
    }
    for (size_t i = 0; i < n; ++i) {
      if (b_[i, 0] < 0) {
        b_[i, 0] *= -1;
      }
    }

    // add artificial variables
    CSCMatrix<Field> A_new(A_);
    auto c_new = Matrix<Field>(1, d + n, 0);
    auto bfs_new = Matrix<Field>(d + n, 1, 0);

    for (size_t i = 0; i < n; ++i) {
      c_new[0, i + d] = -1;
    }
    for (size_t i = 0; i < n; ++i) {
      A_new.add_column();
      A_new.push_to_last_column(i, 1);
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
    if (FieldTraits<Field>::is_strictly_negative(solution.value)) {
      return std::nullopt;
    }

    // Case 2: try to eliminate artificial variables from basic variables (if
    // there are any) using pivot operation

    std::ranges::sort(solution.basic_variables);

    // Case 2 (a): all basic variables are real
    if (solution.basic_variables.back() < d) {
      return BFS<Field>{solution.point[{0, d}, 0], solution.basic_variables};
    }

    // Case 2 (b): there are artificial variables amongst basic variables.
    // Trying to replace them with a real ones.
    std::vector<size_t> real_basic_vars;
    for (size_t basic_var : solution.basic_variables) {
      if (basic_var < d) {
        real_basic_vars.push_back(basic_var);
      }
    }

    auto row_basis = linalg::complete_row_basis(
        linalg::to_dense_transposed(A_), real_basic_vars);

    if (row_basis.size() < n) {
      // rows of A are linearly dependent
      throw std::runtime_error(
          "Linear dependent rows in A. This case is not implemented yet.");
    }

    return BFS<Field>{solution.point[{0, d}, 0], std::move(row_basis)};
  }

  // same as solve_from, but automatically finds bfs
  LPSolution<Field> solve() {
    auto bfs = find_bfs();

    if (!bfs.has_value()) {
      return NoFeasibleElements{};
    }

    return variant_cast<LPSolution<Field>>(solve_from(*bfs));
  }
};
