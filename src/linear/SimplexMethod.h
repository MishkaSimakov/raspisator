#pragma once
#include <variant>

#include "linear/model/LP.h"
#include "matrix/LU.h"
#include "matrix/Matrix.h"
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

// Double, double toil and trouble;
// Fire burn and caldron bubble.
// - Macbeth
//
// There are a lot of problems with a double. This version of SimplexMethod
// is trying to be numerically stable.

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

  Field estimate_rounding_errors(const Matrix<Field>& tableau,
                                 const std::vector<size_t>& basic_vars) const {
    auto [n, d] = A_.shape();

    auto x = Matrix<Field>(d, 1);
    for (size_t i = 0; i < n; ++i) {
      x[basic_vars[i], 0] = tableau[i, 0];
    }

    auto error = A_ * x - b_;

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

      auto change = linalg::dot(dual_point, A_[{0, n}, i]) - c_[0, i];

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

    size_t n = A_.get_height();

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

  auto get_basic_lup(const std::vector<size_t>& basic_vars) {
    size_t n = basic_vars.size();
    Matrix<Field> temp(n, n);

    for (size_t i = 0; i < n; ++i) {
      temp[{0, n}, i] = A_[{0, n}, basic_vars[i]];
    }

    return linalg::get_lup(temp);
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

  auto get_entering_coordinates(const Matrix<Field>& L, const Matrix<Field>& U,
                                const std::vector<size_t>& P,
                                size_t entering_var) {
    auto Pb =
        linalg::apply_permutation(A_[{0, A_.get_height()}, entering_var], P);

    // solve Ly = Pb
    auto y = linalg::solve_lower(L, Pb, std::true_type{});

    // solve Ux = y
    auto x = linalg::solve_upper(U, y, std::false_type{});

    return x;
  }

  auto get_dual_point(const Matrix<Field>& L, const Matrix<Field>& U,
                      const std::vector<size_t>& P,
                      const std::vector<size_t>& basic_vars) {
    size_t n = A_.get_height();

    Matrix<Field> cb(n, 1);
    for (size_t i = 0; i < n; ++i) {
      cb[i, 0] = c_[0, basic_vars[i]];
    }

    // TODO: lots of unnecessary copying here
    auto z = linalg::solve_lower(linalg::transposed(U), cb, std::false_type{});
    auto y = linalg::solve_upper(linalg::transposed(L), z, std::true_type{});

    std::vector<size_t> transposed_P(n);
    for (size_t i = 0; i < n; ++i) {
      transposed_P[P[i]] = i;
    }

    return linalg::apply_permutation(y, P);
  }

  // finds maximum starting from bfs (basic feasible solution)
  std::variant<FiniteLPSolution<Field>, InfiniteSolution> solve_from(
      BFS<Field> bfs) {
    auto [n, d] = A_.shape();

    if (bfs.point.shape() != std::pair{d, 1}) {
      throw std::invalid_argument("bfs has wrong shape.");
    }

    while (true) {
      auto [L, U, P] = get_basic_lup(bfs.basic_variables);

      std::cout << "basic variables" << std::endl;
      for (auto i : bfs.basic_variables) {
        std::cout << i << " ";
      }
      std::cout << std::endl;
      //
      // std::cout << "L\n" << L << "\nU\n" << U << std::endl;

      // obtain point associated with given basic variables by solving
      // Bu = b
      auto point = get_point(L, U, P);

      // obtain a solution of a dual problem by solving
      // B^T u = c_b
      auto dual_point = get_dual_point(L, U, P, bfs.basic_variables);

      auto [entering_var, change] =
          find_entering_variable(dual_point, bfs.basic_variables);

      std::cout << change << std::endl;

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

      // TODO: pivoting
      // for (size_t j = 0; j < d; ++j) {
      //   if (solution.tableau[i, j + 1] != 0) {
      //     pivot(solution.tableau, solution.basic_variables, j, i);
      //     pivoted = true;
      //   }
      // }

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
  LPSolution<Field> solve() {
    auto bfs = find_bfs();

    std::cout << "found bfs" << std::endl;

    if (!bfs.has_value()) {
      return NoFeasibleElements{};
    }

    return variant_cast<LPSolution<Field>>(solve_from(*bfs));
  }
};
