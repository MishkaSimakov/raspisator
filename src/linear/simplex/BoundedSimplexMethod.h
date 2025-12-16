#pragma once

#include <algorithm>
#include <fstream>
#include <random>
#include <ranges>
#include <variant>

#include "Settings.h"
#include "linear/matrix/LU.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/Rank.h"
#include "linear/matrix/RowBasis.h"
#include "linear/model/LP.h"
#include "linear/sparse/LU.h"
#include "utils/Accumulators.h"
#include "utils/Variant.h"

namespace simplex {

// Double, double toil and trouble;
// Fire burn and caldron bubble.
// - Macbeth
//
// Though this be madness, yet there is method in't.
// - Hamlet

// Solves cx -> max, Ax = b, l <= x <= u
// A is (n, d) matrix, b is (n, 1) matrix, c is (1, d) matrix
// l, u are (d, 1) matrices
// it is assumed that n < d

// This version is specifically tuned for sparse A matrix
template <typename Field>
class BoundedSimplexMethod {
  CSCMatrix<Field> A_;

  Matrix<Field> b_;
  Matrix<Field> c_;

  std::vector<VariableState> variables_;

  linalg::LUPA<Field> lupa_;

  CSCMatrix<Field> L_;
  CSCMatrix<Field> U_;
  std::vector<size_t> P_;

  Settings<Field> settings_;

  mutable std::mt19937 random_engine_;

  Matrix<Field> get_point(Matrix<Field> point, std::vector<size_t> basic_vars,
                          const std::vector<Field>& lower_bounds,
                          const std::vector<Field>& upper_bounds) {
    auto [n, d] = A_.shape();

    Matrix<Field> result(d, 1, 0);

    for (size_t i = 0; i < n; ++i) {
      result[basic_vars[i], 0] = point[i, 0];
    }
    for (size_t i = 0; i < d; ++i) {
      if (variables_[i] == VariableState::AT_LOWER) {
        result[i, 0] = lower_bounds[i];
      } else if (variables_[i] == VariableState::AT_UPPER) {
        result[i, 0] = upper_bounds[i];
      }
    }

    return result;
  }

  SimplexResult<Field> construct_finite_solution(
      Matrix<Field> point, std::vector<size_t> basic_vars,
      const std::vector<Field>& lower_bounds,
      const std::vector<Field>& upper_bounds, size_t iteration) {
    auto result = get_point(point, basic_vars, lower_bounds, upper_bounds);

    Field value = (c_ * result)[0, 0];

    FiniteLPSolution solution{std::move(result), std::move(basic_vars),
                              std::move(value), variables_};

    return SimplexResult<Field>{
        .iterations_count = iteration,
        .solution = solution,
    };
  }

  // largest violation variable selection
  // turned out to be prone to cycling when coupled with strong branching
  std::optional<std::pair<size_t, VariableState>>
  get_dual_leaving_variable_dantzig(
      const Matrix<Field>& point, const std::vector<size_t>& basic_vars,
      const std::vector<Field>& lower_bounds,
      const std::vector<Field>& upper_bounds) const {
    auto [n, d] = A_.shape();

    std::optional<std::pair<size_t, VariableState>> result = std::nullopt;
    Field largest_violation = 0;

    for (size_t i = 0; i < n; ++i) {
      Field lower_violation = lower_bounds[basic_vars[i]] - point[i, 0];
      Field upper_violation = point[i, 0] - upper_bounds[basic_vars[i]];

      if (FieldTraits<Field>::is_strictly_positive(lower_violation -
                                                   largest_violation)) {
        result = {i, VariableState::AT_LOWER};
        largest_violation = lower_violation;
      } else if (FieldTraits<Field>::is_strictly_positive(upper_violation -
                                                          largest_violation)) {
        result = {i, VariableState::AT_UPPER};
        largest_violation = upper_violation;
      }
    }

    return result;
  }

  // Although it seems reasonable to choose always variable with the largest
  // boundaries violation, this approach leads to cycling. To avoid cycling
  // Bland's rule is adopted.
  std::optional<std::pair<size_t, VariableState>>
  get_dual_leaving_variable_bland(
      const Matrix<Field>& point, const std::vector<size_t>& basic_vars,
      const std::vector<Field>& lower_bounds,
      const std::vector<Field>& upper_bounds) const {
    auto [n, d] = A_.shape();

    std::optional<std::pair<size_t, VariableState>> result = std::nullopt;
    size_t smallest_violating_index = 2 * d;

    for (size_t i = 0; i < n; ++i) {
      Field lower_violation = lower_bounds[basic_vars[i]] - point[i, 0];
      Field upper_violation = point[i, 0] - upper_bounds[basic_vars[i]];

      if (FieldTraits<Field>::is_strictly_positive(lower_violation)) {
        if (basic_vars[i] < smallest_violating_index) {
          result = {i, VariableState::AT_LOWER};
          smallest_violating_index = basic_vars[i];
        }
      } else if (FieldTraits<Field>::is_strictly_positive(upper_violation)) {
        if (basic_vars[i] + d < smallest_violating_index) {
          result = {i, VariableState::AT_UPPER};
          smallest_violating_index = basic_vars[i] + d;
        }
      }
    }

    return result;
  }

  // Select a random boundary violating variable. This strategy is used when
  // potential cycling is detected.
  std::optional<std::pair<size_t, VariableState>>
  get_dual_leaving_variable_randomly(
      const Matrix<Field>& point, const std::vector<size_t>& basic_vars,
      const std::vector<Field>& lower_bounds,
      const std::vector<Field>& upper_bounds) const {
    auto [n, d] = A_.shape();

    std::vector<std::pair<size_t, VariableState>> result;

    for (size_t i = 0; i < n; ++i) {
      Field lower_violation = lower_bounds[basic_vars[i]] - point[i, 0];
      Field upper_violation = point[i, 0] - upper_bounds[basic_vars[i]];

      if (FieldTraits<Field>::is_strictly_positive(lower_violation)) {
        result.emplace_back(i, VariableState::AT_LOWER);
      } else if (FieldTraits<Field>::is_strictly_positive(upper_violation)) {
        result.emplace_back(i, VariableState::AT_UPPER);
      }
    }

    if (result.empty()) {
      return std::nullopt;
    }

    std::uniform_int_distribution<size_t> index_dist(0, result.size() - 1);
    size_t index = index_dist(random_engine_);

    return result.at(index);
  }

  Matrix<Field> get_reduced_costs(const CSCMatrix<Field>& L,
                                  const CSCMatrix<Field>& U,
                                  const std::vector<size_t>& P,
                                  const std::vector<size_t>& basic_vars) {
    auto [n, d] = A_.shape();

    Matrix<Field> result(d, 1);

    Matrix<Field> cb(n, 1);
    for (size_t i = 0; i < n; ++i) {
      cb[i, 0] = c_[0, basic_vars[i]];
    }
    linalg::solve_transposed_linear_inplace(L, U, P, cb);

    for (size_t i = 0; i < d; ++i) {
      if (variables_[i] == VariableState::BASIC) {
        result[i, 0] = 0;
        continue;
      }

      result[i, 0] = c_[0, i];

      for (const auto& [row, value] : A_.get_column(i)) {
        result[i, 0] -= value * cb[P[row], 0];
      }
    }

    return result;
  }

  std::optional<size_t> get_dual_entering_variable(
      const CSCMatrix<Field>& L, const CSCMatrix<Field>& U,
      const std::vector<size_t>& P, size_t entering_var,
      VariableState leaving_state, const std::vector<size_t>& basic_vars) {
    auto [n, d] = A_.shape();

    ArgMinimum<Field> min_ratio;

    auto reduced_costs = get_reduced_costs(L, U, P, basic_vars);

    Matrix<Field> e(n, 1, 0);
    e[entering_var, 0] = 1;
    linalg::solve_transposed_linear_inplace(L, U, P, e);

    for (size_t i = 0; i < d; ++i) {
      if (variables_[i] == VariableState::BASIC) {
        continue;
      }

      Field coef = 0;
      for (const auto& [row, value] : A_.get_column(i)) {
        coef += e[P[row], 0] * value;
      }

      if (!FieldTraits<Field>::is_nonzero(coef)) {
        continue;
      }

      // if (!((variables_[i] == VariableState::AT_LOWER &&
      //        !FieldTraits<Field>::is_strictly_positive(reduced_costs[i, 0]))
      //        ||
      //       (variables_[i] == VariableState::AT_UPPER &&
      //        !FieldTraits<Field>::is_strictly_negative(reduced_costs[i,
      //        0])))) {
      //   throw std::runtime_error(
      //       std::format("Current point is not dual feasible! Reduced cost for
      //       "
      //                   "variable #{} has value {}.",
      //                   i, reduced_costs[i, 0]));
      // }
      if (!((variables_[i] == VariableState::AT_LOWER &&
             reduced_costs[i, 0] < Field(1) / Field(1e5)) ||
            (variables_[i] == VariableState::AT_UPPER &&
             reduced_costs[i, 0] > -Field(1) / Field(1e5)))) {
        throw std::runtime_error(
            std::format("Current point is not dual feasible! Reduced cost for "
                        "variable #{} has value {}.",
                        i, reduced_costs[i, 0]));
      }

      Field ratio = reduced_costs[i, 0] / coef;

      if (leaving_state == VariableState::AT_UPPER) {
        ratio *= -1;
      }

      if (leaving_state == VariableState::AT_LOWER) {
        if (variables_[i] == VariableState::AT_LOWER && coef > 0 ||
            variables_[i] == VariableState::AT_UPPER && coef < 0) {
          continue;
        }
      } else {
        if (variables_[i] == VariableState::AT_LOWER && coef < 0 ||
            variables_[i] == VariableState::AT_UPPER && coef > 0) {
          continue;
        }
      }

      min_ratio.record(i, ratio);
    }

    return min_ratio.argmin();
  }

 public:
  BoundedSimplexMethod(CSCMatrix<Field> A, Matrix<Field> b, Matrix<Field> c,
                       Settings<Field> settings = {})
      : A_(std::move(A)),
        b_(std::move(b)),
        c_(std::move(c)),
        variables_(A_.shape().second),
        lupa_(A_),
        L_(A_.shape().first),
        U_(A_.shape().first),
        P_(A.shape().first),
        settings_(settings) {
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

  void dump_state(std::ostream& os, size_t dump_id,
                  const std::vector<Field>& lower,
                  const std::vector<Field>& upper) {
    os << "namespace SimplexDump_" << dump_id << " {\n";

    os << "Matrix<Field> A = {" << linalg::to_dense(A_) << "};\n";
    os << "Matrix<Field> b = {" << b_ << "};\n";
    os << "Matrix<Field> c = {" << c_ << "};\n";

    os << "std::vector<Field> upper = {";
    for (auto value : upper) {
      os << value << ", ";
    }
    os << "};\n";

    os << "std::vector<Field> lower = {";
    for (auto value : lower) {
      os << value << ", ";
    }
    os << "};\n";

    os << "std::vector<VariableState> states = {";
    for (auto state : variables_) {
      if (state == VariableState::BASIC) {
        os << "VariableState::BASIC, ";
      } else if (state == VariableState::AT_LOWER) {
        os << "VariableState::AT_LOWER, ";
      } else {
        os << "VariableState::AT_UPPER, ";
      }
    }
    os << "};\n";
    os << "}\n";

    os << std::flush;
  }

  const std::vector<VariableState>& get_variables_states() const {
    return variables_;
  }

  void setup_warm_start(const std::vector<size_t>& basic_variables) {
    auto [n, d] = A_.shape();

    if (basic_variables.size() != n) {
      throw std::invalid_argument("Wrong basic variables count.");
    }

    lupa_.get_lup(basic_variables, L_, U_, P_);

    auto reduced_costs = get_reduced_costs(L_, U_, P_, basic_variables);

    for (size_t i = 0; i < d; ++i) {
      if (reduced_costs[i, 0] <= 0) {
        variables_[i] = VariableState::AT_LOWER;
      } else {
        variables_[i] = VariableState::AT_UPPER;
      }
    }

    for (size_t basic_var : basic_variables) {
      variables_[basic_var] = VariableState::BASIC;
    }
  }

  void setup_warm_start(const std::vector<VariableState>& variables) {
    variables_ = variables;
  }

  void set_max_iterations(std::optional<size_t> max_iterations) {
    settings_.max_iterations = max_iterations;
  }

  // solves the problem using dual bounded simplex method
  SimplexResult<Field> dual(const std::vector<Field>& lower_bounds,
                            const std::vector<Field>& upper_bounds) {
    auto [n, d] = A_.shape();

    std::vector<size_t> basic_vars;
    for (size_t i = 0; i < variables_.size(); ++i) {
      if (variables_[i] == VariableState::BASIC) {
        basic_vars.push_back(i);
      }
    }

    if (basic_vars.size() != n) {
      throw std::invalid_argument("Wrong number of basic variables.");
    }

    size_t iteration = 0;
    size_t anticycling_iterations = 0;

    while (true) {
      if ((iteration + 1) % 10'000 == 0) {
        // {
        //   size_t since_epoch =
        //       std::chrono::system_clock::now().time_since_epoch() /
        //       std::chrono::milliseconds(1);
        //
        //   std::ofstream os(std::format("simplex_core_dump_{}.h",
        //   since_epoch)); dump_state(os, since_epoch, lower_bounds,
        //   upper_bounds);
        // }
        //
        // throw std::runtime_error("Cycling!");

        std::cout << "Cycling!" << std::endl;
        anticycling_iterations = 5;
      }

      lupa_.get_lup(basic_vars, L_, U_, P_);

      for (size_t i : P_) {
        if (i >= L_.shape().first) {
          std::ofstream os("simplex_core_dump_4.h");
          dump_state(os, 4, lower_bounds, upper_bounds);
          throw std::runtime_error("Wrong!!!");
        }
      }

      Matrix<Field> b(b_);
      for (size_t col = 0; col < d; ++col) {
        if (variables_[col] == VariableState::AT_LOWER) {
          for (const auto& [row, value] : A_.get_column(col)) {
            b[row, 0] -= value * lower_bounds[col];
          }
        } else if (variables_[col] == VariableState::AT_UPPER) {
          for (const auto& [row, value] : A_.get_column(col)) {
            b[row, 0] -= value * upper_bounds[col];
          }
        }
      }
      auto point = linalg::solve_linear(L_, U_, P_, b);

      if (settings_.max_iterations && iteration >= settings_.max_iterations) {
        auto full_point =
            get_point(point, basic_vars, lower_bounds, upper_bounds);
        Field value = (c_ * full_point)[0, 0];

        return SimplexResult<Field>{
            .iterations_count = iteration,
            .solution = ReachedIterationsLimit{value},
        };
      }

      std::optional<std::pair<size_t, VariableState>> leaving;
      if (anticycling_iterations > 0) {
        leaving = get_dual_leaving_variable_randomly(
            point, basic_vars, lower_bounds, upper_bounds);

        --anticycling_iterations;
      } else {
        leaving = get_dual_leaving_variable_dantzig(point, basic_vars,
                                                    lower_bounds, upper_bounds);
      }

      if (!leaving.has_value()) {
        return construct_finite_solution(std::move(point),
                                         std::move(basic_vars), lower_bounds,
                                         upper_bounds, iteration);
      }

      auto entering = get_dual_entering_variable(L_, U_, P_, leaving->first,
                                                 leaving->second, basic_vars);

      if (!entering.has_value()) {
        return SimplexResult<Field>{
            .iterations_count = iteration,
            .solution = NoFeasibleElements{},
        };
      }

      // pivot
      variables_[*entering] = VariableState::BASIC;
      variables_[basic_vars[leaving->first]] = leaving->second;
      basic_vars[leaving->first] = *entering;

      ++iteration;
    }
  }
};

}  // namespace simplex
