#pragma once

#include <algorithm>
#include <chrono>
#include <fstream>
#include <random>
#include <ranges>
#include <unordered_map>
#include <variant>

#include "CyclingDetector.h"
#include "Settings.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/NPY.h"
#include "linear/matrix/RowBasis.h"
#include "linear/model/LP.h"
#include "linear/sparse/LU.h"
#include "utils/Accumulators.h"
#include "utils/Hashers.h"

namespace simplex {

// Double, double toil and trouble;
// Fire burn and caldron bubble.
// - from Macbeth
//
// Though this be madness, yet there is method in't.
// - from Hamlet

// Solves cx -> max, Ax = b, l <= x <= u
// A is (n, d) matrix, b is (n, 1) matrix, c is (1, d) matrix
// l and u are (d, 1) matrices (possibly with infinite values)
// it is assumed that n < d

// This version is specifically tuned for sparse A matrix
template <typename Field>
class BoundedSimplexMethod {
  CSCMatrix<Field> A_;

  Matrix<Field> b_;
  Matrix<Field> c_;

  std::vector<VariableState> variables_;

  linalg::LUPA<Field> lupa_;

  Settings<Field> settings_;

  mutable std::mt19937 random_engine_;

  CyclingDetector<Field> cycling_;

  Matrix<Field> get_point(Matrix<Field> point, std::vector<size_t> basic_vars,
                          const Bounds<Field>& bounds) {
    auto [n, d] = A_.shape();

    Matrix<Field> result(d, 1, 0);

    for (size_t i = 0; i < n; ++i) {
      result[basic_vars[i], 0] = point[i, 0];
    }
    for (size_t i = 0; i < d; ++i) {
      if (variables_[i] == VariableState::AT_LOWER) {
        result[i, 0] = *bounds[i].lower;
      } else if (variables_[i] == VariableState::AT_UPPER) {
        result[i, 0] = *bounds[i].upper;
      }
    }

    return result;
  }

  SimplexResult<Field> construct_finite_solution(Matrix<Field> point,
                                                 std::vector<size_t> basic_vars,
                                                 const Bounds<Field>& bounds,
                                                 size_t iteration) {
    auto result = get_point(point, basic_vars, bounds);

    Field value = (c_ * result)[0, 0];

    FiniteLPSolution solution{std::move(result), std::move(value), variables_};

    return SimplexResult<Field>{
        .iterations_count = iteration,
        .solution = solution,
    };
  }

  // largest violation variable selection
  // turned out to be prone to cycling when coupled with strong branching
  std::optional<std::pair<size_t, VariableState>>
  get_dual_leaving_variable_dantzig(const Matrix<Field>& point,
                                    const std::vector<size_t>& basic_vars,
                                    const Bounds<Field>& bounds) const {
    auto [n, d] = A_.shape();

    ArgMaximum<BoundViolation<Field>, std::less<BoundViolation<Field>>>
        max_violation;

    for (size_t i = 0; i < n; ++i) {
      max_violation.record(i, bounds[basic_vars[i]].get_violation(point[i, 0]));
    }

    if (!max_violation.max()) {
      return std::nullopt;
    }

    switch (max_violation.max()->type) {
      case BoundViolationType::NONE:
        return std::nullopt;
      case BoundViolationType::VIOLATE_LOWER_BOUND:
        return std::pair{*max_violation.argmax(), VariableState::AT_LOWER};
      case BoundViolationType::VIOLATE_UPPER_BOUND:
        return std::pair{*max_violation.argmax(), VariableState::AT_UPPER};
      default:
        std::unreachable();
    }
  }

  // Select a random boundary violating variable. This strategy is used when
  // potential cycling is detected.
  std::optional<std::pair<size_t, VariableState>>
  get_dual_leaving_variable_randomly(const Matrix<Field>& point,
                                     const std::vector<size_t>& basic_vars,
                                     const Bounds<Field>& bounds) const {
    auto [n, d] = A_.shape();

    std::vector<std::pair<size_t, VariableState>> result;

    for (size_t i = 0; i < n; ++i) {
      auto violation = bounds[basic_vars[i]].get_violation(point[i, 0]);

      if (violation.type == BoundViolationType::VIOLATE_LOWER_BOUND) {
        result.emplace_back(i, VariableState::AT_LOWER);
      } else if (violation.type == BoundViolationType::VIOLATE_UPPER_BOUND) {
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

  Matrix<Field> get_reduced_costs(const std::vector<size_t>& basic_vars) {
    auto [n, d] = A_.shape();

    Matrix<Field> cb(n, 1);
    for (size_t i = 0; i < n; ++i) {
      cb[i, 0] = c_[0, basic_vars[i]];
    }
    cb = lupa_.solve_linear_transposed(std::move(cb));

    Matrix<Field> result(d, 1);
    for (size_t i = 0; i < d; ++i) {
      result[i, 0] = c_[0, i];

      for (const auto& [row, value] : A_.get_column(i)) {
        result[i, 0] -= value * cb[row, 0];
      }
    }

    return result;
  }

  std::optional<size_t> get_dual_entering_variable(
      size_t entering_var, VariableState leaving_state,
      const std::vector<size_t>& basic_vars) {
    auto [n, d] = A_.shape();

    ArgMinimum<Field> min_ratio;

    auto reduced_costs = get_reduced_costs(basic_vars);

    auto inverse_row = lupa_.get_row(entering_var);

    for (size_t i = 0; i < d; ++i) {
      if (variables_[i] == VariableState::BASIC) {
        continue;
      }

      Field coef = 0;
      for (const auto& [row, value] : A_.get_column(i)) {
        if (FieldTraits<Field>::is_nonzero(inverse_row[row, 0])) {
          coef += inverse_row[row, 0] * value;
        }
      }

      if (!FieldTraits<Field>::is_nonzero(coef)) {
        continue;
      }

      if (!FieldTraits<Field>::is_nonzero(reduced_costs[i, 0])) {
        reduced_costs[i, 0] = 0;
      }

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
        if (variables_[i] == VariableState::AT_LOWER && coef > Field(0) ||
            variables_[i] == VariableState::AT_UPPER && coef < Field(0)) {
          continue;
        }
      } else {
        if (variables_[i] == VariableState::AT_LOWER && coef < Field(0) ||
            variables_[i] == VariableState::AT_UPPER && coef > Field(0)) {
          continue;
        }
      }

      min_ratio.record(i, ratio);

      // if (*min_ratio.argmin() == i) {
      //   std::cout << reduced_costs[i, 0] << "/" << coef << " ";
      // }
    }

    // std::cout << std::endl;
    // std::cout << "min ratio: " << *min_ratio.min() << std::endl;

    return min_ratio.argmin();
  }

  SimplexResult<Field> dual_implementation(const Bounds<Field>& bounds) {
    auto [n, d] = A_.shape();

    std::vector<size_t> basic_vars;
    for (size_t i = 0; i < variables_.size(); ++i) {
      if (variables_[i] == VariableState::BASIC) {
        basic_vars.push_back(i);
      } else if (variables_[i] == VariableState::AT_LOWER && !bounds[i].lower ||
                 variables_[i] == VariableState::AT_UPPER && !bounds[i].upper) {
        return SimplexResult<Field>{
            .iterations_count = 0,
            .solution = Unbounded{},
        };
      }
    }

    cycling_.clear();
    lupa_.set_columns(basic_vars);

    if (basic_vars.size() != n) {
      throw std::invalid_argument("Wrong number of basic variables.");
    }

    size_t iteration = 0;
    size_t anticycling_iterations = 0;

    auto last_time = std::chrono::high_resolution_clock::now();
    size_t iterations_since_last_time = 0;

    while (true) {
      // std::cout << iteration << std::endl;

      {
        // Matrix<Field> B(n, n);
        // for (size_t i = 0; i < n; ++i) {
        //   for (auto [row, value] : A_.get_column(basic_vars[i])) {
        //     B[row, i] = value;
        //   }
        // }
        //
        // logging::log_npy(B, std::format("iterations/{}.npy", iteration));
        // logging::log_npy(lupa_.get_inverse(),
        //                  std::format("iterations/{}_inv.npy", iteration));
      }

      Field objective = 0;
      Matrix<Field> b(b_);
      for (size_t col = 0; col < d; ++col) {
        if (variables_[col] == VariableState::AT_LOWER) {
          for (const auto& [row, value] : A_.get_column(col)) {
            b[row, 0] -= value * *bounds[col].lower;
          }

          objective += *bounds[col].lower * c_[0, col];
        } else if (variables_[col] == VariableState::AT_UPPER) {
          for (const auto& [row, value] : A_.get_column(col)) {
            b[row, 0] -= value * *bounds[col].upper;
          }

          objective += *bounds[col].upper * c_[0, col];
        }
      }

      auto point = lupa_.solve_linear(b);

      for (size_t i = 0; i < n; ++i) {
        objective += point[i, 0] * c_[0, basic_vars[i]];
      }

      if (cycling_.record(iteration, basic_vars, objective) ==
          CyclingState::HAS_CYCLING) {
        anticycling_iterations = 5;
      }

      ++iterations_since_last_time;
      auto curr_time = std::chrono::high_resolution_clock::now();
      if (curr_time - last_time > std::chrono::seconds{1}) {
        double speed =
            static_cast<double>(iterations_since_last_time) /
            std::chrono::duration<double>(curr_time - last_time).count();

        std::println("{:.1f} itr/s, objective: {}, size: {}, elapsed: {}",
                     speed, objective, lupa_.size(), curr_time - last_time);

        last_time = curr_time;
        iterations_since_last_time = 0;
      }

      if (settings_.max_iterations && iteration >= settings_.max_iterations) {
        return SimplexResult<Field>{
            .iterations_count = iteration,
            .solution = ReachedIterationsLimit{objective},
        };
      }

      std::optional<std::pair<size_t, VariableState>> leaving;
      if (anticycling_iterations > 0) {
        leaving = get_dual_leaving_variable_randomly(point, basic_vars, bounds);

        --anticycling_iterations;
      } else {
        leaving = get_dual_leaving_variable_dantzig(point, basic_vars, bounds);
      }

      if (!leaving.has_value()) {
        return construct_finite_solution(
            std::move(point), std::move(basic_vars), bounds, iteration);
      }

      auto entering = get_dual_entering_variable(leaving->first,
                                                 leaving->second, basic_vars);

      if (!entering.has_value()) {
        return SimplexResult<Field>{
            .iterations_count = iteration,
            .solution = NoFeasibleElements{},
        };
      }

      // pivot
      lupa_.change_column(leaving->first, *entering);

      // std::println("{} -> {}", basic_vars[leaving->first], *entering);

      variables_[*entering] = VariableState::BASIC;
      variables_[basic_vars[leaving->first]] = leaving->second;
      basic_vars[leaving->first] = *entering;

      ++iteration;
    }
  }

  void dump_state(const Bounds<Field>& bounds,
                  const std::vector<VariableState>& initial_variables) {
    size_t dump_id = std::chrono::system_clock::now().time_since_epoch() /
                     std::chrono::milliseconds(1);
    std::string dump_name = std::format("simplex_core_dump_{}.h", dump_id);

    std::ofstream os(dump_name);

    os << "namespace SimplexDump_" << dump_id << " {\n";

    os << "Matrix<Field> A = {" << linalg::to_dense(A_) << "};\n";
    os << "Matrix<Field> b = {" << b_ << "};\n";
    os << "Matrix<Field> c = {" << c_ << "};\n";

    // os << "std::vector<Field> upper = {";
    // for (auto value : upper) {
    //   os << value << ", ";
    // }
    // os << "};\n";
    //
    // os << "std::vector<Field> lower = {";
    // for (auto value : lower) {
    //   os << value << ", ";
    // }
    // os << "};\n";

    os << "std::vector<VariableState> last_states = {";
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

    os << "std::vector<VariableState> init_states = {";
    for (auto state : initial_variables) {
      if (state == VariableState::BASIC) {
        os << "VariableState::BASIC, ";
      } else if (state == VariableState::AT_LOWER) {
        os << "VariableState::AT_LOWER, ";
      } else {
        os << "VariableState::AT_UPPER, ";
      }
    }
    os << "};";

    os << "}\n";

    os << std::flush;

    std::println("Registered failed simplex run into {}.", dump_name);

    if constexpr (std::same_as<Field, double>) {
      std::ofstream inv_os("b_inverse.npy");
      linalg::to_npy(inv_os, lupa_.get_inverse());
    }
  }

 public:
  BoundedSimplexMethod(CSCMatrix<Field> A, Matrix<Field> b, Matrix<Field> c,
                       Settings<Field> settings = {})
      : A_(std::move(A)),
        b_(std::move(b)),
        c_(std::move(c)),
        variables_(A_.shape().second),
        lupa_(A_),
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

  void set_max_iterations(std::optional<size_t> max_iterations) {
    settings_.max_iterations = max_iterations;
  }

  std::vector<VariableState> get_initial_states(
      const std::vector<size_t>& basic_variables) {
    auto [n, d] = A_.shape();

    if (basic_variables.size() != n) {
      throw std::invalid_argument(std::format(
          "Wrong basic variables count: {} != {}", basic_variables.size(), n));
    }

    std::vector<VariableState> states(d);

    lupa_.set_columns(basic_variables);

    auto reduced_costs = get_reduced_costs(basic_variables);

    std::cout << reduced_costs << std::endl;
    for (size_t i = 0; i < d; ++i) {
      if (reduced_costs[i, 0] <= Field(0)) {
        states[i] = VariableState::AT_LOWER;
      } else {
        states[i] = VariableState::AT_UPPER;
      }
    }

    for (size_t basic_var : basic_variables) {
      states[basic_var] = VariableState::BASIC;
    }

    return states;
  }

  std::vector<VariableState> get_initial_states() {
    auto basic_variables =
        linalg::get_row_basis(linalg::transposed(linalg::to_dense(A_)));

    return get_initial_states(basic_variables);
  }

  SimplexResult<Field> dual(
      const Bounds<Field>& bounds,
      const std::vector<VariableState>& initial_variables) {
    variables_ = initial_variables;

    try {
      return dual_implementation(bounds);
    } catch (...) {
      dump_state(bounds, initial_variables);
      throw;
    }
  }

  SimplexResult<Field> dual(
      const Bounds<Field>& bounds,
      const std::vector<size_t>& initial_basic_variables) {
    auto initial_variables = get_initial_states(initial_basic_variables);
    variables_ = initial_variables;

    try {
      return dual_implementation(bounds);
    } catch (...) {
      dump_state(bounds, initial_variables);
      throw;
    }
  }

  SimplexResult<Field> dual(const Bounds<Field>& bounds) {
    auto initial_variables = get_initial_states();
    variables_ = initial_variables;

    try {
      return dual_implementation(bounds);
    } catch (...) {
      dump_state(bounds, initial_variables);
      throw;
    }
  }
};

}  // namespace simplex
