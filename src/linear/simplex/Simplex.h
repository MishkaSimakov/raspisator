#pragma once

#include <algorithm>
#include <chrono>
#include <fstream>
#include <random>
#include <ranges>
#include <unordered_map>
#include <variant>

#include "Accountant.h"
#include "CyclingDetector.h"
#include "Dual.h"
#include "Settings.h"
#include "SimplexCoreDump.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/NPY.h"
#include "linear/matrix/RowBasis.h"
#include "linear/model/LP.h"
#include "linear/sparse/LU.h"
#include "utils/Accumulators.h"

namespace simplex {

template <typename Field>
struct Tolerances {
  Field feasibility{0};
  Field pivot{0};
};

template <>
struct Tolerances<double> {
  double feasibility{1e-7};
  double pivot{1e-7};
};

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
template <typename Field, typename Accountant = EmptyAccountant<Field>>
class Simplex {
  CSCMatrix<Field> A_;
  Matrix<Field> b_;
  Matrix<Field> c_;

  IterationState<Field> state_;

  Settings<Field> settings_;
  Accountant accountant_;

  Tolerances<Field> tolerances_;

  static Matrix<Field> get_point(const IterationState<Field>& state) {
    auto [n, d] = state.problem_shape();

    Matrix<Field> result(d, 1, 0);

    for (size_t i = 0; i < n; ++i) {
      result[state.basic_variables[i], 0] = state.basic_point[i, 0];
    }
    for (size_t i = 0; i < d; ++i) {
      if (state.variables_states[i] == VariableState::AT_LOWER) {
        result[i, 0] = *(*state.bounds)[i].lower;
      } else if (state.variables_states[i] == VariableState::AT_UPPER) {
        result[i, 0] = *(*state.bounds)[i].upper;
      }
    }

    return result;
  }

  template <typename T>
  static SimplexResult<Field> construct_result(
      const IterationState<Field>& state) {
    auto point = get_point(state);

    if constexpr (std::same_as<T, FiniteLPSolution<Field>>) {
      return SimplexResult<Field>{
          .iterations_count = state.iteration_index,
          .solution = FiniteLPSolution{std::move(point), state.objective,
                                       state.variables_states},
      };
    }
    if constexpr (std::same_as<T, NoFeasibleElements>) {
      return SimplexResult<Field>{
          .iterations_count = state.iteration_index,
          .solution = NoFeasibleElements{},
      };
    }
    if constexpr (std::same_as<T, ReachedIterationsLimit<Field>>) {
      return SimplexResult<Field>{
          .iterations_count = state.iteration_index,
          .solution = ReachedIterationsLimit<Field>{state.objective},
      };
    }
    if constexpr (std::same_as<T, Unbounded>) {
      return SimplexResult<Field>{
          .iterations_count = state.iteration_index,
          .solution = Unbounded{},
      };
    }
  }

  Matrix<Field> get_basic_costs(const std::vector<size_t>& basic_vars) const {
    auto [n, d] = A_.shape();

    Matrix<Field> result(n, 1);
    for (size_t i = 0; i < n; ++i) {
      result[i, 0] = c_[0, basic_vars[i]];
    }

    return result;
  }

  Matrix<Field> get_reduced_costs(Matrix<Field> basic_costs) const {
    auto [n, d] = A_.shape();

    basic_costs = state_.lupa.solve_linear_transposed(std::move(basic_costs));

    Matrix<Field> result(d, 1);
    for (size_t i = 0; i < d; ++i) {
      result[i, 0] = c_[0, i];

      for (const auto& [row, value] : A_.get_column(i)) {
        result[i, 0] -= value * basic_costs[row, 0];
      }
    }

    return result;
  }

  std::optional<size_t> get_dual_entering_variable(
      LeavingVariable leaving, const IterationState<Field>& state) const {
    auto [n, d] = state.problem_shape();

    ArgMinimum<Field> min_ratio;

    auto reduced_costs =
        get_reduced_costs(get_basic_costs(state.basic_variables));

    auto inverse_row = state.lupa.get_row(leaving.index);

    for (size_t i = 0; i < d; ++i) {
      if (state.variables_states[i] == VariableState::BASIC) {
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

      if (!((state.variables_states[i] == VariableState::AT_LOWER &&
             reduced_costs[i, 0] < Field(1) / Field(1e5)) ||
            (state.variables_states[i] == VariableState::AT_UPPER &&
             reduced_costs[i, 0] > -Field(1) / Field(1e5)))) {
        throw std::runtime_error(
            std::format("Current point is not dual feasible! Reduced cost for "
                        "variable #{} has value {}.",
                        i, reduced_costs[i, 0]));
      }

      Field ratio = reduced_costs[i, 0] / coef;

      if (leaving.new_state == VariableState::AT_UPPER) {
        ratio *= -1;
      }

      if (leaving.new_state == VariableState::AT_LOWER) {
        if (state.variables_states[i] == VariableState::AT_LOWER &&
                coef > Field(0) ||
            state.variables_states[i] == VariableState::AT_UPPER &&
                coef < Field(0)) {
          continue;
        }
      } else {
        if (state.variables_states[i] == VariableState::AT_LOWER &&
                coef < Field(0) ||
            state.variables_states[i] == VariableState::AT_UPPER &&
                coef > Field(0)) {
          continue;
        }
      }

      min_ratio.record(i, ratio);
    }

    return min_ratio.argmin();
  }

  Matrix<Field> get_rhs(const Bounds<Field>& bounds,
                        const std::vector<VariableState>& states) const {
    auto [n, d] = A_.shape();

    Matrix<Field> rhs(b_);
    for (size_t col = 0; col < d; ++col) {
      if (states[col] == VariableState::AT_LOWER) {
        for (const auto& [row, value] : A_.get_column(col)) {
          rhs[row, 0] -= value * *bounds[col].lower;
        }
      } else if (states[col] == VariableState::AT_UPPER) {
        for (const auto& [row, value] : A_.get_column(col)) {
          rhs[row, 0] -= value * *bounds[col].upper;
        }
      }
    }

    return std::move(rhs);
  }

  Field get_objective(const IterationState<Field>& state) const {
    auto [n, d] = A_.shape();
    KahanSum<Field> objective;

    for (size_t col = 0; col < d; ++col) {
      if (state.variables_states[col] == VariableState::AT_LOWER) {
        objective.add(c_[0, col] * *(*state.bounds)[col].lower);
      } else if (state.variables_states[col] == VariableState::AT_UPPER) {
        objective.add(c_[0, col] * *(*state.bounds)[col].upper);
      }
    }

    for (size_t i = 0; i < state.basic_variables.size(); ++i) {
      objective.add(c_[0, state.basic_variables[i]] * state.basic_point[i, 0]);
    }

    return objective.sum();
  }

  bool should_stop(const IterationState<Field>& state) const {
    return settings_.max_iterations &&
           state.iteration_index >= settings_.max_iterations;
  }

  void initialize_state(const Bounds<Field>& bounds,
                        const std::vector<VariableState>& states) {
    auto [n, d] = A_.shape();

    state_.iteration_index = 0;
    state_.last_cycling_iteration = std::nullopt;

    state_.basic_variables.clear();
    state_.variables_states = states;
    state_.bounds = &bounds;

    for (size_t i = 0; i < d; ++i) {
      if (states[i] == VariableState::BASIC) {
        state_.basic_variables.push_back(i);
      } else if (states[i] == VariableState::AT_LOWER && !bounds[i].lower ||
                 states[i] == VariableState::AT_UPPER && !bounds[i].upper) {
        throw std::invalid_argument(
            "Given point is not valid. Variable is set to boundary, but does "
            "not have it.");
      }
    }

    if (state_.basic_variables.size() != n) {
      throw std::invalid_argument("Wrong number of basic variables.");
    }

    state_.cycling.clear();
    state_.lupa.set_columns(state_.basic_variables);
  }

  SimplexResult<Field> dual_implementation(
      const Bounds<Field>& bounds, const std::vector<VariableState>& states) {
    auto [n, d] = A_.shape();

    if (settings_.is_strict) {
      if (!is_dual_feasible(bounds, states)) {
        throw std::invalid_argument(
            "Given initial state is not dual feasible.");
      }
    }

    // initialize simplex state
    initialize_state(bounds, states);

    DualLeavingVariable<Field> leaving_finder;

    while (true) {
      auto rhs = get_rhs(bounds, state_.variables_states);
      state_.basic_point = state_.lupa.solve_linear(rhs);

      state_.objective = get_objective(state_);

      if (state_.cycling.record(state_.iteration_index, state_.variables_states,
                                state_.objective) ==
          CyclingState::HAS_CYCLING) {
        state_.last_cycling_iteration = state_.iteration_index;
      }

      accountant_.iteration(state_);

      if (should_stop(state_)) {
        return construct_result<ReachedIterationsLimit<Field>>(state_);
      }

      auto leaving = leaving_finder.get(state_);
      if (!leaving) {
        return construct_result<FiniteLPSolution<Field>>(state_);
      }

      auto entering = get_dual_entering_variable(*leaving, state_);
      if (!entering) {
        return construct_result<NoFeasibleElements>(state_);
      }

      // pivot
      state_.lupa.change_column(leaving->index, *entering);

      state_.variables_states[*entering] = VariableState::BASIC;
      state_.variables_states[state_.basic_variables[leaving->index]] =
          leaving->new_state;
      state_.basic_variables[leaving->index] = *entering;

      ++state_.iteration_index;
    }
  }

  std::optional<size_t> get_primal_entering_variable(
      const IterationState<Field>& state) {
    auto [n, d] = state.problem_shape();

    auto reduced_costs =
        get_reduced_costs(get_basic_costs(state.basic_variables));

    ArgMaximum<Field> max_cost;

    for (size_t i = 0; i < d; ++i) {
      if (state.variables_states[i] == VariableState::BASIC) {
        continue;
      }

      Field cost = reduced_costs[i, 0];
      if (state.variables_states[i] == VariableState::AT_UPPER) {
        cost *= -1;
      }

      if (FieldTraits<Field>::is_strictly_positive(cost)) {
        if (state.last_cycling_iteration &&
            *state.last_cycling_iteration + 10 > state.iteration_index) {
          return i;
        }

        max_cost.record(i, cost);
      }
    }

    // logging::log_value(*max_cost.max(), "max_cost.txt");

    return max_cost.argmax();
  }

  struct NoLeaving {};
  struct ToggleBound {
    VariableState new_state;
  };
  std::variant<LeavingVariable, NoLeaving, ToggleBound>
  get_primal_leaving_variable(size_t entering,
                              const IterationState<Field>& state) {
    auto [n, d] = A_.shape();

    Matrix<Field> column(n, 1, 0);
    for (auto [row, value] : A_.get_column(entering)) {
      column[row, 0] = value;
    }

    column = state.lupa.solve_linear(std::move(column));

    // theta is maximum basic variable change so that the point would not become
    // primal infeasible
    auto get_variable_theta = [&](size_t index,
                                  Field epsilon = 0) -> std::optional<Field> {
      if (!FieldTraits<Field>::is_nonzero(column[index, 0])) {
        return std::nullopt;
      }

      Field beta = -column[index, 0];
      if (state.variables_states[entering] == VariableState::AT_UPPER) {
        beta *= -1;
      }

      Field alpha = state.basic_point[index, 0];

      auto bound = (*state.bounds)[state.basic_variables[index]];
      Field current_theta;

      if (beta > 0 && bound.upper) {
        if (alpha > *bound.upper) {
          return epsilon / beta;
        }

        return (*bound.upper - alpha + epsilon) / beta;
      }
      if (beta < 0 && bound.lower) {
        if (alpha < *bound.lower) {
          return -epsilon / beta;
        }

        return (*bound.lower - alpha - epsilon) / beta;
      }

      return std::nullopt;
    };

    // Harris’ ratio test
    // There are 2 steps:
    // 1. Determine theta_max
    // 2. Filter variables using theta_max and determine theta_chosen
    Minimum<Field, std::less<>> theta_max;

    for (size_t i = 0; i < n; ++i) {
      theta_max.record(get_variable_theta(i, 1e-10));
    }

    std::optional<size_t> leaving_id = std::nullopt;

    auto entering_state = state.variables_states[entering];
    auto entering_bound = (*state.bounds)[entering];

    if (theta_max.min()) {
      // logging::log_value(*theta_max.min(), "theta_max.txt");

      ArgMaximum<Field, std::less<>> max_pivot;

      for (size_t i = 0; i < n; ++i) {
        auto current_theta = get_variable_theta(i);

        if (current_theta && *current_theta <= *theta_max.min()) {
          max_pivot.record(i, FieldTraits<Field>::abs(column[i, 0]));
        }
      }

      assert(max_pivot.argmax().has_value());
      // logging::log_value(*max_pivot.max(), "max_pivot.txt");

      Field leaving_theta = *get_variable_theta(*max_pivot.argmax());

      // logging::log_value(leaving_theta, "leaving_theta.txt");
      // logging::log_value(state.basic_point[*max_pivot.argmax(), 0],
                         // "leaving_value.txt");

      Field new_entering_value = entering_state == VariableState::AT_LOWER
                                     ? *entering_bound.lower + leaving_theta
                                     : *entering_bound.upper - leaving_theta;

      if (entering_bound.is_inside(new_entering_value)) {
        leaving_id = *max_pivot.argmax();
      }
    }

    if (!leaving_id) {
      if (entering_state == VariableState::AT_LOWER && entering_bound.upper) {
        return ToggleBound{VariableState::AT_UPPER};
      }
      if (entering_state == VariableState::AT_UPPER && entering_bound.lower) {
        return ToggleBound{VariableState::AT_LOWER};
      }

      return NoLeaving{};
    }

    Field coef = -column[*leaving_id, 0];
    if (state.variables_states[entering] == VariableState::AT_UPPER) {
      coef *= -1;
    }

    auto leaving_state =
        coef > 0 ? VariableState::AT_UPPER : VariableState::AT_LOWER;

    return LeavingVariable{*leaving_id, leaving_state};
  }

  SimplexResult<Field> primal_implementation(
      const Bounds<Field>& bounds, const std::vector<VariableState>& states) {
    auto [n, d] = A_.shape();

    initialize_state(bounds, states);

    while (true) {
      auto rhs = get_rhs(bounds, state_.variables_states);
      state_.basic_point = state_.lupa.solve_linear(rhs);
      state_.objective = get_objective(state_);

      for (size_t i = 0; i < n; ++i) {
        auto value = state_.basic_point[i, 0];
        auto bound = bounds[state_.basic_variables[i]];

        if (!bound.is_inside(value, tolerances_.feasibility)) {
          throw std::runtime_error(
              std::format("During simplex run point became primal infeasible. "
                          "Value: {}, bounds: {}",
                          value, bound));
        }
      }

      if (state_.cycling.record(state_.iteration_index, state_.variables_states,
                                state_.objective) ==
          CyclingState::HAS_CYCLING) {
        state_.last_cycling_iteration = state_.iteration_index;
      }

      accountant_.iteration(state_);

      if (should_stop(state_)) {
        return construct_result<ReachedIterationsLimit<Field>>(state_);
      }

      auto entering = get_primal_entering_variable(state_);
      if (!entering) {
        return construct_result<FiniteLPSolution<Field>>(state_);
      }

      auto leaving_state = get_primal_leaving_variable(*entering, state_);
      if (std::holds_alternative<NoLeaving>(leaving_state)) {
        return construct_result<Unbounded>(state_);
      }

      if (std::holds_alternative<ToggleBound>(leaving_state)) {
        // logging::log(std::format("toggle bound: {}, objective = {}\n",
                                 // *entering, state_.objective),
                     // "simplex_log.txt");

        state_.variables_states[*entering] =
            std::get<ToggleBound>(leaving_state).new_state;
      } else {
        auto leaving = std::get<LeavingVariable>(leaving_state);

        // logging::log(std::format("pivot: {} -> {}, objective = {}\n",
                                 // state_.basic_variables[leaving.index],
                                 // *entering, state_.objective),
                     // "simplex_log.txt");

        state_.lupa.change_column(leaving.index, *entering);

        state_.variables_states[*entering] = VariableState::BASIC;
        state_.variables_states[state_.basic_variables[leaving.index]] =
            leaving.new_state;
        state_.basic_variables[leaving.index] = *entering;
      }

      ++state_.iteration_index;
    }
  }

 public:
  Simplex(CSCMatrix<Field> A, Matrix<Field> b, Matrix<Field> c,
          Settings<Field> settings = {})
      : A_(std::move(A)),
        b_(std::move(b)),
        c_(std::move(c)),
        state_(A_),
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

  // Algorithm for finding initial dual feasible point. It is fast, but may
  // fail. It is guaranteed to work when all variables have both upper and lower
  // bounds.
  std::optional<std::vector<VariableState>> try_get_dual_feasible(
      const Bounds<Field>& bounds) {
    auto [n, d] = A_.shape();

    auto basic_variables =
        linalg::get_row_basis(linalg::transposed(linalg::to_dense(A_)));

    return try_get_dual_feasible(bounds, basic_variables);
  }

  // Algorithm for finding initial dual feasible point. It is fast, but may
  // fail. It is guaranteed to work when all variables have both upper and lower
  // bounds.
  std::optional<std::vector<VariableState>> try_get_dual_feasible(
      const Bounds<Field>& bounds, const std::vector<size_t>& basic_variables) {
    auto [n, d] = A_.shape();

    if (basic_variables.size() != n) {
      throw std::invalid_argument(std::format(
          "Wrong basic variables count: {} != {}", basic_variables.size(), n));
    }

    std::vector<VariableState> states(d);

    state_.lupa.set_columns(basic_variables);

    auto reduced_costs = get_reduced_costs(get_basic_costs(basic_variables));

    for (size_t i = 0; i < d; ++i) {
      if (reduced_costs[i, 0] <= Field(0)) {
        if (!bounds[i].lower) {
          return std::nullopt;
        }

        states[i] = VariableState::AT_LOWER;
      } else {
        if (!bounds[i].upper) {
          return std::nullopt;
        }

        states[i] = VariableState::AT_UPPER;
      }
    }

    for (size_t basic_var : basic_variables) {
      states[basic_var] = VariableState::BASIC;
    }

    return states;
  }

  // Algorithm for finding initial dual feasible point. It runs simplex method
  // for a modified problem. This may take some time. This function returns
  // std::nullopt if and only if the original problem is infeasible.
  // Algorithm is described here:
  // https://www.cs.princeton.edu/courses/archive/fall18/cos521/Lectures/lec17.pdf
  // Important: It is required that either lower or upper bound exists
  std::optional<std::vector<VariableState>> try_get_primal_feasible(
      const Bounds<Field>& bounds) {
    auto [n, d] = A_.shape();

    auto row_basis =
        linalg::get_row_basis(linalg::transposed(linalg::to_dense(A_)));

    // add additional variable t ( + t (b - Al))
    std::vector<Field> l(n, 0);
    for (size_t i = 0; i < d; ++i) {
      if (!bounds[i].lower && !bounds[i].upper) {
        throw std::runtime_error(
            "There is a variable without both lower and upper bound. Current "
            "method does not work with it.");
      }

      Field coef = bounds[i].lower ? *bounds[i].lower : *bounds[i].upper;

      for (auto [row, value] : A_.get_column(i)) {
        l[row] += value * coef;
      }
    }

    for (size_t i = 0; i < n; ++i) {
      l[i] = b_[i, 0] - l[i];
    }

    A_.add_column(l);

    // modify cost vector
    auto old_c = std::move(c_);
    c_ = Matrix<Field>(1, d + 1);
    c_[0, d] = -1;
    // c_ = old_c.get_extended(1, d + 1, -1);

    // add bound for t variable
    Bounds<Field> extended_bounds(d + 1);
    for (size_t i = 0; i < d; ++i) {
      extended_bounds[i] = bounds[i];
    }

    extended_bounds[d] = Bound<Field>{0, 1};

    // initial variables states
    std::vector<VariableState> extended_states(d + 1);
    for (size_t i = 0; i < d; ++i) {
      if (bounds[i].lower) {
        extended_states[i] = VariableState::AT_LOWER;
      } else {
        extended_states[i] = VariableState::AT_UPPER;
      }
    }

    for (size_t i = 0; i < n; ++i) {
      extended_states[row_basis[i]] = VariableState::BASIC;
    }

    extended_states[d] = VariableState::AT_UPPER;

    // solve LP problem
    auto result = primal_implementation(extended_bounds, extended_states);

    // return everything as it was
    A_.pop_back_column();
    c_ = std::move(old_c);

    if (std::holds_alternative<FiniteLPSolution<Field>>(result.solution)) {
      auto solution = std::get<FiniteLPSolution<Field>>(result.solution);

      if (FieldTraits<Field>::is_nonzero(solution.point[d, 0])) {
        return std::nullopt;
      }

      solution.variables.pop_back();
      return std::move(solution.variables);
    }

    return std::nullopt;
  }

  // This function does not check whether matrix formed by basic columns is
  // invertible.
  bool is_dual_feasible(const Bounds<Field>& bounds,
                        const std::vector<VariableState>& states) {
    auto [n, d] = A_.shape();

    std::vector<size_t> basic_variables;
    for (size_t i = 0; i < states.size(); ++i) {
      if (states[i] == VariableState::BASIC) {
        basic_variables.push_back(i);
      }
    }

    if (basic_variables.size() != n) {
      return false;
    }

    state_.lupa.set_columns(basic_variables);
    auto reduced_costs = get_reduced_costs(get_basic_costs(basic_variables));

    for (size_t i = 0; i < states.size(); ++i) {
      if ((states[i] == VariableState::AT_LOWER &&
           (!bounds[i].lower ||
            FieldTraits<Field>::is_strictly_positive(reduced_costs[i, 0]))) ||
          (states[i] == VariableState::AT_UPPER &&
           (!bounds[i].upper ||
            FieldTraits<Field>::is_strictly_negative(reduced_costs[i, 0])))) {
        return false;
      }
    }

    return true;
  }

  // This function does not check whether matrix formed by basic columns is
  // invertible.
  bool is_primal_feasible(const Bounds<Field>& bounds,
                          const std::vector<VariableState>& states) {
    auto [n, d] = A_.shape();

    std::vector<size_t> basic_variables;
    for (size_t i = 0; i < states.size(); ++i) {
      if (states[i] == VariableState::BASIC) {
        basic_variables.push_back(i);
      }
    }

    if (basic_variables.size() != n) {
      return false;
    }

    state_.lupa.set_columns(basic_variables);
    Matrix<Field> b(b_);
    for (size_t col = 0; col < d; ++col) {
      if (states[col] == VariableState::AT_LOWER) {
        if (!bounds[col].lower) {
          return false;
        }

        for (const auto& [row, value] : A_.get_column(col)) {
          b[row, 0] -= value * *bounds[col].lower;
        }
      } else if (states[col] == VariableState::AT_UPPER) {
        if (!bounds[col].upper) {
          return false;
        }

        for (const auto& [row, value] : A_.get_column(col)) {
          b[row, 0] -= value * *bounds[col].upper;
        }
      }
    }

    auto point = state_.lupa.solve_linear(b);
    for (size_t i = 0; i < n; ++i) {
      if (!bounds[basic_variables[i]].is_inside(point[i, 0])) {
        return false;
      }
    }

    return true;
  }

  // Point associated with the given states must be dual feasible
  SimplexResult<Field> dual(const Bounds<Field>& bounds,
                            const std::vector<VariableState>& states) {
    try {
      return dual_implementation(bounds, states);
    } catch (...) {
      SimplexCoreDump<Field>(A_, b_, c_).dump_state(state_);
      throw;
    }
  }

  // Point associated with the given states must be primal feasible
  SimplexResult<Field> primal(const Bounds<Field>& bounds,
                              const std::vector<VariableState>& states) {
    try {
      return primal_implementation(bounds, states);
    } catch (...) {
      SimplexCoreDump<Field>(A_, b_, c_).dump_state(state_);
      throw;
    }
  }
};

}  // namespace simplex
