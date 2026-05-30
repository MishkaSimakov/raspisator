#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <random>
#include <ranges>
#include <unordered_map>
#include <variant>

#include "Accountant.h"
#include "Config.h"
#include "CyclingDetector.h"
#include "Dual.h"
#include "SimplexCoreDump.h"
#include "Tolerance.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/NPY.h"
#include "linear/matrix/Norms.h"
#include "linear/matrix/RowBasis.h"
#include "linear/model/LP.h"
#include "linear/sparse/LU.h"
#include "pricing/primal/MostInfeasible.h"
#include "utils/Accumulators.h"
#include "utils/Variant.h"

namespace simplex {

struct NoLeaving {};
struct NoEntering {};

struct ToggleBound {
  size_t variable_index;
  VariableState new_state;  // should be either AT_UPPER or AT_LOWER
};

struct ChangeBasicVariable {
  size_t leaving_index;      // index of leaving variable in basic_variables
  size_t leaving_variable;   // index of leaving variable
  size_t entering_variable;  // index of entering variable

  VariableState entering_old_state;  // should be either AT_UPPER or AT_LOWER
  VariableState leaving_new_state;   // should be either AT_UPPER or AT_LOWER
};

using IterationAction =
    std::variant<NoLeaving, NoEntering, ToggleBound, ChangeBasicVariable>;

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
  const CSCMatrix<Field> A_;
  const Matrix<Field> b_;
  const Matrix<Field> c_;

  IterationState<Field> state_;

  Config<Field> config_;
  Accountant accountant_;

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
      } else if (state.variables_states[i] == VariableState::NONBASIC_FREE) {
        result[i, 0] = 0;
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

    const auto reduced_costs =
        get_reduced_costs(get_basic_costs(state.basic_variables));

    const auto inverse_row = state.lupa.get_row(leaving.index);

    for (size_t i = 0; i < d; ++i) {
      if (state.variables_states[i] == VariableState::BASIC) {
        continue;
      }

      Field coef = 0;
      for (const auto& [row, value] : A_.get_column(i)) {
        coef += inverse_row[row, 0] * value;
      }

      if (!FieldTraits<Field>::is_nonzero(coef)) {
        continue;
      }

      // TODO: think about drop tolerance
      const Field cost = FieldTraits<Field>::is_nonzero(reduced_costs[i, 0])
                             ? reduced_costs[i, 0]
                             : 0;

      if (!((state.variables_states[i] == VariableState::AT_LOWER &&
             cost < Field(1) / Field(1e5)) ||
            (state.variables_states[i] == VariableState::AT_UPPER &&
             cost > -Field(1) / Field(1e5)))) {
        throw std::runtime_error(
            std::format("Current point is not dual feasible! Reduced cost for "
                        "variable #{} has value {}.",
                        i, reduced_costs[i, 0]));
      }

      Field ratio = cost / coef;

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

    return min_ratio->index;
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
      // TODO: change this to switch, so that new VariableStates can be handled
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
    return config_.max_iterations &&
           state.iteration_index >= config_.max_iterations;
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

    if (config_.is_strict) {
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

  struct PrimalEnteringVariable {
    size_t index;
    Field reduced_cost;
  };

  IterationAction get_primal_leaving_variable(
      PrimalEnteringVariable entering, const IterationState<Field>& state) {
    using std::abs;

    const auto [n, d] = A_.shape();

    Matrix<Field> column(n, 1, 0);
    for (auto [row, value] : A_.get_column(entering.index)) {
      column[row, 0] = value;
    }

    column = state.lupa.solve_linear(std::move(column));

    // theta is maximum entering variable change so that the point would not
    // become primal infeasible for variable i,
    // change = | new_value - old_value |
    const auto get_variable_theta =
        [&](const size_t i, const Field epsilon = 0) -> std::optional<Field> {
      if (abs(column[i, 0]) < config_.tolerance.pivot) {
        return std::nullopt;
      }

      // x_i = \alpha + s \beta x_j, where
      // x_j is entering variable,
      // s = -sign(entering reduced cost)
      const Field alpha = state.basic_point[i, 0];
      const Field beta =
          entering.reduced_cost > 0 ? -column[i, 0] : column[i, 0];

      const auto bound = (*state.bounds)[state.basic_variables[i]];

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

    // Harris' ratio test
    // There are 2 steps:
    // 1. Determine theta_max
    // 2. Filter variables using theta_max and determine theta_chosen
    Minimum<Field> min_theta_bound;

    for (size_t i = 0; i < n; ++i) {
      min_theta_bound.record(
          get_variable_theta(i, FieldTraits<Field>::tolerance));
    }

    std::optional<size_t> leaving_id = std::nullopt;

    const auto entering_state = state.variables_states[entering.index];
    const auto entering_bound = (*state.bounds)[entering.index];

    if (min_theta_bound.has_value()) {
      const Field theta_max = *min_theta_bound;
      // logging::log_value(*theta_max.min(), "theta_max.txt");

      ArgMaximum<Field> max_pivot;

      for (size_t i = 0; i < n; ++i) {
        auto current_theta = get_variable_theta(i);

        if (current_theta && *current_theta <= theta_max) {
          max_pivot.record(i, abs(column[i, 0]));
        }
      }

      // logging::log_value(*max_pivot.max(), "max_pivot.txt");

      const Field leaving_theta = *get_variable_theta(max_pivot->index);

      // logging::log_value(leaving_theta, "leaving_theta.txt");
      // logging::log_value(state.basic_point[*max_pivot.argmax(), 0],
      // "leaving_value.txt");

      Field new_entering_value;

      switch (entering_state) {
        case VariableState::AT_LOWER:
          new_entering_value = *entering_bound.lower + leaving_theta;
          break;
        case VariableState::AT_UPPER:
          new_entering_value = *entering_bound.upper - leaving_theta;
          break;
        case VariableState::NONBASIC_FREE:
          new_entering_value =
              entering.reduced_cost > 0 ? leaving_theta : -leaving_theta;
          break;
        default:
          throw std::runtime_error("Unexpected variable state.");
      }

      if (entering_bound.is_inside(new_entering_value)) {
        leaving_id = max_pivot->index;
      }
    }

    if (!leaving_id) {
      if (entering_state == VariableState::AT_LOWER && entering_bound.upper) {
        return ToggleBound{entering.index, VariableState::AT_UPPER};
      }
      if (entering_state == VariableState::AT_UPPER && entering_bound.lower) {
        return ToggleBound{entering.index, VariableState::AT_LOWER};
      }

      return NoLeaving{};
    }

    const auto leaving_state =
        entering.reduced_cost * column[*leaving_id, 0] < 0
            ? VariableState::AT_UPPER
            : VariableState::AT_LOWER;

    return ChangeBasicVariable{
        .leaving_index = *leaving_id,
        .leaving_variable = state_.basic_variables[*leaving_id],
        .entering_variable = entering.index,
        .entering_old_state = state.variables_states[entering.index],
        .leaving_new_state = leaving_state,
    };
  }

  bool lost_primal_feasibility(const Bounds<Field>& bounds) const {
    auto [n, d] = A_.shape();

    for (size_t i = 0; i < n; ++i) {
      auto value = state_.basic_point[i, 0];
      auto bound = bounds[state_.basic_variables[i]];

      if (!bound.is_inside(value, config_.tolerance.feasibility)) {
        std::println("infeasible: variable {} with value {} not in {}",
                     state_.basic_variables[i], value, bound);
        return true;
      }
    }

    return false;
  }

  IterationAction get_rollback(IterationAction action) {
    return std::visit(
        Overload{
            [](ToggleBound action) -> IterationAction {
              return ToggleBound{
                  .variable_index = action.variable_index,
                  .new_state = action.new_state == VariableState::AT_LOWER
                                   ? VariableState::AT_UPPER
                                   : VariableState::AT_LOWER,
              };
            },
            [](ChangeBasicVariable action) -> IterationAction {
              return ChangeBasicVariable{
                  .leaving_index = action.leaving_index,
                  .leaving_variable = action.entering_variable,
                  .entering_variable = action.leaving_variable,
                  .entering_old_state = action.leaving_new_state,
                  .leaving_new_state = action.entering_old_state,
              };
            },
            [](auto /* action */) -> IterationAction { std::unreachable(); }},
        action);
  }

  void apply_action(IterationAction action) {
    return std::visit(
        Overload{[this](ToggleBound action) {
                   state_.variables_states[action.variable_index] =
                       action.new_state;

                   // std::println("  toggle bound: {}", action.variable_index);
                 },
                 [this](ChangeBasicVariable action) {
                   state_.lupa.change_column(action.leaving_index,
                                             action.entering_variable);

                   state_.variables_states[action.entering_variable] =
                       VariableState::BASIC;
                   state_.variables_states
                       [state_.basic_variables[action.leaving_index]] =
                       action.leaving_new_state;
                   state_.basic_variables[action.leaving_index] =
                       action.entering_variable;

                   // std::println("  basis change: {} -> {}",
                   // action.leaving_variable, action.entering_variable);
                 },
                 [](auto /* action */) { std::unreachable(); }},
        action);
  }

  // It is guaranteed, that after execution of this method, if finite LP
  // solution was found, then inside LUPA basic variables would be selected as
  // columns.
  SimplexResult<Field> primal_implementation(
      const Bounds<Field>& bounds, const std::vector<VariableState>& states) {
    auto [n, d] = A_.shape();

    std::vector<IterationAction> history;

    initialize_state(bounds, states);

    if (!config_.primal_pricing) {
      throw std::runtime_error(
          "Primal pricing must be specified in simplex config.");
    }

    while (true) {
      auto rhs = get_rhs(bounds, state_.variables_states);
      state_.basic_point = state_.lupa.solve_linear(rhs);

      // if constexpr (std::same_as<Field, double>) {
      //   const auto dense_submatrix =
      //       linalg::to_dense(A_).get_columns(state_.basic_variables);
      //
      //   Matrix<double> x(n, 1, 1);
      //   const auto Ax = dense_submatrix * x;
      //
      //   const auto x_lupa = state_.lupa.solve_linear(Ax);
      //
      //   const auto residue = x - x_lupa;
      //
      //   std::println("|| residue || = {}", linalg::inf_norm(residue));
      //   logging::log_value(linalg::inf_norm(residue), "residue_norm.csv");
      // }

      if (lost_primal_feasibility(bounds)) {
        assert(!history.empty());

        std::println("rollback due to infeasibility");

        state_.lupa.refactorize();

        const auto rollback = get_rollback(history.back());

        history.pop_back();
        apply_action(rollback);

        continue;
      }

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

      const auto reduced_costs_matrix =
          get_reduced_costs(get_basic_costs(state_.basic_variables));

      // temporary: transform matrix to vector
      std::vector<Field> reduced_costs(d);
      for (size_t i = 0; i < d; ++i) {
        reduced_costs[i] = reduced_costs_matrix[i, 0];
      }

      std::vector<Field> basic_point(n);
      for (size_t i = 0; i < n; ++i) {
        basic_point[i] = state_.basic_point[i, 0];
      }

      std::vector<Bound<Field>> bounds_vector(d);
      for (size_t i = 0; i < d; ++i) {
        bounds_vector[i] = bounds[i];
      }

      auto entering =
          config_.primal_pricing->get_primal_entering(detail::State<Field>{
              .iteration = state_.iteration_index,
              .objective = state_.objective,
              .basic_point = basic_point,
              .bounds = bounds_vector,
              .reduced_cost = reduced_costs,
              .states = state_.variables_states,
              .basic_vars = state_.basic_variables,
              .tolerance = config_.tolerance,
          });

      if (!entering) {
        return construct_result<FiniteLPSolution<Field>>(state_);
      }

      auto action = get_primal_leaving_variable(
          PrimalEnteringVariable{
              .index = *entering,
              .reduced_cost = reduced_costs[*entering],
          },
          state_);
      if (std::holds_alternative<NoLeaving>(action)) {
        return construct_result<Unbounded>(state_);
      }

      history.push_back(action);
      apply_action(action);

      ++state_.iteration_index;
    }
  }

 public:
  Simplex(CSCMatrix<Field> A, Matrix<Field> b, Matrix<Field> c,
          Config<Field> settings = {})
      : A_(std::move(A)),
        b_(std::move(b)),
        c_(std::move(c)),
        state_(A_),
        config_(std::move(settings)) {
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
    config_.max_iterations = max_iterations;
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

  // Returns primal feasible basis or std::nullopt if problem is infeasible.
  // Note: algorithm is taken from
  // https://people.orie.cornell.edu/dpw/orie6300/Lectures/lec12.pdf
  std::optional<std::vector<VariableState>> get_primal_feasible(
      Bounds<Field> bounds) {
    const auto [n, old_d] = A_.shape();
    const size_t new_d = old_d + n;

    std::vector<VariableState> states(new_d);
    Bounds<Field> new_bounds(new_d);

    for (size_t i = 0; i < old_d; ++i) {
      if (bounds[i].lower) {
        states[i] = VariableState::AT_LOWER;
      } else if (bounds[i].upper) {
        states[i] = VariableState::AT_UPPER;
      } else {
        states[i] = VariableState::NONBASIC_FREE;
      }

      new_bounds[i] = bounds[i];
    }

    for (size_t i = old_d; i < new_d; ++i) {
      states[i] = VariableState::BASIC;
      new_bounds[i] = {0, std::nullopt};
    }

    // add additional slack variable to each constraint
    const auto rhs = get_rhs(bounds, states);

    auto new_A = A_;
    for (size_t i = 0; i < n; ++i) {
      new_A.add_column();

      if (rhs[i, 0] > 0) {
        new_A.push_to_last_column(i, 1);
      } else {
        new_A.push_to_last_column(i, -1);
      }
    }

    auto new_c = Matrix<Field>(1, new_d, 0);

    for (size_t i = old_d; i < new_d; ++i) {
      new_c[0, i] = -1;
    }

    auto helper = Simplex(
        new_A, b_, new_c,
        {
            .primal_pricing = std::make_unique<PrimalMostInfeasible<Field>>(),
        });
    const auto result = helper.primal(new_bounds, states);

    if (!result.is_feasible()) {
      throw std::runtime_error(
          "Something went wrong in primal implementation.");
    }

    auto solution = std::get<FiniteLPSolution<Field>>(result.solution);

    // Case 1
    if (FieldTraits<Field>::is_strictly_negative(solution.value)) {
      // Problem is infeasible
      // TODO: indicate infeasibility in some way
      return std::nullopt;
    }

    // Case 2: try to eliminate artificial variables from basic variables (if
    // there are any) using pivot operation
    for (size_t basic_index = 0; basic_index < n; ++basic_index) {
      const size_t i = helper.state_.basic_variables[basic_index];

      if (i < old_d) {
        continue;
      }

      const auto row = helper.state_.lupa.get_row(basic_index);

      // try to find replacement for i among non-artificial variables
      bool found = false;

      for (size_t j = 0; j < old_d; ++j) {
        if (solution.variables[j] == VariableState::BASIC) {
          continue;
        }

        Field coef = 0;

        for (const auto [index, value] : A_.get_column(j)) {
          coef += row[index, 0] * value;
        }

        if (FieldTraits<Field>::is_nonzero(coef)) {
          // change i -> j in basis
          solution.variables[j] = VariableState::BASIC;
          helper.state_.lupa.change_column(basic_index, j);
          found = true;

          break;
        }
      }

      if (!found) {
        // Problem contains linearly dependent rows.
        // TODO: indicate this to user
        return std::nullopt;
      }
    }

    solution.variables.resize(old_d);
    return solution.variables;
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
