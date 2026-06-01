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
#include "Feasibility.h"
#include "Math.h"
#include "SimplexCoreDump.h"
#include "Tolerance.h"
#include "linear/matrix/Matrix.h"
#include "linear/matrix/NPY.h"
#include "linear/matrix/Norms.h"
#include "linear/matrix/RowBasis.h"
#include "linear/model/LP.h"
#include "linear/sparse/LU.h"
#include "pricing/primal/MostInfeasible.h"
#include "ratio/primal/Harris.h"
#include "utils/Accumulators.h"
#include "utils/Variant.h"

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
template <typename Field, typename Accountant = EmptyAccountant<Field>>
class Simplex {
  const CSCMatrix<Field> A_;
  const Matrix<Field> b_;
  const Matrix<Field> c_;

  IterationState<Field> state_;

  Config<Field> config_;
  Accountant accountant_;

  template <typename F>
    requires std::invocable<F> &&
             std::same_as<std::invoke_result_t<F>, std::optional<std::string>>
  void validate(F validator) const {
    if (config_.validate_input) {
      if (auto result = validator()) {
        throw std::runtime_error(*result);
      }
    }
  }

  template <typename T>
  static SimplexResult<Field> construct_result(
      const IterationState<Field>& state) {
    auto point =
        detail::get_point_from_basis(*state.bounds, state.variables_states,
                                     state.basic_variables, state.basic_point);

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

  std::optional<size_t> get_dual_entering_variable(
      LeavingVariable leaving, const IterationState<Field>& state) const {
    auto [n, d] = state.problem_shape();

    ArgMinimum<Field> min_ratio;

    const auto simplex_multipliers = state_.lupa.solve_linear_transposed(
        detail::get_basic_cost(c_, state_.basic_variables));

    const auto reduced_costs =
        detail::get_reduced_cost(A_, c_, simplex_multipliers);

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

    return min_ratio.has_value() ? std::optional{min_ratio->index}
                                 : std::nullopt;
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

    if (config_.validate_input) {
      if (!is_dual_feasible(A_, b_, c_, bounds, states)) {
        throw std::invalid_argument(
            "Given initial state is not dual feasible.");
      }
    }

    // initialize simplex state
    initialize_state(bounds, states);

    DualLeavingVariable<Field> leaving_finder;

    while (true) {
      auto rhs =
          detail::get_adjusted_rhs(A_, b_, bounds, state_.variables_states);
      state_.basic_point = state_.lupa.solve_linear(rhs);

      state_.objective =
          detail::get_objective(c_, *state_.bounds, state_.variables_states,
                                state_.basic_variables, state_.basic_point);

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

  bool lost_primal_feasibility(const Bounds<Field>& bounds) const {
    auto [n, d] = A_.shape();

    for (size_t i = 0; i < n; ++i) {
      auto value = state_.basic_point[i, 0];
      auto bound = bounds[state_.basic_variables[i]];

      if (!bound.contains(value, config_.tolerance.feasibility)) {
        std::println("infeasible: variable {} with value {} not in {}",
                     state_.basic_variables[i], value, bound);
        return true;
      }
    }

    return false;
  }

  // It is guaranteed, that after execution of this method, if finite LP
  // solution was found, then inside LUPA basic variables would be selected as
  // columns.
  SimplexResult<Field> primal_implementation(
      const Bounds<Field>& bounds, const std::vector<VariableState>& states) {
    auto [n, d] = A_.shape();

    std::vector<Bound<Field>> bounds_vector(d);
    for (size_t i = 0; i < d; ++i) {
      bounds_vector[i] = bounds[i];
    }

    initialize_state(bounds, states);

    if (!config_.primal_pricing) {
      throw std::runtime_error(
          "Primal pricing must be specified in simplex config.");
    }

    while (true) {
      auto rhs =
          detail::get_adjusted_rhs(A_, b_, bounds, state_.variables_states);
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
        std::println("rollback due to infeasibility");

        state_.lupa.refactorize();
        continue;
      }

      state_.objective =
          detail::get_objective(c_, *state_.bounds, state_.variables_states,
                                state_.basic_variables, state_.basic_point);

      if (state_.cycling.record(state_.iteration_index, state_.variables_states,
                                state_.objective) ==
          CyclingState::HAS_CYCLING) {
        state_.last_cycling_iteration = state_.iteration_index;
      }

      accountant_.iteration(state_);

      if (should_stop(state_)) {
        return construct_result<ReachedIterationsLimit<Field>>(state_);
      }

      const auto simplex_multipliers = state_.lupa.solve_linear_transposed(
          detail::get_basic_cost(c_, state_.basic_variables));

      const auto reduced_costs_matrix =
          detail::get_reduced_cost(A_, c_, simplex_multipliers);

      // temporary: transform matrix to vector
      std::vector<Field> reduced_costs(d);
      for (size_t i = 0; i < d; ++i) {
        reduced_costs[i] = reduced_costs_matrix[i, 0];
      }

      std::vector<Field> basic_point(n);
      for (size_t i = 0; i < n; ++i) {
        basic_point[i] = state_.basic_point[i, 0];
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

      auto action =
          detail::primal_ratio_test(A_, *entering, reduced_costs[*entering],
                                    state_, config_.tolerance.pivot);

      if (std::holds_alternative<detail::Unbounded>(action)) {
        return construct_result<Unbounded>(state_);
      }

      std::visit(Overload{
                     [this](detail::ChangeBasis action) {
                       change_basis(action.leaving_index,
                                    action.entering_variable, action.new_state);
                     },
                     [this](detail::ToggleBound action) {
                       change_bound(action.variable, action.new_state);
                     },
                     [](auto) { throw std::runtime_error("Invalid action."); },
                 },
                 action);

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

  void change_basis(size_t leaving_index, size_t entering_variable,
                    VariableState leaving_state) {
    state_.lupa.change_column(leaving_index, entering_variable);

    state_.variables_states[entering_variable] = VariableState::BASIC;

    state_.variables_states[state_.basic_variables[leaving_index]] =
        leaving_state;
    state_.basic_variables[leaving_index] = entering_variable;
  }

  void change_bound(size_t variable_index, VariableState new_bound) {
    validate([&] -> std::optional<std::string> {
      if (new_bound != VariableState::AT_LOWER ||
          new_bound != VariableState::AT_UPPER) {
        return "new_bound must be either AT_LOWER or AT_UPPER.";
      }

      if (new_bound == VariableState::AT_LOWER &&
          !(*state_.bounds)[variable_index].lower) {
        return "new_bound is set to AT_LOWER, but variable doesn't have a "
               "lower bound.";
      }

      if (new_bound == VariableState::AT_UPPER &&
          !(*state_.bounds)[variable_index].upper) {
        return "new_bound is set to AT_UPPER, but variable doesn't have an "
               "upper bound.";
      }

      return std::nullopt;
    });

    state_.variables_states[variable_index] = new_bound;
  }

  // current basis getters
  std::vector<size_t> get_basic_vars() const { return state_.basic_variables; }

  std::vector<VariableState> get_states() const {
    return state_.variables_states;
  }

  Matrix<Field> get_tableau_row(size_t row) const {
    return state_.lupa.get_row(row);
  }
};

}  // namespace simplex
