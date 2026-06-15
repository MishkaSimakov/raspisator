#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <optional>
#include <random>
#include <ranges>
#include <unordered_map>
#include <variant>

#include "Accountant.h"
#include "Config.h"
#include "CyclingDetector.h"
#include "Feasibility.h"
#include "SimplexCoreDump.h"
#include "SimplexMath.h"
#include "Tolerance.h"
#include "linalg/Matrix.h"
#include "linalg/NPY.h"
#include "linalg/Norm.h"
#include "linalg/RowBasis.h"
#include "linalg/lu/LUPA.h"
#include "problem/StandardLP.h"
#include "ratio/primal/Harris.h"
#include "simplex/Result.h"
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
template <typename Field>
class Simplex {
  const problem::StandardLP<Field>* problem_ = nullptr;

  // current matrix factorization
  std::optional<linalg::LUPA<Field>> lupa_;

  // current simplex position in problem space
  std::vector<size_t> basic_vars_;
  std::vector<VariableState> var_states_;
  Vector<Field> basic_point_;

  Vector<Field> reduced_cost_;
  Field objective_;

  size_t iteration_;
  bool intentional_repeat_;

  Config<Field> config_;

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
  SimplexResult<Field> construct_result() const {
    auto point = get_point();
    Field objective = linalg::dot(problem_->cost, point);

    if constexpr (std::same_as<T, FiniteLPSolution<Field>>) {
      return SimplexResult<Field>{
          .iterations_count = iteration_,
          .solution =
              FiniteLPSolution{std::move(point), objective, var_states_},
      };
    }
    if constexpr (std::same_as<T, NoFeasibleElements>) {
      return SimplexResult<Field>{
          .iterations_count = iteration_,
          .solution = NoFeasibleElements{},
      };
    }
    if constexpr (std::same_as<T, ReachedIterationsLimit<Field>>) {
      return SimplexResult<Field>{
          .iterations_count = iteration_,
          .solution = ReachedIterationsLimit<Field>{objective},
      };
    }
    if constexpr (std::same_as<T, Unbounded>) {
      return SimplexResult<Field>{
          .iterations_count = iteration_,
          .solution = Unbounded{},
      };
    }
  }

  std::optional<size_t> get_dual_entering_variable(
      LeavingVariable leaving, const Vector<Field>& reduced_cost,
      const Vector<Field>& leaving_row) const {
    auto [n, d] = problem_->matrix.shape();

    ArgMinimum<Field> min_ratio;

    for (size_t i = 0; i < d; ++i) {
      if (var_states_[i] == VariableState::BASIC) {
        continue;
      }

      Field coef = 0;
      for (const auto& [row, value] : problem_->matrix.get_column(i)) {
        coef += leaving_row[row] * value;
      }

      if (!FieldTraits<Field>::is_nonzero(coef)) {
        continue;
      }

      // TODO: think about drop tolerance
      const Field cost =
          FieldTraits<Field>::is_nonzero(reduced_cost[i]) ? reduced_cost[i] : 0;

      Field ratio = cost / coef;

      if (leaving.new_state == VariableState::AT_UPPER) {
        ratio *= -1;
      }

      if (leaving.new_state == VariableState::AT_LOWER) {
        if (var_states_[i] == VariableState::AT_LOWER && coef > Field(0) ||
            var_states_[i] == VariableState::AT_UPPER && coef < Field(0)) {
          continue;
        }
      } else {
        if (var_states_[i] == VariableState::AT_LOWER && coef < Field(0) ||
            var_states_[i] == VariableState::AT_UPPER && coef > Field(0)) {
          continue;
        }
      }

      min_ratio.record(i, ratio);
    }

    return min_ratio.has_value() ? std::optional{min_ratio->index}
                                 : std::nullopt;
  }

  void initialize_state(const std::vector<VariableState>& states) {
    auto [n, d] = problem_->matrix.shape();

    basic_vars_.clear();
    var_states_ = states;

    iteration_ = 0;
    intentional_repeat_ = false;

    for (size_t i = 0; i < d; ++i) {
      if (states[i] == VariableState::BASIC) {
        basic_vars_.push_back(i);
      }
    }

    lupa_->set_columns(basic_vars_);
    recalculate_incremental();
  }

  std::optional<std::string> validate_states(
      const std::vector<VariableState>& states) const {
    const auto [n, d] = problem_->matrix.shape();

    if (states.size() != d) {
      return std::format("Wrong states size: {} != {}.", states.size(), d);
    }

    size_t basic_count = 0;

    for (size_t i = 0; i < d; ++i) {
      switch (states[i]) {
        case VariableState::AT_LOWER:
          if (!problem_->var_bounds[i].lower) {
            return std::format(
                "Variable {} is set to AT_LOWER, but doesn't have a lower "
                "bound.",
                problem_->var_name(i));
          }
          break;
        case VariableState::AT_UPPER:
          if (!problem_->var_bounds[i].upper) {
            return std::format(
                "Variable {} is set to AT_UPPER, but doesn't have an upper "
                "bound.",
                problem_->var_name(i));
          }
          break;
        case VariableState::NONBASIC_FREE:
          if (!problem_->var_bounds[i].is_free()) {
            return std::format(
                "Variable {} is set to NONBASIC_FREE, but isn't free.",
                problem_->var_name(i));
          }
          break;
        case VariableState::BASIC:
          ++basic_count;
          break;
        default:
          throw std::runtime_error("Unknown variable state.");
      }
    }

    if (basic_count != n) {
      return std::format("Wrong basic variables count: {} != {}.", basic_count,
                         n);
    }

    return std::nullopt;
  }

  bool violate_primal_bounds() const {
    for (size_t i = 0; i < basic_point_.size(); ++i) {
      const auto& bound = problem_->var_bounds[basic_vars_[i]];

      if (!bound.contains(basic_point_[i], config_.tolerance.feasibility)) {
        std::println("infeasible: variable {} with value {} not in {}",
                     basic_vars_[i], basic_point_[i], bound);
        return true;
      }
    }

    return false;
  }

  bool violate_dual_bounds(const Vector<Field>& reduced_cost) const {
    using std::abs;
    const auto [n, d] = problem_->matrix.shape();

    for (size_t i = 0; i < d; ++i) {
      switch (var_states_[i]) {
        case VariableState::AT_LOWER:
          if (reduced_cost[i] > config_.tolerance.feasibility) {
            return true;
          }
          break;
        case VariableState::AT_UPPER:
          if (reduced_cost[i] < -config_.tolerance.feasibility) {
            return true;
          }
          break;
        case VariableState::NONBASIC_FREE:
          if (abs(reduced_cost[i]) > config_.tolerance.feasibility) {
            return true;
          }
          break;
        case VariableState::BASIC:
          break;
        default:
          throw std::runtime_error("Unknown variable state.");
      }
    }

    return false;
  }

  enum class IterationResult {
    FEASIBLE,
    INFEASIBLE,
    UNBOUNDED,
    REFACTORIZE_REPEAT,
    MOVED,
  };

  IterationResult primal_iteration() {
    using std::abs;

    basic_point_ =
        lupa_->solve_linear(detail::get_adjusted_rhs(*problem_, var_states_));

    if (violate_primal_bounds()) {
      if (lupa_->get_changes_since_refactorization() > 0) {
        return IterationResult::REFACTORIZE_REPEAT;
      }

      throw std::runtime_error(
          "Point became primal infeasible. Not implemented.");
    }

    StateView<Field> state_view{
        .problem = *problem_,
        .lupa = *lupa_,
        .iteration = iteration_,
        .objective = objective_,
        .basic_point = basic_point_,
        .states = var_states_,
        .basic_vars = basic_vars_,
        .intentional_repeat = intentional_repeat_,
        .tolerance = config_.tolerance,
    };

    if (config_.accountant) {
      config_.accountant->iteration(state_view);
    }

    auto entering =
        config_.primal_pricing->get_primal_entering(state_view, reduced_cost_);

    if (!entering) {
      return IterationResult::FEASIBLE;
    }

    const Vector pivot_col =
        lupa_->solve_linear(problem_->matrix.get_column_as_matrix(*entering));

    auto action = detail::primal_ratio_test(
        *problem_, state_view, *entering, reduced_cost_[*entering], pivot_col);

    if (auto* move = std::get_if<ChangeBasisMove<Field>>(&action)) {
      // suspicious pivot, refactorize and try again
      if (abs(pivot_col[move->leaving_index]) <
              config_.tolerance.suspicious_pivot &&
          lupa_->get_changes_since_refactorization() > 0) {
        std::cout << iteration_ << " suspicious pivot" << std::endl;

        return IterationResult::REFACTORIZE_REPEAT;
      }

      const Vector pivot_row = linalg::transpose(problem_->matrix) *
                               lupa_->get_row(move->leaving_index);

      // notify primal pricing
      config_.primal_pricing->move(*move, state_view, pivot_row, pivot_col);

      // incremental objective update
      objective_ += reduced_cost_[move->entering_variable] * move->step_length;

      // incremental reduced cost update
      const Field entering_reduced_cost =
          reduced_cost_[move->entering_variable];

      reduced_cost_[basic_vars_[move->leaving_index]] =
          -entering_reduced_cost / pivot_row[move->entering_variable];

      for (size_t i = 0; i < problem_->matrix.cols(); ++i) {
        if (i == basic_vars_[move->leaving_index]) {
          continue;
        }

        if (var_states_[i] == VariableState::BASIC ||
            i == move->entering_variable) {
          reduced_cost_[i] = 0;
        } else {
          reduced_cost_[i] -= entering_reduced_cost * pivot_row[i] /
                              pivot_row[move->entering_variable];
        }
      }

      change_basis(move->leaving_index, move->entering_variable,
                   move->new_state);

      return IterationResult::MOVED;
    }

    if (auto* move = std::get_if<ToggleBoundMove<Field>>(&action)) {
      config_.primal_pricing->move(*move, state_view);

      // reduced cost doesn't change
      objective_ += reduced_cost_[move->variable] * move->step_length;

      change_bound(move->variable, move->new_state);

      return IterationResult::MOVED;
    }

    return IterationResult::UNBOUNDED;
  }

  void recalculate_incremental() {
    // recalculate incremental reduced cost
    const Vector pi =
        lupa_->solve_linear_transposed(Vector(problem_->cost[basic_vars_]));

    reduced_cost_ = problem_->cost - linalg::transpose(problem_->matrix) * pi;

    // recalculate incremental objective value
    const Vector basic_point =
        lupa_->solve_linear(detail::get_adjusted_rhs(*problem_, var_states_));

    objective_ = detail::get_objective(problem_->cost, problem_->var_bounds,
                                       var_states_, basic_vars_, basic_point);
  }

  void primal_refactorize() {
    lupa_->refactorize();

    recalculate_incremental();

    config_.primal_pricing->post_refactorization(*problem_, *lupa_, var_states_,
                                                 basic_vars_);
  }

  // It is guaranteed, that after execution of this method, if finite LP
  // solution was found, then inside LUPA basic variables would be selected as
  // columns.
  SimplexResult<Field> primal_implementation(
      const std::vector<VariableState>& states) {
    using std::abs;

    if (!config_.primal_pricing) {
      throw std::logic_error(
          "Primal pricing must be specified in simplex config.");
    }

    auto [n, d] = problem_->matrix.shape();

    initialize_state(states);

    config_.primal_pricing->init(*problem_, *lupa_, var_states_, basic_vars_);

    std::cout << "starting, objective is " << objective_ << std::endl;

    //
    while (true) {
      if (config_.max_iterations && iteration_ > *config_.max_iterations) {
        return construct_result<ReachedIterationsLimit<Field>>();
      }

      const auto result = primal_iteration();

      ++iteration_;

      switch (result) {
        case IterationResult::FEASIBLE:
          return construct_result<FiniteLPSolution<Field>>();
        case IterationResult::UNBOUNDED:
          return construct_result<Unbounded>();
        case IterationResult::INFEASIBLE:
          return construct_result<NoFeasibleElements>();
        case IterationResult::MOVED:
          intentional_repeat_ = false;
          break;
        case IterationResult::REFACTORIZE_REPEAT:
          primal_refactorize();
          intentional_repeat_ = true;
          break;
        default:
          throw std::runtime_error("Unknown iteration result");
      }

      if (lupa_->get_changes_since_refactorization() > 250) {
        primal_refactorize();
      }
    }
  }

  SimplexResult<Field> dual_implementation(
      const std::vector<VariableState>& states) {
    if (!config_.dual_pricing) {
      throw std::runtime_error(
          "Dual pricing must be specified in simplex config.");
    }

    auto [n, d] = problem_->matrix.shape();

    initialize_state(states);
    size_t iteration = 0;

    while (true) {
      Vector rhs = detail::get_adjusted_rhs(*problem_, var_states_);
      basic_point_ = lupa_->solve_linear(rhs);

      Field objective =
          detail::get_objective(problem_->cost, problem_->var_bounds,
                                var_states_, basic_vars_, basic_point_);

      StateView<Field> state_view{
          .problem = *problem_,
          .lupa = *lupa_,
          .iteration = iteration,
          .objective = objective,
          .basic_point = basic_point_,
          .states = var_states_,
          .basic_vars = basic_vars_,
          .tolerance = config_.tolerance,
      };

      if (config_.accountant) {
        config_.accountant->iteration(state_view);
      }

      if (config_.max_iterations && iteration > *config_.max_iterations) {
        return construct_result<ReachedIterationsLimit<Field>>();
      }

      auto leaving = config_.dual_pricing->get_dual_leaving(state_view);
      if (!leaving) {
        return construct_result<FiniteLPSolution<Field>>();
      }

      const Vector pi =
          lupa_->solve_linear_transposed(Vector(problem_->cost[basic_vars_]));
      const Vector reduced_cost =
          problem_->cost - linalg::transpose(problem_->matrix) * pi;

      if (violate_dual_bounds(reduced_cost)) {
        throw std::runtime_error("Not implemented.");
      }

      const Vector leaving_row = lupa_->get_row(leaving->index);

      auto entering =
          get_dual_entering_variable(*leaving, reduced_cost, leaving_row);
      if (!entering) {
        return construct_result<NoFeasibleElements>();
      }

      change_basis(leaving->index, *entering, leaving->new_state);

      ++iteration;
    }
  }

 public:
  explicit Simplex(Config<Field> config = {}) : config_(std::move(config)) {}

  // config setters
  void set_problem(const problem::StandardLP<Field>& problem) {
    validate([&] -> std::optional<std::string> {
      problem.validate();
      return std::nullopt;
    });

    problem_ = &problem;
    lupa_.emplace(problem.matrix);
  }

  void set_validate_input(bool value) { config_.validate_input = value; }

  void set_max_iterations(std::optional<size_t> max_iterations) {
    config_.max_iterations = max_iterations;
  }

  template <typename T, typename... Args>
  void set_accountant(Args&&... args) {
    config_.accountant = std::make_unique<T>(std::forward<Args>(args)...);
  }

  template <typename T, typename... Args>
  void set_primal_pricing(Args&&... args) {
    config_.primal_pricing = std::make_unique<T>(std::forward<Args>(args)...);
  }

  template <typename T, typename... Args>
  void set_dual_pricing(Args&&... args) {
    config_.dual_pricing = std::make_unique<T>(std::forward<Args>(args)...);
  }

  // Point associated with the given states must be dual feasible
  SimplexResult<Field> dual(const std::vector<VariableState>& states) {
    validate([&] -> std::optional<std::string> {
      if (!is_dual_feasible(*problem_, states, config_.tolerance.feasibility)) {
        return "Initial point is not dual feasible.";
      }

      return validate_states(states);
    });

    try {
      return dual_implementation(states);
    } catch (...) {
      dump_state(*problem_, var_states_);
      throw;
    }
  }

  // Point associated with the given states must be primal feasible
  SimplexResult<Field> primal(const std::vector<VariableState>& states) {
    validate([&] -> std::optional<std::string> {
      if (!is_primal_feasible(*problem_, states,
                              config_.tolerance.feasibility)) {
        return "Initial point is not primal feasible.";
      }

      return validate_states(states);
    });

    try {
      return primal_implementation(states);
    } catch (...) {
      dump_state(*problem_, var_states_);
      throw;
    }
  }

  // traversal in problem space
  void change_basis(size_t leaving_index, size_t entering_variable,
                    VariableState leaving_state) {
    lupa_->change_column(leaving_index, entering_variable);

    var_states_[entering_variable] = VariableState::BASIC;

    var_states_[basic_vars_[leaving_index]] = leaving_state;
    basic_vars_[leaving_index] = entering_variable;
  }

  void change_bound(size_t variable_index, VariableState new_bound) {
    validate([&] -> std::optional<std::string> {
      if (new_bound != VariableState::AT_LOWER &&
          new_bound != VariableState::AT_UPPER) {
        return "new_bound must be either AT_LOWER or AT_UPPER.";
      }

      if (new_bound == VariableState::AT_LOWER &&
          !problem_->var_bounds[variable_index].lower) {
        return "new_bound is set to AT_LOWER, but variable doesn't have a "
               "lower bound.";
      }

      if (new_bound == VariableState::AT_UPPER &&
          !problem_->var_bounds[variable_index].upper) {
        return "new_bound is set to AT_UPPER, but variable doesn't have an "
               "upper bound.";
      }

      return std::nullopt;
    });

    var_states_[variable_index] = new_bound;
  }

  // current basis getters
  std::vector<size_t> get_basic_vars() const { return basic_vars_; }

  Vector<Field> get_point() const {
    const size_t n = var_states_.size();

    Vector<Field> result(n);

    for (size_t i = 0; i < basic_vars_.size(); ++i) {
      result[basic_vars_[i]] = basic_point_[i];
    }
    for (size_t i = 0; i < n; ++i) {
      switch (var_states_[i]) {
        case VariableState::AT_LOWER:
          result[i] = *problem_->var_bounds[i].lower;
          break;
        case VariableState::AT_UPPER:
          result[i] = *problem_->var_bounds[i].upper;
          break;
        case VariableState::NONBASIC_FREE:
          result[i] = 0;
          break;
        case VariableState::BASIC:
          break;
        default:
          throw std::runtime_error("Unknown variable state.");
      }
    }

    return result;
  }

  std::vector<VariableState> get_states() const { return var_states_; }

  Vector<Field> get_tableau_row(size_t row) const {
    return lupa_->get_row(row);
  }
};

}  // namespace simplex
