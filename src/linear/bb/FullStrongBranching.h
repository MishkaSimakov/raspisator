#pragma once

#include <vector>

#include "BaseAccountant.h"
#include "Settings.h"
#include "linear/bb/NodeRelativeLocation.h"
#include "linear/model/MILP.h"
#include "linear/simplex/BoundedSimplexMethod.h"
#include "utils/Accumulators.h"

template <typename Field, typename Accountant = BaseAccountant<Field>>
class FullStrongBranchingBranchAndBound {
  size_t total_nodes_count_;

  // Node is pushed into waiting list when feasible solution in it is found and
  // branching variable is calculated
  struct Node {
    size_t id;

    std::vector<Field> upper;
    std::vector<Field> lower;

    std::vector<VariableState> variables_states;
    size_t branch_variable;

    Field value;
    Field expected;

    Field branch_value;
    Field fractional;
  };

  simplex::BoundedSimplexMethod<Field> lp_solver_;

  std::vector<Node> waiting_;
  std::optional<Node> lifo_slot_;
  std::optional<FiniteMILPSolution<Field>> incumbent_;

  std::vector<VariableType> variables_;

  BranchAndBoundSettings<Field> settings_;

  // used to calculate average simplex iterations count
  ArithmeticMean<double> simplex_iterations_;

  Accountant accountant_;

  static bool is_integer(const Field& value) {
    auto floored = FieldTraits<Field>::floor(value);

    return !FieldTraits<Field>::is_nonzero(floored - value);
  }

  bool is_candidate_incumbent(const Matrix<Field>& point) const {
    for (size_t i = 0; i < variables_.size(); ++i) {
      if (variables_[i] == VariableType::INTEGER && !is_integer(point[i, 0])) {
        return false;
      }
    }

    return true;
  }

  bool try_prune(const FiniteLPSolution<Field>& solution, const Node& node) {
    if (is_candidate_incumbent(solution.point)) {
      if (!incumbent_ || incumbent_->value < solution.value) {
        incumbent_ = FiniteMILPSolution<Field>(solution.point, solution.value);
      }

      accountant_.pruned_by_integrality(node.id, solution);

      return true;
    }

    if (incumbent_ && !FieldTraits<Field>::is_strictly_positive(
                          solution.value - incumbent_->value)) {
      accountant_.pruned_by_bounds(node.id, solution);
      return true;
    }

    return false;
  }

  void register_simplex_fail(const std::vector<VariableState>& states,
                             const std::vector<Field>& lower,
                             const std::vector<Field>& upper) {
    size_t since_epoch = std::chrono::system_clock::now().time_since_epoch() /
                         std::chrono::milliseconds(1);
    std::string dump_name = std::format("simplex_core_dump_{}.h", since_epoch);

    lp_solver_.setup_warm_start(states);

    std::ofstream os(dump_name);
    lp_solver_.dump_state(os, since_epoch, lower, upper);

    std::println("Registered failed simplex run into {}.", dump_name);
  }

  Field merge_score(std::optional<Field> left_score,
                    std::optional<Field> right_score) {
    if (!left_score && !right_score) {
      return -1000;  // -inf
    }

    if (!left_score) {
      return -100 + *right_score;
    }

    if (!right_score) {
      return -100 + *left_score;
    }

    return (1 - settings_.score_factor) * std::min(*left_score, *right_score) +
           settings_.score_factor * std::max(*left_score, *right_score);
  }

  std::optional<Field> score_via_simplex(
      size_t index, const Node& node, const FiniteLPSolution<Field>& solution) {
    // full strong branching
    Field fractional = FieldTraits<Field>::fractional(solution.point[index, 0]);
    Field branch_value = FieldTraits<Field>::floor(solution.point[index, 0]);

    // calculate left node
    std::optional<Field> lower_score;
    auto tight_upper = node.upper;
    tight_upper[index] = branch_value;

    try {
      lp_solver_.setup_warm_start(node.variables_states);
      auto run_result = lp_solver_.dual(node.lower, tight_upper);

      if (run_result.is_feasible()) {
        lower_score =
            std::get<FiniteLPSolution<Field>>(run_result.solution).value -
            solution.value;
      }

      accountant_.strong_branching_simplex_run(node.id, index, run_result);
    } catch (...) {
      register_simplex_fail(node.variables_states, node.lower, tight_upper);
      throw;
    }

    // calculate right node
    std::optional<Field> upper_score;
    auto tight_lower = node.lower;
    tight_lower[index] = branch_value + 1;

    try {
      lp_solver_.setup_warm_start(node.variables_states);
      auto run_result = lp_solver_.dual(tight_lower, node.upper);

      if (run_result.is_feasible()) {
        upper_score =
            std::get<FiniteLPSolution<Field>>(run_result.solution).value -
            solution.value;
      }

      accountant_.strong_branching_simplex_run(node.id, index, run_result);
    } catch (...) {
      register_simplex_fail(node.variables_states, tight_lower, node.upper);
      throw;
    }

    return merge_score(lower_score, upper_score);
  }

  size_t find_branching_variable(const Node& node,
                                 const FiniteLPSolution<Field>& solution) {
    ArgMinimum<Field> score;

    // calculate pseudo cost score for each candidate
    for (size_t i = 0; i < variables_.size(); ++i) {
      if (variables_[i] == VariableType::INTEGER &&
          !is_integer(solution.point[i, 0])) {
        auto current_score = score_via_simplex(i, node, solution);
        score.record(i, current_score);
      }
    }

    if (!score.argmin()) {
      // return first violating variable
      for (size_t i = 0; i < variables_.size(); ++i) {
        if (variables_[i] == VariableType::INTEGER &&
            !is_integer(solution.point[i, 0])) {
          return i;
        }
      }
    }

    return *score.argmin();
  }

  void log_bounds(const Node& node,
                  const FiniteLPSolution<Field>& current_solution) {
    std::string lb = "*";

    if (incumbent_) {
      std::stringstream ss;
      ss << incumbent_->value;
      lb = ss.str();
    }

    std::println("#{}; LB: {}; UB: {}", node.id, lb, current_solution.value);
    std::cout << simplex_iterations_.mean() << std::endl;
  }

  void record_simplex_run(const Node& parent, NodeRelativeLocation location,
                          const SimplexResult<Field>& run_result) {
    assert(location != NodeRelativeLocation::ROOT);

    simplex_iterations_.record(run_result.iterations_count);
  }

  void try_push_to_waiting(Node node, const SimplexResult<Field>& run_result) {
    // prune
    if (std::holds_alternative<NoFeasibleElements>(run_result.solution)) {
      accountant_.pruned_by_infeasibility(node.id);
      return;
    }

    auto solution = std::get<FiniteLPSolution<Field>>(run_result.solution);

    if (try_prune(solution, node)) {
      return;
    }

    log_bounds(node, solution);

    // branch
    node.value = solution.value;
    node.variables_states = solution.variables;

    // TODO: expected value
    node.expected = node.value;

    node.branch_variable = find_branching_variable(node, solution);
    node.branch_value =
        FieldTraits<Field>::floor(solution.point[node.branch_variable, 0]);
    node.fractional =
        solution.point[node.branch_variable, 0] - node.branch_value;

    accountant_.set_branching_variable(node.id, solution, node.branch_variable,
                                       node.branch_value);

    if (!lifo_slot_) {
      lifo_slot_ = std::move(node);
    } else if (lifo_slot_->expected > node.expected) {
      waiting_.push_back(std::move(*lifo_slot_));
      lifo_slot_ = std::move(node);
    } else {
      waiting_.push_back(std::move(node));
    }
  }

  void calculate_node(const Node& parent, NodeRelativeLocation location) {
    assert(location != NodeRelativeLocation::ROOT);

    Node child;
    child.id = ++total_nodes_count_;

    accountant_.set_child(parent.id, child.id, location);

    // tighten bounds
    child.lower = parent.lower;
    child.upper = parent.upper;

    if (location == NodeRelativeLocation::LEFT_CHILD) {
      child.upper[parent.branch_variable] = parent.branch_value;
    } else {
      child.lower[parent.branch_variable] = parent.branch_value + 1;
    }

    // run simplex method for subproblem
    lp_solver_.setup_warm_start(parent.variables_states);
    SimplexResult<Field> run_result;

    try {
      run_result = lp_solver_.dual(child.lower, child.upper);
    } catch (...) {
      register_simplex_fail(parent.variables_states, child.lower, child.upper);
      throw;
    }

    accountant_.simplex_run(child.id, run_result);
    record_simplex_run(parent, location, run_result);

    try_push_to_waiting(std::move(child), run_result);
  }

  static Matrix<Field> perturbed_costs(
      const MILPProblem<Field>& problem,
      const BranchAndBoundSettings<Field>& settings) {
    if (settings.perturbation == PerturbationMode::DISABLED) {
      return problem.c;
    }

    auto [n, d] = problem.A.shape();
    Matrix<Field> result = problem.c;

    if (settings.perturbation == PerturbationMode::CONSTANT) {
      for (size_t i = 0; i < d; ++i) {
        if (!FieldTraits<Field>::is_nonzero(problem.c[0, i])) {
          result[0, i] = settings.perturbation_value;
        }
      }
    } else if (settings.perturbation ==
               PerturbationMode::FOR_INTEGER_SOLUTION) {
      std::vector<size_t> perturbed;

      for (size_t i = 0; i < d; ++i) {
        if (FieldTraits<Field>::is_nonzero(problem.c[0, i])) {
          continue;
        }

        if (FieldTraits<Field>::is_nonzero(problem.lower_bounds[i])) {
          continue;
        }

        if (!FieldTraits<Field>::is_nonzero(problem.lower_bounds[i] -
                                            problem.upper_bounds[i])) {
          continue;
        }

        result[0, i] = static_cast<Field>(1) / problem.upper_bounds[i];
        perturbed.push_back(i);
      }

      for (size_t i : perturbed) {
        result[0, i] /= static_cast<Field>(2 * perturbed.size());
      }
    }

    return result;
  }

  std::optional<Node> pop_node() {
    if (lifo_slot_) {
      auto result = std::move(*lifo_slot_);
      lifo_slot_ = std::nullopt;

      return result;
    }

    if (waiting_.empty()) {
      return std::nullopt;
    }

    // choose node with the lowest expected value
    ArgMinimum<Field> min_expected;
    for (size_t i = 0; i < waiting_.size(); ++i) {
      min_expected.record(i, waiting_[i].expected);
    }

    auto result = std::move(waiting_[*min_expected.argmin()]);
    waiting_.erase(waiting_.begin() + *min_expected.argmin());

    return result;
    // auto result = std::move(waiting_.back());
    // waiting_.pop_back();
    //
    // return result;
  }

 public:
  explicit FullStrongBranchingBranchAndBound(
      Matrix<Field> A, Matrix<Field> b, Matrix<Field> c,
      std::vector<Field> lower, std::vector<Field> upper,
      std::vector<VariableType> variables,
      BranchAndBoundSettings<Field> settings = {})
      : total_nodes_count_(0),
        lp_solver_(CSCMatrix(A), b, c),
        variables_(variables),
        settings_(settings) {
    Node root;

    root.id = ++total_nodes_count_;
    root.lower = lower;
    root.upper = upper;

    auto basic_variables = linalg::get_row_basis(linalg::transposed(A));
    lp_solver_.setup_warm_start(basic_variables);
    root.variables_states = lp_solver_.get_variables_states();

    accountant_.set_root(root.id);
    lifo_slot_ = std::move(root);
  }

  const Accountant& get_accountant() const { return accountant_; }

  MILPSolution<Field> solve() {
    // solve for root
    std::optional<Node> root = pop_node();
    assert(root.has_value());

    auto run_result = lp_solver_.dual(root->lower, root->upper);
    try_push_to_waiting(std::move(*root), run_result);

    // main branch and bound cycle
    std::optional<Node> parent;
    while ((parent = pop_node())) {
      calculate_node(*parent, NodeRelativeLocation::LEFT_CHILD);
      calculate_node(*parent, NodeRelativeLocation::RIGHT_CHILD);

      if (settings_.max_nodes && total_nodes_count_ > settings_.max_nodes) {
        return ReachedNodesLimit{};
      }
    }

    if (incumbent_) {
      return *incumbent_;
    }

    return NoFiniteSolution{};
  }
};
