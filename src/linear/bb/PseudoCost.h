#pragma once

#include <vector>

#include "BaseAccountant.h"
#include "Settings.h"
#include "linear/model/MILP.h"
#include "linear/simplex/BoundedSimplexMethod.h"

template <typename Field, typename Accountant = BaseAccountant<Field>>
class PseudoCostBranchAndBound {
  constexpr static double kInitialPseudoCost = 1;
  constexpr static double kReliabilityParameter = 4;
  constexpr static double kInfeasibilityScore = -1;

  size_t total_nodes_count_;

  struct Node {
    size_t id;

    std::vector<Field> upper;
    std::vector<Field> lower;

    std::vector<VariableState> variables_states;

    // parent and branching information
    Field parent_value;
    size_t branch_variable;
    Field fractional_part;
    bool is_lower;

    Field expected_value;
  };

  simplex::BoundedSimplexMethod<Field> lp_solver_;

  std::vector<Node> waiting_;
  std::optional<FiniteMILPSolution<Field>> incumbent_;

  std::vector<VariableType> variables_;

  std::vector<Field> upper_pseudocosts_;
  std::vector<Field> lower_pseudocosts_;
  std::vector<size_t> branch_count_;

  BranchAndBoundSettings<Field> settings_;

  // used to calculate average simplex iterations count
  size_t total_simplex_iterations_;
  size_t total_simplex_runs_;

  Accountant accountant_;

  size_t average_simplex_iterations() const {
    double average = static_cast<double>(total_simplex_iterations_) /
                     static_cast<double>(total_simplex_runs_);

    return std::ceil(average);
  }

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

  // double score_branching_variable(const Node& node, const Matrix<Field>&
  // point,
  //                                 size_t index) const {
  //   // most infeasible branching
  //   // branch on the variable, whose fractional value is closest to 0.5
  //   Field f = point[index, 0] - FieldTraits<Field>::floor(point[index, 0]);
  //   return FieldTraits<Field>::abs(0.5 - f);
  // }

  double merge_score(double left_score, double right_score) {
    constexpr double kScoreFactor = 1. / 6;

    return (1 - kScoreFactor) * std::min(left_score, right_score) +
           kScoreFactor * std::max(left_score, right_score);
  }

  double score_via_simplex(const Node& node,
                           const FiniteLPSolution<Field>& solution,
                           size_t index) {
    // full strong branching
    // run simplex method for gamma iterations for each subproblem,
    // score is the change of the objective
    auto score_function = Overload{
        [](NoFeasibleElements) { return kInfeasibilityScore; },
        [&solution](const FiniteLPSolution<Field>& subproblem_solution) {
          return static_cast<double>(subproblem_solution.value -
                                     solution.value);
        },
    };

    double fractional = static_cast<double>(
        solution.point[index, 0] -
        FieldTraits<Field>::floor(solution.point[index, 0]));

    // calculate left node
    double left_delta;
    auto tight_upper = node.upper;
    tight_upper[index] = FieldTraits<Field>::floor(solution.point[index, 0]);

    try {
      lp_solver_.setup_warm_start(node.variables_states);
      auto run_result = lp_solver_.dual(node.lower, tight_upper);
      left_delta = std::visit(score_function, run_result.solution);

      lower_pseudocosts_[index] += left_delta / fractional;

      accountant_.strong_branching_simplex_run(node.id, index, run_result);
    } catch (...) {
      lp_solver_.setup_warm_start(node.variables_states);

      std::ofstream os(std::format("simplex_core_dump_{}.h", 0));
      lp_solver_.dump_state(os, 0, node.lower, tight_upper);
      throw;
    }

    // calculate right node
    auto tight_lower = node.lower;
    tight_lower[index] =
        FieldTraits<Field>::floor(solution.point[index, 0]) + 1;

    double right_delta;
    try {
      lp_solver_.setup_warm_start(node.variables_states);
      auto run_result = lp_solver_.dual(tight_lower, node.upper);
      right_delta = std::visit(score_function, run_result.solution);

      upper_pseudocosts_[index] += right_delta / (1 - fractional);

      accountant_.strong_branching_simplex_run(node.id, index, run_result);
    } catch (...) {
      lp_solver_.setup_warm_start(node.variables_states);

      std::ofstream os(std::format("simplex_core_dump_{}.h", 0));
      lp_solver_.dump_state(os, 0, tight_lower, node.upper);
      throw;
    }

    ++branch_count_[index];

    return merge_score(left_delta, right_delta);
  }

  // size_t find_branching_variable(const Node& node, const Matrix<Field>&
  // point) {
  //   std::optional<std::pair<size_t, double>> min_score;
  //
  //   for (size_t i = 0; i < variables_.size(); ++i) {
  //     if (variables_[i] == VariableType::INTEGER && !is_integer(point[i, 0]))
  //     {
  //       double score = score_branching_variable(node, point, i);
  //
  //       if (!min_score || score < min_score->second) {
  //         min_score = {i, score};
  //       }
  //     }
  //   }
  //
  //   assert(min_score.has_value());
  //   return min_score->first;
  // }

  bool is_reliable(size_t variable_id) {
    return branch_count_[variable_id] >= kReliabilityParameter;
  }

  double score_via_pseudocost(size_t index, const Matrix<Field>& point) {
    double fractional = static_cast<double>(
        point[index, 0] - FieldTraits<Field>::floor(point[index, 0]));

    double lower_cost =
        fractional * lower_pseudocosts_[index] / branch_count_[index];
    double upper_cost =
        (1 - fractional) * upper_pseudocosts_[index] / branch_count_[index];

    return merge_score(lower_cost, upper_cost);
  }

  size_t find_branching_variable(const Node& node,
                                 const FiniteLPSolution<Field>& solution) {
    std::optional<std::pair<size_t, double>> min_score;

    // calculate pseudo cost score for each candidate
    for (size_t i = 0; i < variables_.size(); ++i) {
      if (variables_[i] == VariableType::INTEGER &&
          !is_integer(solution.point[i, 0])) {
        std::cout << branch_count_[i] << " ";

        double score = is_reliable(i) ? score_via_pseudocost(i, solution.point)
                                      : score_via_simplex(node, solution, i);

        if (!min_score || score < min_score->second) {
          min_score = {i, score};
        }
      }
    }

    std::cout << std::endl;
    assert(min_score.has_value());
    return min_score->first;
  }

  void partition(const FiniteLPSolution<Field>& solution, const Node& parent) {
    size_t branch_variable = find_branching_variable(parent, solution);
    Field branch_value =
        FieldTraits<Field>::floor(solution.point[branch_variable, 0]);

    accountant_.set_branching_variable(parent.id, solution, branch_variable,
                                       branch_value);

    Node left;
    left.id = ++total_nodes_count_;
    left.variables_states = solution.variables;
    left.lower = parent.lower;
    left.upper = parent.upper;
    left.parent_value = solution.value;

    left.upper[branch_variable] = branch_value;

    accountant_.set_left_child(parent.id, left.id);
    waiting_.push_back(left);

    Node right;
    right.id = ++total_nodes_count_;
    right.variables_states = solution.variables;
    right.lower = parent.lower;
    right.upper = parent.upper;
    right.parent_value = solution.value;

    right.lower[branch_variable] = branch_value + 1;

    accountant_.set_right_child(parent.id, right.id);
    waiting_.push_back(right);
  }

  void log_bounds(const Node& node,
                  const FiniteLPSolution<Field>& current_solution) {
    std::string lb = incumbent_ ? std::to_string(incumbent_->value) : "*";

    std::println("#{}; LB: {}; UB: {}", node.id, lb, current_solution.value);
  }

 public:
  explicit PseudoCostBranchAndBound(MILPProblem<Field> problem,
                                    BranchAndBoundSettings<Field> settings)
      : total_nodes_count_(0),
        lp_solver_(CSCMatrix(problem.A), problem.b, problem.c),
        variables_(problem.variables_),
        upper_pseudocosts_(problem.variables_.size(), kInitialPseudoCost),
        lower_pseudocosts_(problem.variables_.size(), kInitialPseudoCost),
        branch_count_(problem.variables_.size(), 0),
        settings_(settings) {
    Node root;

    root.id = ++total_nodes_count_;
    root.lower = problem.lower_bounds;
    root.upper = problem.upper_bounds;

    auto basic_variables = linalg::get_row_basis(linalg::transposed(problem.A));
    lp_solver_.setup_warm_start(basic_variables);
    root.variables_states = lp_solver_.get_variables_states();

    accountant_.set_root(root.id);
    waiting_.emplace_back(root);
  }

  const Accountant& get_accountant() const { return accountant_; }

  MILPSolution<Field> solve() {
    while (!waiting_.empty()) {
      if (settings_.max_nodes && total_nodes_count_ > settings_.max_nodes) {
        return ReachedNodesLimit{};
      }

      Node current = waiting_.back();
      waiting_.pop_back();
      lp_solver_.setup_warm_start(current.variables_states);

      auto run_result = lp_solver_.dual(current.lower, current.upper);
      accountant_.simplex_run(current.id, run_result);

      ++total_simplex_runs_;
      total_simplex_iterations_ += run_result.iterations_count;

      if (std::holds_alternative<NoFeasibleElements>(run_result.solution)) {
        accountant_.pruned_by_infeasibility(current.id);
        continue;
      }

      auto solution = std::get<FiniteLPSolution<Field>>(run_result.solution);

      if (try_prune(solution, current)) {
        continue;
      }

      log_bounds(current, solution);

      partition(solution, current);
    }

    if (incumbent_) {
      return *incumbent_;
    }

    return NoFiniteSolution{};
  }
};
