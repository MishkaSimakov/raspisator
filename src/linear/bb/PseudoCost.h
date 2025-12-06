#pragma once

#include <queue>
#include <vector>

#include "Settings.h"
#include "linear/BoundedSimplexMethod.h"
#include "linear/bb/Accountant.h"
#include "linear/model/MILP.h"

template <typename Field, typename Accountant = BaseAccountant<Field>>
class PseudoCostBranchAndBound {
  size_t total_nodes_count_;

  struct Node {
    size_t id;

    std::vector<Field> upper;
    std::vector<Field> lower;

    std::vector<VariableState> variables_states;

    Field expected_value;
  };

  BoundedSimplexMethod<Field> lp_solver_;

  std::vector<Node> waiting_;
  std::optional<FiniteMILPSolution<Field>> incumbent_;

  std::vector<VariableType> variables_;
  std::vector<Field> upper_pseudocosts_;
  std::vector<Field> lower_pseudocosts_;

  BranchAndBoundSettings<Field> settings_;

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

  // double score_branching_variable(const Node& node, const Matrix<Field>& point,
  //                                 size_t index) const {
  //   // most infeasible branching
  //   // branch on the variable, whose fractional value is closest to 0.5
  //   Field f = point[index, 0] - FieldTraits<Field>::floor(point[index, 0]);
  //   return FieldTraits<Field>::abs(0.5 - f);
  // }

  double merge_score(std::optional<double> left_score,
                     std::optional<double> right_score) {
    if (!left_score && !right_score) {
      // -infinity
      return -1e10;
    }

    if (!left_score && right_score) {
      return *right_score;
    }

    if (left_score && !right_score) {
      return *left_score;
    }

    constexpr double kScoreFactor = 1. / 6;

    return (1 - kScoreFactor) * std::min(*left_score, *right_score) +
           kScoreFactor * std::max(*left_score, *right_score);
  }

  double score_branching_variable(const Node& node, const Matrix<Field>&
  point,
                                  size_t index) {
    // full strong branching
    auto score_function = Overload{
        [](NoFeasibleElements) -> std::optional<Field> { return std::nullopt;
        },
        [](const FiniteLPSolution<Field>& solution) -> std::optional<Field> {
          return solution.value;
        },
    };

    // calculate left node
    auto tight_upper = node.upper;
    tight_upper[index] = FieldTraits<Field>::floor(point[index, 0]);

    lp_solver_.setup_warm_start(node.variables_states);
    std::optional<double> left_score =
        std::visit(score_function, lp_solver_.dual(node.lower, tight_upper));

    // calculate right node
    auto tight_lower = node.lower;
    tight_lower[index] = FieldTraits<Field>::floor(point[index, 0]) + 1;

    lp_solver_.setup_warm_start(node.variables_states);
    std::optional<double> right_score =
        std::visit(score_function, lp_solver_.dual(tight_lower, node.upper));

    return merge_score(left_score, right_score);
  }

  size_t find_branching_variable(const Node& node, const Matrix<Field>& point) {
    std::optional<std::pair<size_t, double>> min_score;

    for (size_t i = 0; i < variables_.size(); ++i) {
      if (variables_[i] == VariableType::INTEGER && !is_integer(point[i, 0])) {
        double score = score_branching_variable(node, point, i);

        if (!min_score || score < min_score->second) {
          min_score = {i, score};
        }
      }
    }

    assert(min_score.has_value());
    return min_score->first;
  }

  void partition(const FiniteLPSolution<Field>& solution, const Node& parent) {
    size_t branch_variable = find_branching_variable(parent, solution.point);
    Field branch_value =
        FieldTraits<Field>::floor(solution.point[branch_variable, 0]);

    accountant_.set_branching_variable(parent.id, solution, branch_variable,
                                       branch_value);

    Node left;
    left.id = ++total_nodes_count_;
    left.variables_states = solution.variables;
    left.lower = parent.lower;
    left.upper = parent.upper;

    left.upper[branch_variable] = branch_value;

    accountant_.set_left_child(parent.id, left.id);
    waiting_.push_back(left);

    Node right;
    right.id = ++total_nodes_count_;
    right.variables_states = solution.variables;
    right.lower = parent.lower;
    right.upper = parent.upper;

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
      auto relaxed_solution = lp_solver_.dual(current.lower, current.upper);

      if (std::holds_alternative<NoFeasibleElements>(relaxed_solution)) {
        accountant_.pruned_by_infeasibility(current.id);
        continue;
      }

      if (try_prune(std::get<FiniteLPSolution<Field>>(relaxed_solution),
                    current)) {
        continue;
      }

      log_bounds(current, std::get<FiniteLPSolution<Field>>(relaxed_solution));

      partition(std::get<FiniteLPSolution<Field>>(relaxed_solution), current);
    }

    if (incumbent_) {
      return *incumbent_;
    }

    return NoFiniteSolution{};
  }
};
