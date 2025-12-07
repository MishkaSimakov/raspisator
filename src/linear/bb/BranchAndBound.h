#pragma once

#include <optional>
#include <queue>

#include "../simplex/BoundedSimplexMethod.h"
#include "BranchAndBoundTree.h"
#include "linear/bb/Settings.h"
#include "linear/matrix/Matrix.h"
#include "linear/model/LP.h"
#include "linear/model/MILP.h"

template <typename Field>
struct NodeComparator {
  bool operator()(const InteriorNode<Field>* left,
                  const InteriorNode<Field>* right) const {
    return left->value < right->value;
  }
};

// TODO: implement more efficient version for binary variables
template <typename Field, typename LPSolver>
class BranchAndBound {
  const MILPProblem<Field> problem_;
  const BranchAndBoundSettings<Field> settings_;

  BranchAndBoundTree<Field> tree_;

  std::priority_queue<InteriorNode<Field>*, std::vector<InteriorNode<Field>*>,
                      NodeComparator<Field>>
      queue_;

  std::optional<std::pair<Field, Matrix<Field>>> lower_bound_;

  simplex::BoundedSimplexMethod<Field> lp_solver_;

  static Matrix<Field> floor_solution(const Matrix<Field>& point) {
    auto result = Matrix<Field>::zeros_like(point);

    for (size_t i = 0; i < point.get_height(); ++i) {
      result[i, 0] = FieldTraits<Field>::floor(point[i, 0]);
    }

    return result;
  }

  // returns branching variable index
  size_t find_branching(const Matrix<Field>& solution) const {
    size_t max_fractional_index = 0;
    Field max_fractional_value = 0;

    for (size_t i = 0; i < problem_.variables_.size(); ++i) {
      if (problem_.variables_[i] != VariableType::INTEGER) {
        continue;
      }

      Field fractional =
          solution[i, 0] - FieldTraits<Field>::floor(solution[i, 0]);

      if (fractional > max_fractional_value) {
        max_fractional_value = fractional;
        max_fractional_index = i;
      }
    }

    return max_fractional_index;
  }

  std::pair<std::vector<Field>, std::vector<Field>> accumulate_bounds(
      const InteriorNode<Field>* parent, NodeRelativeLocation type) {
    if (parent == nullptr) {
      return {problem_.lower_bounds, problem_.upper_bounds};
    }

    std::vector<Field> lower_bounds = problem_.lower_bounds;
    std::vector<Field> upper_bounds = problem_.upper_bounds;

    auto [path, directions] = tree_.get_path_to(parent);

    directions.push_back(type);

    for (size_t i = 0; i < path.size(); ++i) {
      const InteriorNode<Field>* node = path[i];

      if (directions[i] == NodeRelativeLocation::LEFT_CHILD) {
        upper_bounds[node->branching_variable] = node->branching_value;
      } else {
        lower_bounds[node->branching_variable] = node->branching_value + 1;
      }
    }

    return {lower_bounds, upper_bounds};
  }

  void calculate_node(InteriorNode<Field>* parent, NodeRelativeLocation type) {
    auto [lower_bounds, upper_bounds] = accumulate_bounds(parent, type);

    if (parent == nullptr) {
      auto columns_basis =
          linalg::get_row_basis(linalg::transposed(problem_.A));
      lp_solver_.setup_warm_start(columns_basis);
    } else {
      lp_solver_.setup_warm_start(parent->variables);
    }

    auto relaxed_solution = lp_solver_.dual(lower_bounds, upper_bounds);

    // prune branch if there is no feasible solution
    if (std::holds_alternative<NoFeasibleElements>(relaxed_solution)) {
      tree_.add_node(parent, type, UnfeasibleNode<Field>{});
      return;
    }

    auto& finite_solution = std::get<FiniteLPSolution<Field>>(relaxed_solution);
    auto& point = finite_solution.point;
    auto& value = finite_solution.value;

    size_t br_variable = find_branching(point);

    Field floored = FieldTraits<Field>::floor(point[br_variable, 0]);

    // prune branch if all values are integer
    if (!FieldTraits<Field>::is_strictly_positive(point[br_variable, 0] -
                                                  floored)) {
      // possibly update lower bound
      if (!lower_bound_ || value > lower_bound_->first) {
        size_t d = problem_.A.get_width();

        lower_bound_ = std::pair{value, point[{0, d}, 0]};
      }

      tree_.add_node(parent, type, IntegerSolutionNode{value});
      return;
    }

    // prune branch if upper bound is worse than current best solution
    if (lower_bound_ && value <= lower_bound_->first) {
      tree_.add_node(parent, type, BoundedNode{value});
      return;
    }

    auto& child = tree_.add_node(parent, type, InteriorNode<Field>{});

    child.variables = finite_solution.variables;
    child.branching_variable = br_variable;
    child.branching_value = floored;
    child.value = value;

    queue_.push(&child);
  }

 public:
  explicit BranchAndBound(MILPProblem<Field> problem,
                          BranchAndBoundSettings<Field> settings = {})
      : problem_(problem),
        settings_(settings),
        lp_solver_(CSCMatrix(problem.A), problem.b, problem.c) {}

  MILPSolution<Field> solve() {
    // solve problem
    calculate_node(nullptr, NodeRelativeLocation::ROOT);

    while (!queue_.empty()) {
      // node with maximum upper bound
      auto* current_node = queue_.top();
      queue_.pop();

      if (lower_bound_) {
        std::println("LB: {}; UB: {}", lower_bound_->first,
                     current_node->value);
      } else {
        std::println("LB: *; UB: {}", current_node->value);
      }

      calculate_node(current_node, NodeRelativeLocation::LEFT_CHILD);
      calculate_node(current_node, NodeRelativeLocation::RIGHT_CHILD);

      if (settings_.max_nodes.has_value() &&
          tree_.size() > *settings_.max_nodes) {
        return ReachedNodesLimit{};
      }
    }

    if (!lower_bound_) {
      return NoFiniteSolution{};
    }

    return FiniteMILPSolution{lower_bound_->second, lower_bound_->first};
  }

  const BranchAndBoundTree<Field>& get_tree() const { return tree_; }
};
