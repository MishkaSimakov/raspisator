#pragma once

#include <optional>
#include <queue>
#include <ranges>
#include <sstream>
#include <unordered_map>

#include "linear/model/LP.h"
#include "linear/model/MILP.h"
#include "matrix/Matrix.h"

// Branch and Bound nodes
template <typename Field>
struct Node {
  bool calculated = false;

  Field upper_bound;

  size_t branching_variable;
  Field branching_value;

  Node* parent = nullptr;

  Node* left = nullptr;
  Node* right = nullptr;
};

template <typename Field>
struct NodeComparator {
  bool operator()(const Node<Field>* left, const Node<Field>* right) const {
    return left->upper_bound < right->upper_bound;
  }
};

// draws branch and bound tree using graphviz
template <typename Field>
class GraphvizBuilder {
 public:
  std::string build(const Node<Field>* root) {
    ss_ << "digraph G {\n";
    ss_ << "  node [shape=box];\n";

    if (root != nullptr) {
      add_node(root);
    }

    ss_ << "}\n";
    return ss_.str();
  }

 private:
  std::stringstream ss_;
  std::unordered_map<const Node<Field>*, std::string> ids_;

  std::string getId(const Node<Field>* ptr) {
    auto it = ids_.find(ptr);
    if (it != ids_.end()) {
      return it->second;
    }

    std::string id = "n" + std::to_string(ids_.size());
    ids_[ptr] = id;
    return id;
  }

  void add_node(const Node<Field>* node) {
    std::string id = getId(node);

    ss_ << std::format("  {} [label=\"UB={}\\nx_{} <> {}\"];\n", id,
                       node->upper_bound, node->branching_variable,
                       node->branching_value);

    if (node->left) {
      ss_ << std::format("  {} -> {} [label=\"left\"];\n", id,
                         getId(node->left));
      add_node(node->left);
    }

    if (node->right) {
      ss_ << std::format("  {} -> {} [label=\"right\"];\n", id,
                         getId(node->right));
      add_node(node->right);
    }
  }
};

// TODO: implement more efficient version for binary variables
template <typename Field, LPSolver<Field> LPSolver>
class BranchAndBound {
  const MILPProblem<Field> problem_;

  // stores branch and bound tree nodes
  // root node is nodes_[0]
  std::deque<Node<Field>> nodes_;

  // maps integer variables indices to corresponding slack variables
  std::unordered_map<size_t, std::pair<size_t, size_t>> slack_variables_;

  std::priority_queue<Node<Field>*, std::vector<Node<Field>*>,
                      NodeComparator<Field>>
      queue_;

  std::optional<std::pair<Field, Matrix<Field>>> lower_bound_;

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

  enum class NodeType { LEFT_CHILD, RIGHT_CHILD, ROOT };

  // returns an array of variable shifts that were applied
  // they are used to revert constraints
  std::tuple<Matrix<Field>, Matrix<Field>, Matrix<Field>, std::vector<Field>>
  apply_constraints(const Node<Field>* node, NodeType type) {
    auto [n, d] = problem_.A.shape();
    std::vector<Field> shifts(d, 0);

    if (type == NodeType::ROOT) {
      return {problem_.A, problem_.b, problem_.c, shifts};
    }

    struct UpperBound {
      size_t slack_variable_index;
      Field upper_bound;
    };

    std::unordered_map<size_t, UpperBound> upper_bounds;
    // lower bounds are stored in shifts

    const auto* curr = node;
    NodeType curr_type = type;

    while (true) {
      size_t br_variable = curr->branching_variable;

      if (curr_type == NodeType::LEFT_CHILD) {
        // apply left child constraints
        // x_i <= a iff x_i + lambda_i = a iff y_i + lambda_i = a + shifts[i]
        auto itr = upper_bounds.find(br_variable);
        if (itr == upper_bounds.end()) {
          auto bound = UpperBound{upper_bounds.size(), curr->branching_value};
          upper_bounds.emplace(br_variable, bound);
        } else {
          itr->second.upper_bound =
              std::min(itr->second.upper_bound, curr->branching_value);
        }
      } else if (curr_type == NodeType::RIGHT_CHILD) {
        // apply right child constraints
        // x_i >= a + 1 iff x_i - a - 1 >= 0 iff y_i >= 0,
        // where y_i = x_i - a - 1
        shifts[curr->branching_variable] = std::min(
            shifts[curr->branching_variable], -curr->branching_value - 1);
      }

      auto* next = curr->parent;

      if (next == nullptr) {
        break;
      }

      curr_type =
          next->left == curr ? NodeType::LEFT_CHILD : NodeType::RIGHT_CHILD;
      curr = next;
    }

    // apply upper bounds
    auto A = problem_.A.get_extended(n + upper_bounds.size(),
                                     d + upper_bounds.size(), 0);
    auto b = problem_.b.get_extended(n + upper_bounds.size(), 1, 0);
    auto c = problem_.c.get_extended(1, d + upper_bounds.size(), 0);

    for (const auto& [var_index, info] : upper_bounds) {
      size_t slack_index = info.slack_variable_index;
      A[n + slack_index, d + slack_index] = 1;
      A[n + slack_index, var_index] = 1;

      b[n + slack_index, 0] = info.upper_bound;
    }

    // apply shifts
    for (size_t i = 0; i < shifts.size(); ++i) {
      for (size_t j = 0; j < A.get_height(); ++j) {
        b[j, 0] += A[j, i] * shifts[i];
      }
    }

    return {A, b, c, shifts};
  }

  enum class NodeCalculationResult {
    INFINITE_SOLUTION,
    NO_FEASIBLE_ELEMENTS,
    INTEGER_SOLUTION,
    WORSE_THAN_LOWER_BOUND,
    NORMAL
  };

  NodeCalculationResult calculate_node(Node<Field>* parent, NodeType type) {
    auto [A, b, c, shifts] = apply_constraints(parent, type);

    // std::cout << "c:\n" << c << std::endl;
    // std::cout << "A | b:\n" << linalg::hstack(A, b) << std::endl;

    auto relaxed_solution = LPSolver(A, b, c).solve();

    // prune branch if there is no feasible solution
    if (std::holds_alternative<InfiniteSolution>(relaxed_solution)) {
      return NodeCalculationResult::INFINITE_SOLUTION;
    }
    if (std::holds_alternative<NoFeasibleElements>(relaxed_solution)) {
      return NodeCalculationResult::NO_FEASIBLE_ELEMENTS;
    }

    auto& finite_solution = std::get<FiniteLPSolution<Field>>(relaxed_solution);

    // std::cout << linalg::transposed(finite_solution.point) << std::endl;
    for (size_t i = 0; i < shifts.size(); ++i) {
      finite_solution.point[i, 0] -= shifts[i];
    }

    auto br_variable = find_branching(finite_solution.point);

    Field total_shift = 0;
    for (size_t i = 0; i < shifts.size(); ++i) {
      total_shift += shifts[i] * problem_.c[0, i];
    }

    Field floored =
        FieldTraits<Field>::floor(finite_solution.point[br_variable, 0]);
    Field fractional = finite_solution.point[br_variable, 0] - floored;

    Field upper_bound = finite_solution.value - total_shift;

    // prune branch if all values are integer
    if (fractional == 0 &&
        (!lower_bound_ || lower_bound_->first < upper_bound)) {
      size_t d = problem_.A.get_width();

      Matrix<Field> trimmed_solution(d, 1);
      for (size_t i = 0; i < d; ++i) {
        trimmed_solution[i, 0] = finite_solution.point[i, 0];
      }

      lower_bound_ = std::pair{upper_bound, trimmed_solution};
      return NodeCalculationResult::INTEGER_SOLUTION;
    }

    // prune branch if upper bound is worse than current best solution
    if (lower_bound_ && upper_bound <= lower_bound_->first) {
      return NodeCalculationResult::WORSE_THAN_LOWER_BOUND;
    }

    Node<Field>& child = nodes_.emplace_back();

    child.parent = parent;
    child.branching_variable = br_variable;
    child.branching_value = floored;
    child.upper_bound = upper_bound;

    if (type == NodeType::LEFT_CHILD) {
      parent->left = &child;
    } else if (type == NodeType::RIGHT_CHILD) {
      parent->right = &child;
    }

    queue_.push(&child);

    return NodeCalculationResult::NORMAL;
  }

 public:
  explicit BranchAndBound(MILPProblem<Field> problem) : problem_(problem) {}

  MILPSolution<Field> solve() {
    calculate_node(nullptr, NodeType::ROOT);

    while (!queue_.empty()) {
      std::cout << GraphvizBuilder<Field>().build(&nodes_[0]) << std::endl;

      // node with maximum upper bound
      auto* current_node = queue_.top();
      queue_.pop();

      current_node->calculated = true;

      std::cout << "left child" << std::endl;
      calculate_node(current_node, NodeType::LEFT_CHILD);
      std::cout << "right child" << std::endl;
      calculate_node(current_node, NodeType::RIGHT_CHILD);
    }

    if (!lower_bound_) {
      return NoFiniteSolution{};
    }

    return FiniteMILPSolution{lower_bound_->second, lower_bound_->first};
  }

  const Node<Field>* get_root() const { return &nodes_[0]; }
};
