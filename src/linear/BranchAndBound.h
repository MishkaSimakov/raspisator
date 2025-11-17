#pragma once
#include <fmt/format.h>

#include <queue>
#include <ranges>
#include <sstream>
#include <unordered_map>

#include "BranchAndBound.h"
#include "Matrix.h"
#include "SimplexMethod.h"
#include "utils/Variant.h"

// c x -> max, s.t. Ax = b, x >= 0
// and all integer_indices variables are integer
template <typename Field>
struct MILPProblem {
  Matrix<Field> A;
  Matrix<Field> b;
  Matrix<Field> c;

  std::vector<size_t> integer_indices;
};

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
  std::unordered_map<const void*, std::string> ids_;

  std::string getId(const void* ptr) {
    auto it = ids_.find(ptr);
    if (it != ids_.end()) {
      return it->second;
    }

    std::string id = "n" + std::to_string(ids_.size());
    ids_[ptr] = id;
    return id;
  }

  void add_node(Node<Field>* node_ptr) {
    std::string id = getId(node_ptr);

    std::visit(
        Overload{[id, this](const Node<Field>& node) {
                   ss_ << fmt::format(
                       "  {} [label=\"UB={}\\nLB={}\\nx_{} <> {}\"];\n", id,
                       node.upper_bound, node.lower_bound,
                       node.branching_variable, node.branching_value);

                   if (node.left) {
                     ss_ << fmt::format("  {} -> {};\n", id, getId(node.left));
                     add_node(node.left);
                   }

                   if (node.right) {
                     ss_ << fmt::format("  {} -> {};\n", id, getId(node.right));
                     add_node(node.right);
                   }
                 },
                 [id, this](const InfeasibleNode&) {
                   ss_ << "  " << id
                       << " [label=\"Infeasible\", style=filled, "
                          "fillcolor=lightgray];\n";
                 }},
        *node_ptr);
  }
};

template <typename Field, typename LPSolver>
class BranchAndBound {
  MILPProblem<Field> problem_;

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
  static size_t find_branching(const Matrix<Field>& solution) {
    size_t max_fractional_index = 0;
    Field max_fractional_value = 0;

    for (size_t i = 0; i < solution.get_height(); ++i) {
      Field fractional =
          solution[i, 0] - FieldTraits<Field>::floor(solution[i, 0]);

      if (fractional > max_fractional_value) {
        max_fractional_value = fractional;
        max_fractional_index = i;
      }
    }

    return max_fractional_index;
  }

  // for each integer variable adds new slack variable that will be used
  // during branch and bound algorithm
  void add_slack_variables() {
    auto [n, d] = problem_.A.shape();
    size_t integer_cnt = problem_.integer_indices.size();

    auto new_A = problem_.A.get_extended(n + integer_cnt, d + integer_cnt, 0);
    auto new_b = problem_.b.get_extended(n + integer_cnt, 1, 0);
    auto new_c = problem_.c.get_extended(1, d + integer_cnt, 0);

    for (size_t i = 0; i < integer_cnt; ++i) {
      new_A[n + i, d + i] = -1;
      new_A[n + i, problem_.integer_indices[i]] = 1;
    }

    problem_.A = std::move(new_A);
    problem_.b = std::move(new_b);
    problem_.c = std::move(new_c);

    for (size_t i = 0; i < integer_cnt; ++i) {
      slack_variables_[problem_.integer_indices[i]] = std::pair{i + n, i + d};
    }
  }

  void shift_variable(size_t index, Field delta) {
    for (size_t i = 0; i < problem_.A.get_height(); ++i) {
      problem_.b[i, 0] += problem_.A[i, index] * delta;
    }
  }

  // Each variable can be shifted. shifts[i] = a means that x_i was replaced
  // by y_i = x_i + a
  Node<Field>* solve_recursive(std::vector<Field> shifts, size_t depth) {
    if (depth > 1) {
      return &nodes_.emplace_back(InfeasibleNode{});
    }

    std::println("depth: {}", depth);

    auto& node = nodes_.emplace_back(step(problem_, shifts));
    Node<Field>* standard_node = std::get_if<Node<Field>>(&node);

    // we end traversal of a branch in two cases:
    // 1. there are no feasible points
    // 2. lower bound is equal to upper bound
    if (standard_node == nullptr ||
        standard_node->lower_bound == standard_node->upper_bound) {
      return &node;
    }

    auto [slack_var_row, slack_var_col] =
        slack_variables_[standard_node->branching_variable];

    // traverse left child
    // x_i <= a iff x_i + lambda_i = a iff y_i + lambda_i = a + shifts[i]
    {
      // remember old values before changing them
      Field A_coef = problem_.A[slack_var_row, slack_var_col];
      Field b_coef = problem_.b[slack_var_row, 0];

      problem_.A[slack_var_row, slack_var_col] = 1;
      problem_.b[slack_var_row, 0] = standard_node->branching_value +
                                     shifts[standard_node->branching_variable];

      standard_node->left = solve_recursive(shifts, depth + 1);

      problem_.A[slack_var_row, slack_var_col] = A_coef;
      problem_.b[slack_var_row, 0] = b_coef;
    }

    // traverse right child
    // x_i >= a + 1 iff x_i - a - 1 >= 0 iff y_i >= 0, where y_i = x_i - a - 1

    Field delta = -standard_node->branching_value - 1 -
                  shifts[standard_node->branching_variable];
    shift_variable(standard_node->branching_variable, delta);
    shifts[standard_node->branching_variable] += delta;

    standard_node->right = solve_recursive(shifts, depth + 1);

    shift_variable(standard_node->branching_variable, -delta);
    shifts[standard_node->branching_variable] -= delta;

    return &node;
  }

  // returns an array of variable shifts that were applied
  // they are used to revert constraints
  std::vector<Field> apply_constraints(const Node<Field>* node) {
    std::vector<Field> shifts(problem_.A.get_width(), 0);

    const auto* curr = node;
    const auto* prev = node->parent;

    while (prev != nullptr) {
      if (curr == std::get_if<Node<Field>>(prev->left)) {
        // apply left child constraints
        // x_i <= a iff x_i + lambda_i = a iff y_i + lambda_i = a + shifts[i]
        auto [row, col] = slack_variables_[curr->branching_variable];

        problem_.A[row, col] = 1;
        problem_.b[row, 0] = curr->branching_value;
      } else {
        // apply right child constraints
        // x_i >= a + 1 iff x_i - a - 1 >= 0 iff y_i >= 0,
        // where y_i = x_i - a - 1
        // we postpone application of shifts, just remember they must be done
        Field delta = -curr->branching_value - 1;
        shifts[curr->branching_variable] += delta;
      }
    }

    // apply shifts
    for (size_t i = 0; i < shifts.size(); ++i) {
      shift_variable(i, shifts[i]);
    }

    return shifts;
  }

  void remove_constraints(const std::vector<Field>& shifts) {
    for (size_t i = 0; i < shifts.size(); ++i) {
      shift_variable(i, -shifts[i]);
    }

    for (auto [row, col] : slack_variables_ | std::views::values) {
      problem_.A[row, col] = -1;
      problem_.b[row, 0] = 0;
    }
  }

  void calculate_node(Node<Field>* node, const std::vector<Field>& shifts) {
    auto relaxed_solution =
        LPSolver(problem_.A, problem_.b, problem_.c).solve();
    FiniteSolution<Field>* finite_solution =
        std::get_if<FiniteSolution<Field>>(&relaxed_solution);

    if (finite_solution == nullptr) {
      *node = InfeasibleNode{};
      return;
    }

    auto br_variable = find_branching(finite_solution->point);

    Field total_shift = 0;
    for (size_t i = 0; i < shifts.size(); ++i) {
      total_shift += shifts[i] * problem_.c[0, i];
    }

    Field floored =
        FieldTraits<Field>::floor(finite_solution->point[br_variable, 0]);
    Field fractional = finite_solution->point[br_variable, 0] - floored;

    Node<Field> result;

    result.branching_variable = br_variable;
    result.branching_value = floored - shifts[br_variable];
    result.upper_bound = finite_solution->value - total_shift;

    // if all values are integer we can update lower bound
    if (fractional == 0 &&
        (!lower_bound_ || lower_bound_->first < result.upper_bound)) {
      lower_bound_ = std::pair{result.upper_bound, finite_solution->point};
      return;
    }

    *node = std::move(result);
    queue_.push(std::get_if<Node<Field>>(node));
  }

 public:
  explicit BranchAndBound(MILPProblem<Field> problem) : problem_(problem) {}

  std::pair<Field, Matrix<Field>> solve() {
    add_slack_variables();

    std::vector<Field> shifts(problem_.A.get_width(), 0);

    auto& root = nodes_.emplace_back(Node<Field>{});
    calculate_node(&root, shifts);

    while (true) {
      // node with maximum upper bound
      auto* current_node = queue_.top();
      queue_.pop();

      if (lower_bound_ && current_node->upper_bound <= lower_bound_->first) {
        return *lower_bound_;
      }

      // left branch
      {
        current_node->left =
            &nodes_.emplace_back(Node<Field>{.parent = current_node});

        auto shifts =
            apply_constraints(std::get_if<Node<Field>>(current_node->left));
        calculate_node(current_node->left, shifts);
        remove_constraints(shifts);
      }

      // right branch
      {
        current_node->right =
            &nodes_.emplace_back(Node<Field>{.parent = current_node});

        auto shifts =
            apply_constraints(std::get_if<Node<Field>>(current_node->right));
        calculate_node(current_node->right, shifts);
        remove_constraints(shifts);
      }
    }
  }

  const Node<Field>* get_root() const { return &nodes_[0]; }
};
