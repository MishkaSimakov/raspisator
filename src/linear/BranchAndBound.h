#pragma once

#include <optional>
#include <queue>
#include <ranges>
#include <sstream>
#include <unordered_map>

#include "linear/CheckBFS.h"
#include "linear/model/LP.h"
#include "linear/model/MILP.h"
#include "matrix/Matrix.h"
#include "matrix/RowBasis.h"

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

  Matrix<Field> solution;
  std::vector<size_t> basic_variables;
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

    std::string fill =
        node->calculated ? ", style=filled, fillcolor=lightblue" : "";

    ss_ << std::format("  {} [label=\"UB={}\\nx_{} <> {}\" {}];\n", id,
                       node->upper_bound, node->branching_variable,
                       node->branching_value, fill);

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
  enum class ConstraintType { UPPER, LOWER };

  struct Constraint {
    ConstraintType type;
    size_t variable_index;
    Field value;
  };

  // returns path from the root to the given node
  static std::vector<const Node<Field>*> get_path_from_root(
      const Node<Field>* node) {
    std::vector<const Node<Field>*> result;

    while (node != nullptr) {
      result.push_back(node);
      node = node->parent;
    }

    std::ranges::reverse(result);

    return result;
  }

  // returned vector contains constraints from root to child
  std::vector<Constraint> accumulate_constraints(const Node<Field>* parent,
                                                 NodeType type) {
    std::unordered_map<size_t, size_t> upper_constraints;
    std::unordered_map<size_t, size_t> lower_constraints;

    std::vector<Constraint> result;
    auto path = get_path_from_root(parent);

    for (size_t i = 0; i < path.size(); ++i) {
      const Node<Field>* node = path[i];

      bool is_left_child;
      if (i + 1 == path.size()) {
        is_left_child = type == NodeType::LEFT_CHILD;
      } else {
        is_left_child = path[i + 1] == node->left;
      }

      if (is_left_child) {
        if (upper_constraints.contains(node->branching_variable)) {
          result[upper_constraints.at(node->branching_variable)].value =
              node->branching_value;
        } else {
          upper_constraints[node->branching_variable] = result.size();
          result.emplace_back(ConstraintType::UPPER, node->branching_variable,
                              node->branching_value);
        }
      } else {
        if (lower_constraints.contains(node->branching_variable)) {
          result[lower_constraints.at(node->branching_variable)].value =
              node->branching_value + 1;
        } else {
          lower_constraints[node->branching_variable] = result.size();
          result.emplace_back(ConstraintType::LOWER, node->branching_variable,
                              node->branching_value + 1);
        }
      }
    }

    return result;
  }

  // returns an array of variable shifts that were applied
  // they are used to revert constraints
  std::tuple<Matrix<Field>, Matrix<Field>, Matrix<Field>, std::vector<Field>,
             std::optional<BFS<Field>>>
  apply_constraints(const Node<Field>* node, NodeType type) {
    auto [n, d] = problem_.A.shape();
    std::vector<Field> shifts(d, 0);

    if (type == NodeType::ROOT) {
      auto bfs = LPSolver(problem_.A, problem_.b, problem_.c).find_bfs();
      return {problem_.A, problem_.b, problem_.c, shifts, *bfs};
    }

    auto constraints = accumulate_constraints(node, type);

    size_t upper_count = 0;
    for (const Constraint& constraint : constraints) {
      if (constraint.type == ConstraintType::UPPER) {
        ++upper_count;
      }
    }

    auto A = problem_.A.get_extended(n + upper_count, d + upper_count, 0);
    auto b = problem_.b.get_extended(n + upper_count, 1, 0);
    auto c = problem_.c.get_extended(1, d + upper_count, 0);

    size_t slack_index = 0;

    for (const Constraint& constraint : constraints) {
      if (constraint.type == ConstraintType::UPPER) {
        A[n + slack_index, d + slack_index] = 1;
        A[n + slack_index, constraint.variable_index] = 1;

        b[n + slack_index, 0] = constraint.value;

        ++slack_index;
      }
    }

    for (const Constraint& constraint : constraints) {
      if (constraint.type == ConstraintType::LOWER) {
        shifts[constraint.variable_index] = -constraint.value;

        b[{0, b.get_height()}, 0].add_mul(
            A[{0, A.get_height()}, constraint.variable_index],
            -constraint.value);
      }
    }

    // update parent solution to form bfs for a child
    BFS<Field> bfs(node->solution.get_extended(d + upper_count, 1, 0),
                   node->basic_variables);

    size_t negative_index;

    if (type == NodeType::LEFT_CHILD) {
      if (d + upper_count > node->solution.get_height()) {
        // added new column and row
        // if it turns out that this row is linearly dependent with previous
        // rows then this branch is unfeasible
        if (linalg::get_row_basis(A).size() != n + upper_count) {
          return {A, b, c, shifts, std::nullopt};
        }

        // add it to basic variables
        bfs.basic_variables.push_back(d + upper_count - 1);
      }

      size_t slack_index = 0;
      for (const Constraint& constraint : constraints) {
        if (constraint.type == ConstraintType::UPPER) {
          if (constraint.variable_index == node->branching_variable) {
            break;
          }

          ++slack_index;
        }
      }

      bfs.point[d + slack_index, 0] =
          bfs.point[node->branching_variable, 0] - node->branching_value;
      negative_index = d + slack_index;
    } else {
      bfs.point[node->branching_variable, 0] +=
          shifts[node->branching_variable];
      negative_index = node->branching_variable;
    }

    // std::cout << "a:\n" << A << std::endl;
    // std::cout << "b:\n" << b << std::endl;
    // std::cout << "c:\n" << c << std::endl;
    //
    // std::cout << linalg::transposed(bfs.point) << std::endl;
    // for (size_t i : bfs.basic_variables) {
    //   std::cout << i << " ";
    // }
    // std::cout << std::endl;

    auto solver = LPSolver(A, b, c);
    auto reconstructed_bfs = solver.reconstruct_bfs(bfs, negative_index);

    return {A, b, c, shifts, std::move(reconstructed_bfs)};
  }

  enum class NodeCalculationResult {
    INFINITE_SOLUTION,
    NO_FEASIBLE_ELEMENTS,
    INTEGER_SOLUTION,
    WORSE_THAN_LOWER_BOUND,
    NORMAL
  };

  NodeCalculationResult calculate_node(Node<Field>* parent, NodeType type) {
    auto [A, b, c, shifts, bfs] = apply_constraints(parent, type);

    if (!bfs) {
      return NodeCalculationResult::NO_FEASIBLE_ELEMENTS;
    }

    auto relaxed_solution = LPSolver(A, b, c).solve_from(*bfs);

    // prune branch if there is no feasible solution
    if (std::holds_alternative<InfiniteSolution>(relaxed_solution)) {
      return NodeCalculationResult::INFINITE_SOLUTION;
    }

    FiniteLPSolution<Field>& finite_solution =
        std::get<FiniteLPSolution<Field>>(relaxed_solution);

    auto original_point = finite_solution.point;
    for (size_t i = 0; i < shifts.size(); ++i) {
      original_point[i, 0] -= shifts[i];
    }

    auto br_variable = find_branching(original_point);

    Field total_shift = 0;
    for (size_t i = 0; i < shifts.size(); ++i) {
      total_shift += shifts[i] * problem_.c[0, i];
    }

    Field floored = FieldTraits<Field>::floor(original_point[br_variable, 0]);
    Field fractional = original_point[br_variable, 0] - floored;

    Field upper_bound = finite_solution.value - total_shift;

    // prune branch if all values are integer
    if (!FieldTraits<Field>::is_strictly_positive(fractional)) {
      // possibly update lower bound
      if (!lower_bound_ || upper_bound > lower_bound_->first) {
        size_t d = problem_.A.get_width();

        lower_bound_ = std::pair{upper_bound, original_point[{0, d}, 0]};
      }

      return NodeCalculationResult::INTEGER_SOLUTION;
    }

    // prune branch if upper bound is worse than current best solution
    if (lower_bound_ && upper_bound <= lower_bound_->first) {
      return NodeCalculationResult::WORSE_THAN_LOWER_BOUND;
    }

    Node<Field>& child = nodes_.emplace_back();

    child.solution = finite_solution.point;
    child.basic_variables = finite_solution.basic_variables;

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
    // solve problem
    calculate_node(nullptr, NodeType::ROOT);

    while (!queue_.empty()) {
      // node with maximum upper bound
      auto* current_node = queue_.top();
      queue_.pop();

      if (lower_bound_) {
        std::println("LB: {}; UB: {}", lower_bound_->first,
                     current_node->upper_bound);
      } else {
        std::println("LB: *; UB: {}", current_node->upper_bound);
      }

      current_node->calculated = true;

      calculate_node(current_node, NodeType::LEFT_CHILD);
      calculate_node(current_node, NodeType::RIGHT_CHILD);
    }

    if (!lower_bound_) {
      return NoFiniteSolution{};
    }

    return FiniteMILPSolution{lower_bound_->second, lower_bound_->first};
  }

  const Node<Field>* get_root() const { return &nodes_[0]; }
};
