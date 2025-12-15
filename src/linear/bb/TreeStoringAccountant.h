#pragma once

#include <cassert>
#include <optional>
#include <sstream>
#include <unordered_map>
#include <variant>

#include "linear/bb/BaseAccountant.h"

template <typename Field>
struct TreeStoringAccountant : BaseAccountant<Field> {
  struct PrunedByInfeasibility {};
  struct PrunedByIntegrality {
    Field value;
  };
  struct PrunedByBounds {
    Field value;
  };
  struct Branched {
    Field value;

    size_t branch_variable;
    Field branch_value;
  };
  struct Unvisited {};

  struct Node {
    std::variant<PrunedByInfeasibility, PrunedByIntegrality, PrunedByBounds,
                 Branched, Unvisited>
        state;

    std::optional<size_t> left_child = std::nullopt;
    std::optional<size_t> right_child = std::nullopt;

    std::optional<size_t> simplex_iterations = std::nullopt;
    std::vector<std::pair<size_t, size_t>> strong_branching_iterations;
  };

  std::unordered_map<size_t, Node> tree;
  std::optional<size_t> root;

  // graphviz methods
 private:
  static std::string get_node_graphviz_attributes(PrunedByInfeasibility,
                                                  bool minimalistic) {
    if (minimalistic) {
      return "style=filled, fillcolor=red";
    } else {
      return "label=\"unfeasible\", style=filled, fillcolor=red";
    }
  }

  static std::string get_node_graphviz_attributes(
      const PrunedByIntegrality& state, bool minimalistic) {
    if (minimalistic) {
      return "style=filled, fillcolor=green)";
    } else {
      return std::format(
          R"(label="integer\n{}", style=filled, fillcolor=green)", state.value);
    }
  }

  static std::string get_node_graphviz_attributes(const PrunedByBounds& state,
                                                  bool minimalistic) {
    if (minimalistic) {
      return "style=filled, fillcolor=orange";
    } else {
      return std::format(
          R"(label="bounded\n{}", style=filled, fillcolor=orange)",
          state.value);
    }
  }

  static std::string get_node_graphviz_attributes(const Branched& state,
                                                  bool minimalistic) {
    if (minimalistic) {
      return "";
    } else {
      return std::format(R"(label="{}\nx_{} <> {}")", state.value,
                         state.branch_variable, state.branch_value);
    }
  }

  static std::string get_node_graphviz_attributes(Unvisited,
                                                  bool minimalistic) {
    if (minimalistic) {
      return "style=filled, fillcolor=gray";
    } else {
      return "label=\"unvisited\", style=filled, fillcolor=gray";
    }
  }

  void node_to_graphviz(std::stringstream& os, size_t id,
                        bool minimalistic) const {
    const Node& node = tree.at(id);

    std::string attrs = std::visit(
        [minimalistic](const auto& node) {
          return get_node_graphviz_attributes(node, minimalistic);
        },
        node.state);
    std::println(os, "  {} [{}]", id, attrs);

    if (node.left_child) {
      if (minimalistic) {
        std::println(os, "  {} -> {};", id, *node.left_child);
      } else {
        std::println(os, "  {} -> {} [label=\"left\"];", id, *node.left_child);
      }
      node_to_graphviz(os, *node.left_child, minimalistic);
    }

    if (node.right_child) {
      if (minimalistic) {
        std::println(os, "  {} -> {};", id, *node.right_child);
      } else {
        std::println(os, "  {} -> {} [label=\"right\"];", id,
                     *node.right_child);
      }
      node_to_graphviz(os, *node.right_child, minimalistic);
    }
  }

  void iterations_recursive(std::ostream& os, size_t id, size_t depth) const {
    const Node& node = tree.at(id);

    std::string iterations_cnt = node.simplex_iterations
                                     ? std::to_string(*node.simplex_iterations)
                                     : "null";

    std::println(os, "{},{},{}", id, depth, iterations_cnt);

    if (node.left_child) {
      iterations_recursive(os, *node.left_child, depth + 1);
    }

    if (node.right_child) {
      iterations_recursive(os, *node.right_child, depth + 1);
    }
  }

  void strong_branching_iterations_recursive(std::ostream& os, size_t id,
                                             size_t depth) const {
    const Node& node = tree.at(id);

    for (auto [variable, cnt] : node.strong_branching_iterations) {
      std::println(os, "{},{},{},{}", id, depth, variable, cnt);
    }

    if (node.left_child) {
      strong_branching_iterations_recursive(os, *node.left_child, depth + 1);
    }

    if (node.right_child) {
      strong_branching_iterations_recursive(os, *node.right_child, depth + 1);
    }
  }

  void node_to_csv(std::stringstream& os, size_t id) const {
    const Node& node = tree.at(id);

    std::string left =
        node.left_child ? std::to_string(*node.left_child) : "null";
    std::string right =
        node.right_child ? std::to_string(*node.right_child) : "null";

    std::string value = std::visit(
        []<typename T>(const T& state) -> std::string {
          if constexpr (std::same_as<T, Unvisited> ||
                        std::same_as<T, PrunedByInfeasibility>) {
            return "null";
          } else {
            return std::to_string(state.value);
          }
        },
        node.state);

    std::string state = std::visit(
        Overload{
            [](PrunedByInfeasibility) { return "pruned_by_infeasibility"; },
            [](PrunedByIntegrality) { return "pruned_by_integrality"; },
            [](PrunedByBounds) { return "pruned_by_bounds"; },
            [](Branched) { return "branched"; },
            [](Unvisited) { return "unvisited"; },
        },
        node.state);

    std::println(os, "{},{},{},{},{}", id, value, state, left, right);

    if (node.left_child) {
      node_to_csv(os, *node.left_child);
    }

    if (node.right_child) {
      node_to_csv(os, *node.right_child);
    }
  }

 public:
  std::string iterations_to_csv() const {
    std::stringstream os;

    std::println(os, "id,depth,iterations_cnt");

    if (root) {
      iterations_recursive(os, *root, 0);
    }

    return os.str();
  }

  std::string strong_branching_iterations_to_csv() const {
    std::stringstream os;

    std::println(os, "id,depth,variable_id,iterations_cnt");

    if (root) {
      strong_branching_iterations_recursive(os, *root, 0);
    }

    return os.str();
  }

  std::string to_graphviz(bool minimalistic = false) const {
    std::stringstream os;

    std::println(os, "digraph G {{");
    std::println(os, "  node [shape=box];");

    if (root) {
      node_to_graphviz(os, *root, minimalistic);
    }

    std::println(os, "}}");
    return os.str();
  }

  std::string to_csv() const {
    std::stringstream os;

    std::println(os, "id,value,state,left_child,right_child");

    if (root) {
      node_to_csv(os, *root);
    }

    return os.str();
  }

  // methods, called by B & B algorithm
  void set_root(size_t id) {
    assert(!root);
    root = id;
    tree.emplace(id, Node{Unvisited{}});
  }

  void pruned_by_infeasibility(size_t id) {
    assert(std::holds_alternative<Unvisited>(tree.at(id).state));

    tree.at(id).state = PrunedByInfeasibility{};
  }
  void pruned_by_integrality(size_t id,
                             const FiniteLPSolution<Field>& solution) {
    assert(std::holds_alternative<Unvisited>(tree.at(id).state));

    tree.at(id).state = PrunedByIntegrality{solution.value};
  }
  void pruned_by_bounds(size_t id, const FiniteLPSolution<Field>& solution) {
    assert(std::holds_alternative<Unvisited>(tree.at(id).state));

    tree.at(id).state = PrunedByBounds{solution.value};
  }

  void set_branching_variable(size_t id,
                              const FiniteLPSolution<Field>& solution,
                              size_t branch_variable, Field branch_value) {
    assert(std::holds_alternative<Unvisited>(tree.at(id).state));

    tree.at(id).state = Branched{solution.value, branch_variable, branch_value};
  }

  void set_child(size_t parent_id, size_t child_id,
                 NodeRelativeLocation location) {
    assert(location != NodeRelativeLocation::ROOT);

    if (location == NodeRelativeLocation::LEFT_CHILD) {
      assert(!tree.at(parent_id).left_child);
      tree.at(parent_id).left_child = child_id;
    } else {
      assert(!tree.at(parent_id).right_child);
      tree.at(parent_id).right_child = child_id;
    }

    tree.emplace(child_id, Node{Unvisited{}});
  }

  void simplex_run(size_t node_id, const SimplexResult<Field>& result) {
    tree.at(node_id).simplex_iterations = result.iterations_count;
  }
  void strong_branching_simplex_run(size_t node_id, size_t variable_id,
                                    const SimplexResult<Field>& result) {
    tree.at(node_id).strong_branching_iterations.emplace_back(
        variable_id, result.iterations_count);
  }
};
