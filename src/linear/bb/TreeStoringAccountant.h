#pragma once

#include <cassert>
#include <optional>
#include <sstream>
#include <unordered_map>
#include <variant>

#include "linear/bb/Accountant.h"

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
    std::optional<size_t> left_child;
    std::optional<size_t> right_child;

    std::variant<PrunedByInfeasibility, PrunedByIntegrality, PrunedByBounds,
                 Branched, Unvisited>
        state;
  };

  std::unordered_map<size_t, Node> tree;
  std::optional<size_t> root;

  void emplace_unvisited(size_t id) {
    tree.emplace(id, Node{std::nullopt, std::nullopt, Unvisited{}});
  }

  // graphviz methods
 private:
  static std::string get_node_graphviz_attributes(PrunedByInfeasibility) {
    return "label=\"unfeasible\", style=filled, fillcolor=red";
  }

  static std::string get_node_graphviz_attributes(
      const PrunedByIntegrality& state) {
    return std::format(R"(label="integer\n{}", style=filled, fillcolor=green)",
                       state.value);
  }

  static std::string get_node_graphviz_attributes(const PrunedByBounds& state) {
    return std::format(R"(label="bounded\n{}", style=filled, fillcolor=orange)",
                       state.value);
  }

  static std::string get_node_graphviz_attributes(const Branched& state) {
    return std::format(R"(label="{}\nx_{} <> {}")", state.value,
                       state.branch_variable, state.branch_value);
  }

  static std::string get_node_graphviz_attributes(Unvisited) {
    return "label=\"unvisited\", style=filled, fillcolor=gray";
  }

  void node_to_graphviz(std::stringstream& os, size_t id) const {
    const Node& node = tree.at(id);

    std::string attrs = std::visit(
        [](const auto& node) { return get_node_graphviz_attributes(node); },
        node.state);
    std::println(os, "  {} [{}]", id, attrs);

    if (node.left_child) {
      std::println(os, "  {} -> {} [label=\"left\"];", id, *node.left_child);
      node_to_graphviz(os, *node.left_child);
    }

    if (node.right_child) {
      std::println(os, "  {} -> {} [label=\"right\"];", id, *node.right_child);
      node_to_graphviz(os, *node.right_child);
    }
  }

 public:
  std::string to_graphviz() const {
    std::stringstream os;

    std::println(os, "digraph G {{");
    std::println(os, "  node [shape=box];");

    if (root) {
      node_to_graphviz(os, *root);
    }

    std::println(os, "}}");
    return os.str();
  }

  // methods, called by B & B algorithm
  void set_root(size_t id) {
    assert(!root);
    root = id;
    emplace_unvisited(id);
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

  void set_left_child(size_t parent_id, size_t child_id) {
    assert(!tree.at(parent_id).left_child);

    tree.at(parent_id).left_child = child_id;
    emplace_unvisited(child_id);
  }
  void set_right_child(size_t parent_id, size_t child_id) {
    assert(!tree.at(parent_id).right_child);

    tree.at(parent_id).right_child = child_id;
    emplace_unvisited(child_id);
  }
};
