#pragma once

#include "BranchAndBoundTree.h"

// draws branch and bound tree using graphviz
template <typename Field>
class GraphvizBuilder {
 public:
  std::string build(const BranchAndBoundTree<Field>& tree) {
    ss_ << "digraph G {\n";
    ss_ << "  node [shape=box];\n";

    if (tree.get_root() != nullptr) {
      add_node(tree.get_root());
    }

    ss_ << "}\n";
    return ss_.str();
  }

 private:
  std::stringstream ss_;
  std::unordered_map<const void*, std::string> ids_;

  std::string get_id(const void* ptr) {
    auto it = ids_.find(ptr);
    if (it != ids_.end()) {
      return it->second;
    }

    std::string id = "n" + std::to_string(ids_.size());
    ids_[ptr] = id;
    return id;
  }

  void add_node(const Node<Field>* node) {
    std::visit([this](const auto& node) { add_node_specialized(node); }, *node);
  }

  void add_node_specialized(const InteriorNode<Field>& node) {
    std::string id = get_id(&node);

    ss_ << std::format("  {} [label=\"{}\\nx_{} <> {}\"];\n", id, node.value,
                       node.branching_variable, node.branching_value);

    if (node.left) {
      ss_ << std::format("  {} -> {} [label=\"left\"];\n", id,
                         get_id(node.left));
      add_node(node.left);
    }

    if (node.right) {
      ss_ << std::format("  {} -> {} [label=\"right\"];\n", id,
                         get_id(node.right));
      add_node(node.right);
    }
  }

  void add_node_specialized(const IntegerSolutionNode<Field>& node) {
    ss_ << std::format(
        "  {} [label=\"integer\\n{}\", style=filled, fillcolor=green];\n",
        get_id(&node), node.value);
  }

  void add_node_specialized(const BoundedNode<Field>& node) {
    ss_ << std::format(
        "  {} [label=\"bounded\\n{}\", style=filled, fillcolor=orange];\n",
        get_id(&node), node.value);
  }

  void add_node_specialized(const UnfeasibleNode<Field>& node) {
    ss_ << std::format(
        "  {} [label=\"unfeasible\", style=filled, fillcolor=red];\n",
        get_id(&node));
  }
};
