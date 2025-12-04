#pragma once

#include <algorithm>
#include <deque>
#include <variant>

#include "linear/matrix/Matrix.h"

// Branch and Bound nodes
template <typename Field>
struct Node;

template <typename Field>
struct InteriorNode {
  size_t branching_variable;
  Field branching_value;

  InteriorNode* parent = nullptr;

  Node<Field>* left = nullptr;
  Node<Field>* right = nullptr;

  Field value;
  Matrix<Field> solution;
  std::vector<size_t> basic_variables;
};

template <typename Field>
struct UnfeasibleNode {
  InteriorNode<Field>* parent = nullptr;
};

template <typename Field>
struct IntegerSolutionNode {
  Field value;

  InteriorNode<Field>* parent = nullptr;
};

template <typename Field>
struct BoundedNode {
  Field value;

  InteriorNode<Field>* parent = nullptr;
};

template <typename Field>
struct Node : std::variant<InteriorNode<Field>, UnfeasibleNode<Field>,
                           IntegerSolutionNode<Field>, BoundedNode<Field>> {
  using std::variant<InteriorNode<Field>, UnfeasibleNode<Field>,
                     IntegerSolutionNode<Field>, BoundedNode<Field>>::variant;

  void set_parent(InteriorNode<Field>* parent) {
    std::visit([parent](auto& node) { node.parent = parent; }, *this);
  }
};

enum class NodeRelativeLocation { LEFT_CHILD, RIGHT_CHILD, ROOT };

template <typename Field>
class BranchAndBoundTree {
  std::deque<Node<Field>> nodes_;

 public:
  // returns path from given node to the root and directions along this path
  static std::pair<std::vector<const InteriorNode<Field>*>,
                   std::vector<NodeRelativeLocation>>
  get_path_to(const InteriorNode<Field>* node) {
    std::vector<const InteriorNode<Field>*> nodes;
    std::vector<NodeRelativeLocation> directions;

    while (true) {
      nodes.push_back(node);
      auto* parent = node->parent;

      if (parent == nullptr) {
        break;
      }

      directions.push_back(std::get_if<InteriorNode<Field>>(parent->left) ==
                                   node
                               ? NodeRelativeLocation::LEFT_CHILD
                               : NodeRelativeLocation::RIGHT_CHILD);

      node = parent;
    }

    return {std::move(nodes), std::move(directions)};
  }

  template <typename T>
  T& add_node(InteriorNode<Field>* parent, NodeRelativeLocation location,
              T value) {
    if (parent == nullptr && location != NodeRelativeLocation::ROOT) {
      throw std::invalid_argument(
          "Only root node is allowed to be without parent.");
    }

    if (location == NodeRelativeLocation::ROOT && !nodes_.empty()) {
      throw std::invalid_argument("There must be only one root.");
    }

    auto& node = nodes_.emplace_back(value);

    if (location != NodeRelativeLocation::ROOT) {
      node.set_parent(parent);

      if (location == NodeRelativeLocation::LEFT_CHILD) {
        parent->left = &node;
      } else {
        parent->right = &node;
      }
    }

    return std::get<T>(node);
  }

  const Node<Field>* get_root() const {
    if (nodes_.empty()) {
      return nullptr;
    }

    return &nodes_[0];
  }

  size_t size() const { return nodes_.size(); }
};
