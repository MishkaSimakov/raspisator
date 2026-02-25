#pragma once

#include <iostream>
#include <optional>
#include <print>
#include <string>
#include <vector>

struct NodeParameters {
  size_t id;
  std::string label;
  std::optional<std::string> fill;
};

struct EdgeParameters {
  size_t from;
  size_t to;

  std::string label;
};

struct DrawerSettings {
  bool draw_labels = true;
};

class GraphvizDrawer {
  DrawerSettings settings_;

  std::vector<NodeParameters> nodes_;
  std::vector<EdgeParameters> edges_;

  void draw_nodes(std::ostream& os) const {
    for (const auto& node : nodes_) {
      std::string attributes;

      if (settings_.draw_labels && !node.label.empty()) {
        attributes += std::format("label=\"{}\",", node.label);
      }

      if (node.fill) {
        attributes += std::format("style=filled,fillcolor={},", *node.fill);
      }

      std::println(os, "  {} [{}]", node.id, attributes);
    }
  }

  void draw_edges(std::ostream& os) const {
    for (const auto& edge : edges_) {
      std::string attributes;

      if (settings_.draw_labels && !edge.label.empty()) {
        attributes += std::format("label=\"{}\",", edge.label);
      }

      std::println(os, "  {} -> {} [{}];", edge.from, edge.to, attributes);
    }
  }

 public:
  GraphvizDrawer() = default;

  void add_node(NodeParameters parameters) {
    nodes_.push_back(std::move(parameters));
  }

  void add_edge(EdgeParameters parameters) {
    edges_.push_back(std::move(parameters));
  }

  void draw(std::ostream& os) const {
    std::println(os, "digraph G {{");
    std::println(os, "  node [shape=box];");

    draw_nodes(os);

    draw_edges(os);

    std::println(os, "}}");
  }
};
