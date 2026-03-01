#pragma once

#include <format>
#include <string>

#include "String.h"
#include "model/STN.h"
#include "utils/Variant.h"

using namespace std::string_literals;

template <typename Field>
std::vector<std::string> translate_states(const STN<Field>& stn) {
  const auto& states = stn.get_states();

  std::vector<std::string> result(states.size());

  for (size_t i = 0; i < states.size(); ++i) {
    size_t id = states[i].get_id();

    auto props = std::visit(
        Overload{[](InputState<Field>) { return "in"s; },
                 [](OutputState<Field>) { return "out"s; },
                 [](NonStorableState<Field>) { return "ns"s; },
                 [](const NormalState<Field>& state) {
                   return std::format("({}, {}-{})", state.initial_stock,
                                      state.min_level, state.max_level);
                 }},
        states[i]);

    result[i] = std::format("S{} [label=\"State {}\n{}\"];", id, id, props);
  }

  return result;
}

template <typename Field>
std::vector<std::string> translate_tasks(const STN<Field>& stn) {
  const auto& tasks = stn.get_tasks();

  std::vector<std::string> result(tasks.size());

  for (size_t i = 0; i < tasks.size(); ++i) {
    size_t id = tasks[i].get_id();
    std::vector<std::string> units_description_parts;

    for (const auto [unit, props] : stn.get_task_units(tasks[i])) {
      units_description_parts.push_back(
          std::format("{{ Unit {}|({}u, bs: {}-{}) }}", unit->get_id(),
                      props.batch_processing_time, props.batch_min_size,
                      props.batch_max_size));
    }

    std::string units_description = str::join(units_description_parts, "|");

    result[i] = std::format("T{} [label=\"{{ Task {} | {} }}\"];", id, id,
                            units_description);
  }

  return result;
}

template <typename Field>
std::vector<std::string> translate_arcs(const STN<Field>& stn) {
  const auto& tasks = stn.get_tasks();

  std::vector<std::string> result;

  for (size_t i = 0; i < tasks.size(); ++i) {
    size_t id = tasks[i].get_id();

    for (const auto& [input, fraction] : tasks[i].get_inputs()) {
      result.push_back(std::format("S{} -> T{} [label=\"{}\"];",
                                   input->get_id(), id, fraction));
    }

    for (const auto& [output, fraction] : tasks[i].get_outputs()) {
      result.push_back(std::format("T{} -> S{} [label=\"{}\"];", id,
                                   output->get_id(), fraction));
    }
  }

  return result;
}

template <typename Field>
std::string to_graphviz(const STN<Field>& stn) {
  constexpr auto graph_template =
      "digraph process_flow {{\n"
      "rankdir=LR;\n"
      "fontname=\"Helvetica\";\n"
      "node [fontname=\"Helvetica\"];\n\n"
      "// --- States ---\n"
      "node [shape=ellipse, style=filled, fillcolor=white];\n"
      "{}\n\n"
      "// --- Tasks ---\n"
      "node [shape=record, style=filled, fillcolor=\"#f2f2f2\"];\n"
      "{}\n\n"
      "// --- Arcs ---\n"
      "{}\n"
      "}}\n";

  std::vector<std::string> states = translate_states(stn);
  std::vector<std::string> tasks = translate_tasks(stn);
  std::vector<std::string> arcs = translate_arcs(stn);

  return std::format(graph_template, str::join(states, "\n"),
                     str::join(tasks, "\n"), str::join(arcs, "\n"));
}
