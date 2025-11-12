#include "utils/Drawing.h"

#include <fmt/format.h>
#include <fmt/ranges.h>

std::vector<std::string> translate_states(const STN& stn) {
  const auto& states = stn.get_states();

  std::vector<std::string> result(states.size());

  for (size_t i = 0; i < states.size(); ++i) {
    size_t id = states[i].get_id();
    result[i] = fmt::format("S{} [label=\"State {}\n(-,-,-)\"];", id, id);
  }

  return result;
}

std::vector<std::string> translate_tasks(const STN& stn) {
  const auto& tasks = stn.get_tasks();

  std::vector<std::string> result(tasks.size());

  for (size_t i = 0; i < tasks.size(); ++i) {
    size_t id = tasks[i].get_id();
    std::vector<std::string> units_description;

    for (const auto [unit, props] : stn.get_task_units(tasks[i])) {
      units_description.push_back(
          fmt::format("{{ Unit {}|({}u, bs: {}-{}) }}", unit->get_id(),
                      props.batch_processing_time, props.batch_min_size,
                      props.batch_max_size));
    }

    result[i] = fmt::format("T{} [label=\"{{ Task {} | {} }}\"];", id, id,
                            fmt::join(units_description, "|"));
  }

  return result;
}

std::vector<std::string> translate_arcs(const STN& stn) {
  const auto& tasks = stn.get_tasks();

  std::vector<std::string> result;

  for (size_t i = 0; i < tasks.size(); ++i) {
    size_t id = tasks[i].get_id();

    for (const State* input : tasks[i].get_inputs()) {
      result.push_back(fmt::format("S{} -> T{}", input->get_id(), id));
    }

    for (const State* output : tasks[i].get_outputs()) {
      result.push_back(fmt::format("T{} -> S{}", id, output->get_id()));
    }
  }

  return result;
}

std::string to_graphviz(const STN& stn) {
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

  return fmt::format(graph_template, fmt::join(states, "\n"),
                     fmt::join(tasks, "\n"), fmt::join(arcs, "\n"));
}
