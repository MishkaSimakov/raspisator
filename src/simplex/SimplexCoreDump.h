#pragma once

#include <chrono>
#include <fstream>

#include "utils/String.h"

namespace simplex {

namespace detail {

static size_t get_dump_id() {
  return std::chrono::system_clock::now().time_since_epoch() /
         std::chrono::milliseconds(1);
}

static void dump_states(std::ostream& os, std::string_view name,
                        std::span<const VariableState> states) {
  os << "std::vector<simplex::VariableState> " << name << " = {"
     << str::join(states | std::views::transform([](VariableState state) {
                    return "simplex::VariableState::" + to_string(state);
                  }),
                  ", ")
     << "};\n";
}

}  // namespace detail

template <typename Field>
void dump_state(const problem::StandardLP<Field>& problem,
                const std::vector<VariableState>& init_states,
                const std::vector<VariableState>& last_states) {
  const size_t dump_id = detail::get_dump_id();
  std::string dump_name = std::format("simplex_core_dump_{}.h", dump_id);

  std::ofstream os(dump_name);

  if (!os) {
    throw std::runtime_error("Failed to open file for simplex core dump");
  }

  os << "namespace SimplexDump_" << dump_id << " {\n";

  os << "Matrix<Field> A = {" << problem.matrix << "};\n";
  os << "Matrix<Field> b = {" << problem.rhs << "};\n";
  os << "Matrix<Field> c = {" << problem.cost << "};\n";

  std::vector<std::string> string_bounds(problem.var_bounds.size());
  for (size_t i = 0; i < problem.var_bounds.size(); ++i) {
    std::string bound = "Bound<Field>(";

    if (problem.var_bounds[i].lower) {
      bound += std::format("{}", *problem.var_bounds[i].lower);
    } else {
      bound += "std::nullopt";
    }

    bound += ", ";

    if (problem.var_bounds[i].upper) {
      bound += std::format("{}", *problem.var_bounds[i].upper);
    } else {
      bound += "std::nullopt";
    }

    bound += ")";

    string_bounds[i] = bound;
  }

  os << "std::vector<Bound<Field>> bounds = {" << str::join(string_bounds, ", ")
     << "};\n";

  detail::dump_states(os, "init_states", init_states);
  detail::dump_states(os, "last_states", last_states);

  os << "}\n";

  os << std::flush;

  std::println("Registered failed simplex run into {}.", dump_name);
}

}  // namespace simplex
