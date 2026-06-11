#pragma once

#include <chrono>
#include <fstream>

#include "linear/model/LP.h"
#include "utils/String.h"

namespace simplex {

namespace detail {

static size_t get_dump_id() {
  return std::chrono::system_clock::now().time_since_epoch() /
         std::chrono::milliseconds(1);
}

}  // namespace detail

template <typename Field>
void dump_state(const problem::StandardLP<Field>& problem,
                const std::vector<VariableState>& var_states) {
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
    std::string bound = "std::pair{";

    if (problem.var_bounds[i].lower) {
      bound += std::format("{}", *problem.var_bounds[i].lower);
    } else {
      bound += "std::nullopt";
    }

    bound += ",";

    if (problem.var_bounds[i].upper) {
      bound += std::format("{}", *problem.var_bounds[i].upper);
    } else {
      bound += "std::nullopt";
    }

    bound += "}";

    string_bounds[i] = bound;
  }

  os << "Bounds<Field> bounds = {" << str::join(string_bounds, ", ") << "};\n";

  os << "std::vector<VariableState> last_states = {";
  for (auto var : var_states) {
    switch (var) {
      case VariableState::BASIC:
        os << "VariableState::BASIC, ";
        break;
      case VariableState::AT_LOWER:
        os << "VariableState::AT_LOWER, ";
        break;
      case VariableState::AT_UPPER:
        os << "VariableState::AT_UPPER, ";
        break;
      case VariableState::NONBASIC_FREE:
        os << "VariableState::NONBASIC_FREE, ";
        break;
      default:
        throw std::runtime_error("Unknown variable state.");
    }
  }
  os << "};\n";

  os << "};";

  os << "}\n";

  os << std::flush;

  std::println("Registered failed simplex run into {}.", dump_name);
}

}  // namespace simplex
