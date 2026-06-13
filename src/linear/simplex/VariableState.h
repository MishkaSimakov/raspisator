#pragma once

#include <string>

namespace simplex {

enum class VariableState { AT_LOWER, AT_UPPER, NONBASIC_FREE, BASIC };

inline std::string to_string(VariableState state) {
  switch (state) {
    case VariableState::AT_LOWER:
      return "AT_LOWER";
    case VariableState::AT_UPPER:
      return "AT_UPPER";
    case VariableState::NONBASIC_FREE:
      return "NONBASIC_FREE";
    case VariableState::BASIC:
      return "BASIC";
    default:
      throw std::runtime_error("Unknown variable state.");
  }
}

}  // namespace simplex
