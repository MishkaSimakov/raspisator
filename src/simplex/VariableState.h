#pragma once

#include <algorithm>
#include <ranges>
#include <string>

#include "utils/String.h"

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
      std::unreachable();
  }
}

inline std::string to_string(std::span<const VariableState> states) {
  return str::join(states | std::views::transform([](VariableState state) {
                     return to_string(state);
                   }),
                   " ");
}

}  // namespace simplex
