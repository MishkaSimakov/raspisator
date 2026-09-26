#pragma once

#include <format>
#include <sstream>
#include <utility>

template <typename T1, typename T2>
struct std::formatter<std::pair<T1, T2>, char> : std::formatter<std::string> {
  template <class FmtContext>
  auto format(const std::pair<T1, T2>& value, FmtContext& ctx) const {
    std::ostringstream out;
    out << "(" << value.first << ", " << value.second << ")";

    return std::ranges::copy(std::move(out).str(), ctx.out()).out;
  }
};
