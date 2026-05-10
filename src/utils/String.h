#pragma once

#include <algorithm>
#include <cctype>
#include <locale>

namespace str {
// Trim from the start
inline std::string ltrim(std::string s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
          }));

  return s;
}

// Trim from the end
inline std::string rtrim(std::string s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char ch) { return !std::isspace(ch); })
              .base(),
          s.end());

  return s;
}

inline std::string trim(std::string s) { return ltrim(rtrim(s)); }

inline bool is_space(char c) {
  return std::isspace(static_cast<unsigned char>(c)) != 0;
}

inline bool all_spaces(std::string_view s) {
  return std::ranges::all_of(
      s, [](unsigned char c) { return std::isspace(c) != 0; });
}

std::string join(std::ranges::range auto&& range, std::string_view delimiter) {
  std::string result;

  auto itr = std::ranges::begin(range);
  auto end = std::ranges::end(range);

  if (itr == end) {
    return result;
  }

  result += *itr;
  ++itr;

  while (itr != end) {
    result += delimiter;
    result += *itr;

    ++itr;
  }

  return result;
}
}  // namespace str
