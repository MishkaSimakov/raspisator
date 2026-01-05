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
}  // namespace str
