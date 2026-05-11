#pragma once

#include <array>
#include <cassert>
#include <string_view>

#include "Types.h"
#include "mps/ParseError.h"
#include "utils/String.h"

namespace mps::detail {

struct IndicatorRecord {
  SectionType type;

  // Some sections contain additional data after their name. For example, NAME
  // section.
  // This field points at a range from the first non-space symbol after section
  // name to the end of the string.
  std::string_view data;
};

class IndicatorRecordTokenizer {
  constexpr static std::array indicator_record_types = {
      std::pair{"NAME", SectionType::NAME},
      std::pair{"ROWS", SectionType::ROWS},
      std::pair{"COLUMNS", SectionType::COLUMNS},
      std::pair{"RHS", SectionType::RHS},
      std::pair{"BOUNDS", SectionType::BOUNDS},
      std::pair{"RANGES", SectionType::RANGES},
      std::pair{"ENDATA", SectionType::ENDATA},
      std::pair{"OBJECT", SectionType::OBJECT},
  };

  static std::pair<SectionType, std::string_view> parse_type(
      std::string_view line) {
    size_t current = 0;
    while (current < line.size() && !str::is_space(line[current])) {
      ++current;
    }

    std::string_view type_string = line.substr(0, current);
    std::string_view remaining = line.substr(current);

    for (const auto& [name, value] : indicator_record_types) {
      if (type_string == name) {
        return {value, remaining};
      }
    }

    throw ParseError(
        std::format("Unknown indicator record type: '{}'.", type_string));
  }

 public:
  // @record must point at non-empty string that starts with non-space symbol
  static IndicatorRecord parse(std::string_view line) {
    assert(!line.empty() && !str::is_space(line[0]));

    auto [type, remaining] = parse_type(line);

    // consume leading spaces
    while (!remaining.empty() && str::is_space(remaining.front())) {
      remaining.remove_prefix(1);
    }

    return IndicatorRecord{
        .type = type,
        .data = remaining,
    };
  }
};

}  // namespace mps::detail
