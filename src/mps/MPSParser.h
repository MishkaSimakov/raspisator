#pragma once

#include <array>
#include <cassert>
#include <iostream>
#include <memory>

#include "DataRecordTokenizer.h"
#include "IndicatorRecordTokenizer.h"
#include "Types.h"
#include "utils/String.h"

#include "sections/BoundsParser.h"
#include "sections/ColumnsParser.h"
#include "sections/RHSParser.h"
#include "sections/RangesParser.h"
#include "sections/RowsParser.h"
#include "sections/SectionParser.h"

namespace mps {

template <typename Field>
class MPSParser {
  constexpr static size_t sections_count =
      static_cast<size_t>(SectionType::SECTIONS_COUNT);

  static bool should_skip_line(std::string_view line) {
    if (line.empty() || line[0] == '*' || line[0] == '$') {
      return true;
    }

    if (str::all_spaces(line)) {
      return true;
    }

    return false;
  }

  static auto init_parsers() {
    std::array<std::unique_ptr<SectionParser<Field>>, sections_count> parsers;

    parsers[static_cast<size_t>(SectionType::ROWS)] =
        std::make_unique<RowsParser<Field>>();
    parsers[static_cast<size_t>(SectionType::COLUMNS)] =
        std::make_unique<ColumnsParser<Field>>();
    parsers[static_cast<size_t>(SectionType::RHS)] =
        std::make_unique<RHSParser<Field>>();
    parsers[static_cast<size_t>(SectionType::BOUNDS)] =
        std::make_unique<BoundsParser<Field>>();
    parsers[static_cast<size_t>(SectionType::RANGES)] =
        std::make_unique<RangesParser<Field>>();

    return parsers;
  }

 public:
  static MPSParsingState<Field> parse(std::istream& is, Format format) {
    auto parsers = init_parsers();

    MPSParsingState<Field> state;

    std::string line;
    std::optional<SectionType> current_section = std::nullopt;

    while (std::getline(is, line)) {
      if (should_skip_line(line)) {
        continue;
      }

      assert(!line.empty() && "empty line should've been skipped");

      if (line[0] != ' ') {
        // indicator record
        const auto record = IndicatorRecordTokenizer::parse(line);

        if (record.type == SectionType::NAME) {
          state.problem_name = record.data;
          continue;
        }

        if (!record.data.empty()) {
          throw std::runtime_error(std::format(
              "Indicator record of type {} doesn't accept additional data.",
              section_type_to_string(record.type)));
        }

        if (record.type == SectionType::ENDATA) {
          break;
        }

        current_section = record.type;
      } else {
        // data record
        if (!current_section.has_value()) {
          throw std::runtime_error("Data record must be inside section.");
        }

        const auto record = DataRecordTokenizer::parse(
            line, format,
            parsers[static_cast<size_t>(*current_section)]->has_field_1());

        parsers[static_cast<size_t>(*current_section)]->parse(record, state);
      }
    }

    return state;
  }
};

}  // namespace mps
