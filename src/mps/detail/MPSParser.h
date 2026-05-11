#pragma once

#include <array>
#include <cassert>
#include <iostream>
#include <memory>

#include "DataRecordTokenizer.h"
#include "IndicatorRecordTokenizer.h"
#include "Types.h"
#include "mps/Format.h"
#include "utils/String.h"

#include "sections/BoundsParser.h"
#include "sections/ColumnsParser.h"
#include "sections/RHSParser.h"
#include "sections/RangesParser.h"
#include "sections/RowsParser.h"
#include "sections/SectionParser.h"

namespace mps::detail {

template <typename Field>
class MPSParser {
  constexpr static size_t sections_count =
      static_cast<size_t>(SectionType::SECTIONS_COUNT);

  constexpr static SectionType mandatory_sections[] = {
      SectionType::NAME, SectionType::ROWS,   SectionType::COLUMNS,
      SectionType::RHS,  SectionType::ENDATA,
  };

  struct Section {
    std::unique_ptr<SectionParser<Field>> parser{nullptr};
    bool visited{false};

    Section() = default;

    explicit Section(std::unique_ptr<SectionParser<Field>> parser)
        : parser(std::move(parser)) {}
  };

  static bool should_skip_line(std::string_view line) {
    if (line.empty() || line[0] == '*' || line[0] == '$') {
      return true;
    }

    if (str::all_spaces(line)) {
      return true;
    }

    return false;
  }

  static std::array<Section, sections_count> init_parsers() {
    std::array<Section, sections_count> parsers;

    parsers[static_cast<size_t>(SectionType::ROWS)] =
        Section(std::make_unique<RowsParser<Field>>());
    parsers[static_cast<size_t>(SectionType::COLUMNS)] =
        Section(std::make_unique<ColumnsParser<Field>>());
    parsers[static_cast<size_t>(SectionType::RHS)] =
        Section(std::make_unique<RHSParser<Field>>());
    parsers[static_cast<size_t>(SectionType::BOUNDS)] =
        Section(std::make_unique<BoundsParser<Field>>());
    parsers[static_cast<size_t>(SectionType::RANGES)] =
        Section(std::make_unique<RangesParser<Field>>());

    return parsers;
  }

  static void check_mandatory_sections(
      const MPSParsingState<Field>& state,
      const std::array<Section, sections_count>& sections) {
    for (const SectionType section : mandatory_sections) {
      if (!sections[static_cast<size_t>(section)].visited) {
        throw ParseError(std::format("Section '{}' is mandatory.",
                                     section_type_to_string(section)));
      }
    }
  }

  static bool section_has_data(SectionType type) {
    return type == SectionType::NAME || type == SectionType::OBJECT;
  }

 public:
  static MPSParsingState<Field> parse(std::istream& is, Format format) {
    auto sections = init_parsers();

    MPSParsingState<Field> state;

    size_t row_index = 0;
    std::string line;
    std::optional<SectionType> current_section = std::nullopt;

    while (std::getline(is, line)) {
      if (format == Format::FIXED) {
        line.resize(std::max(line.size(), 61uz), ' ');
      }

      ++row_index;

      if (should_skip_line(line)) {
        continue;
      }

      assert(!line.empty() && "empty line should've been skipped");

      try {
        if (line[0] != ' ') {
          // indicator record
          const auto record = IndicatorRecordTokenizer::parse(line);

          // NAME section may be duplicated
          if (record.type == SectionType::NAME) {
            sections[static_cast<size_t>(record.type)].visited = true;
            state.problem_name = record.data;
            continue;
          }

          // check if we visited this type of section before
          if (sections[static_cast<size_t>(record.type)].visited) {
            throw ParseError(std::format("Section '{}' is duplicated.",
                                         section_type_to_string(record.type)));
          }
          sections[static_cast<size_t>(record.type)].visited = true;

          // teardown parser for previous section
          if (current_section != std::nullopt) {
            auto& section = sections[static_cast<size_t>(*current_section)];

            if (section.parser != nullptr) {
              section.parser->teardown();
            }
          }

          if (!section_has_data(record.type) && !record.data.empty()) {
            throw ParseError(std::format(
                "Indicator record of type '{}' doesn't accept additional data.",
                section_type_to_string(record.type)));
          }

          if (record.type == SectionType::ENDATA) {
            break;
          }

          current_section = record.type;
        } else {
          // data record
          if (!current_section.has_value()) {
            throw ParseError("Data record must be inside section.");
          }

          if (*current_section == SectionType::OBJECT) {
            throw ParseError("OBJECT section must not contain data records.");
          }

          auto& section = sections[static_cast<size_t>(*current_section)];

          assert(section.parser != nullptr);

          const auto record = DataRecordTokenizer::parse(
              line, format, section.parser->has_field_1());

          section.parser->parse(record, state);
        }
      } catch (ParseError& error) {
        error.set_line(row_index);
        throw;
      }
    }

    try {
      // teardown parser for the last section
      if (current_section != std::nullopt) {
        auto& section = sections[static_cast<size_t>(*current_section)];

        if (section.parser != nullptr) {
          section.parser->teardown();
        }
      }

      check_mandatory_sections(state, sections);
    } catch (ParseError& error) {
      error.set_line(row_index);
      throw;
    }

    return state;
  }
};

}  // namespace mps::detail
