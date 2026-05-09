#pragma once

#include "DataRecordTokenizer.h"
#include "Types.h"

#include "linear/problem/MILPProblem.h"
#include "sections/RowsParser.h"
#include "sections/SectionParser.h"
#include "utils/String.h"

namespace mps {

template <typename Field>
class MPSReader {
  constexpr static std::array indicator_record_types = {
      std::pair{"NAME", SectionType::NAME},
      std::pair{"ROWS", SectionType::ROWS},
      std::pair{"COLUMNS", SectionType::COLUMNS},
      std::pair{"RHS", SectionType::RHS},
      std::pair{"BOUNDS", SectionType::BOUNDS},
      std::pair{"RANGES", SectionType::RANGES},
      std::pair{"ENDATA", SectionType::ENDATA},
  };

  constexpr static size_t sections_count =
      static_cast<size_t>(SectionType::SECTIONS_COUNT);

  std::array<std::unique_ptr<SectionParser<Field>>, sections_count> parsers_;

  MPSReader() {
    parsers_[SectionType::ROWS] = std::make_unique<RowsParser<Field>>();
  }

  bool should_skip_line(std::string_view line) {
    if (line.empty() || line[0] == '*' || line[0] == '$') {
      return true;
    }

    if (str::all_spaces(line)) {
      return true;
    }

    return false;
  }

  static SectionType parse_indicator_record(std::string_view line) {
    for (const auto& [name, value] : indicator_record_types) {
      if (line == name) {
        return value;
      }
    }

    throw std::runtime_error("Unknown indicator record.");
  }

  void parse_mps(std::istream& is, Format format) {
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
        current_section = parse_indicator_record(line);

        if (current_section == SectionType::ENDATA) {
          break;
        }
      } else {
        // data record
        if (!current_section.has_value()) {
          throw std::runtime_error("Data record must be inside section.");
        }

        const auto record = DataRecordTokenizer::parse(
            line, format, parsers_[current_section]->has_field_1());

        parsers_[current_section]->parse(record, state);
      }
    }
  }

  MILPProblem<Field> generate_problem() const {}

 public:
  static MILPProblem<Field> read(std::istream& is, Format format) {
    auto reader = MPSReader();

    reader.parse_mps(is, format);

    return reader.generate_problem();
  }
};

}  // namespace mps
