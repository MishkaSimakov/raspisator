#pragma once

#include <array>
#include <format>
#include <string>
#include <string_view>

#include "SectionParser.h"
#include "utils/String.h"

namespace mps {

template <typename Field>
class RowsParser final : public SectionParser<Field> {
  static RowSense parse_row_sense(std::string_view sense) {
    if (sense.size() != 1) {
      throw std::runtime_error("Unknown row sense.");
    }

    switch (sense[0]) {
      case 'E':
        return RowSense::EQUAL;
      case 'L':
        return RowSense::LESS_THAN;
      case 'G':
        return RowSense::GREATER_THAN;
      case 'N':
        return RowSense::FREE;
      default:
        throw std::runtime_error("Unknown row sense.");
    }
  }

 public:
  bool has_field_1() const override { return true; }

  void parse(DataRecord record, MPSParsingState<Field>& state) override {
    const size_t index = state.rows.size();
    const auto sense = parse_row_sense(record.fields[0]);
    const auto name = record.fields[1];

    auto [itr, inserted] = state.rows_map.emplace(name, index);
    if (!inserted) {
      throw std::runtime_error(std::format("Duplicated row name: {}.", name));
    }

    state.rows.emplace_back(sense, name);

    for (size_t i = 2; i < 6; ++i) {
      if (!str::all_spaces(record.fields[i])) {
        throw std::runtime_error("Fields 3-6 must be empty in ROWS section.");
      }
    }
  }
};

}  // namespace mps
