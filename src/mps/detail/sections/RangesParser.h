#pragma once

#include <cassert>
#include <format>

#include "SectionParser.h"
#include "linear/FieldTraits.h"

namespace mps::detail {

template <typename Field>
class RangesParser final : public SectionParser<Field> {
  void parse_coefficient(const DataRecord& record,
                         MPSParsingState<Field>& state, size_t index) {
    assert(index == 0 || index == 1);

    const auto row = record.fields[2 + 2 * index];

    const auto itr = state.rows_map.find(row);
    if (itr == state.rows_map.end()) {
      throw std::runtime_error(std::format("Unknown row name: {}.", row));
    }

    if (str::all_spaces(record.fields[3 + 2 * index])) {
      throw std::runtime_error("Coefficient must be specified.");
    }

    if (state.rows[itr->second].range.has_value()) {
      throw std::runtime_error(
          std::format("Range for row {} is specified twice.",
                      state.rows[itr->second].name));
    }

    state.rows[itr->second].range =
        FieldTraits<Field>::from_string(record.fields[3 + 2 * index]);
  }

 public:
  bool has_field_1() const override { return false; }

  void parse(const DataRecord& record, MPSParsingState<Field>& state) override {
    const auto vector_name = record.fields[1];

    if (!state.ranges_vector_name) {
      state.ranges_vector_name = vector_name;
    }

    if (*state.ranges_vector_name != vector_name) {
      return;
    }

    if (str::all_spaces(record.fields[2])) {
      throw std::runtime_error(
          "There should be at least one row specified for RANGES data record.");
    }

    parse_coefficient(record, state, 0);

    if (!str::all_spaces(record.fields[4])) {
      parse_coefficient(record, state, 1);
    }
  }
};

}  // namespace mps::detail
