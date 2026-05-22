#pragma once

#include <cassert>
#include <format>

#include "SectionParser.h"
#include "linear/FieldTraits.h"
#include "mps/ParseError.h"

namespace mps::detail {

template <typename Field>
class RHSParser final : public SectionParser<Field> {
  void parse_coefficient(const DataRecord& record,
                         MPSParsingState<Field>& state, size_t index) {
    assert(index == 0 || index == 1);

    const auto row = record.fields[2 + 2 * index];

    const auto itr = state.rows_map.find(row);
    if (itr == state.rows_map.end()) {
      throw ParseError(std::format("Unknown row name: \"{}\".", row));
    }

    if (str::all_spaces(record.fields[3 + 2 * index])) {
      throw ParseError("Coefficient must be specified.");
    }

    if (state.rows[itr->second].rhs.has_value()) {
      throw ParseError(std::format("RHS for row '{}' is specified twice.",
                                   state.rows[itr->second].name));
    }

    const std::optional<Field> parsed_value = FieldTraits<Field>::from_string(
        str::trim(record.fields[3 + 2 * index]));

    if (!parsed_value) {
      throw ParseError(
          std::format("Failed to parse value in Field {}.", 4 + 2 * index));
    }

    state.rows[itr->second].rhs = *parsed_value;
  }

 public:
  bool has_field_1() const override { return false; }

  void parse(const DataRecord& record, MPSParsingState<Field>& state) override {
    const auto vector_name = record.fields[1];

    if (!state.rhs_vector_name) {
      state.rhs_vector_name = vector_name;
    }

    if (*state.rhs_vector_name != vector_name) {
      return;
    }

    if (str::all_spaces(record.fields[2])) {
      throw ParseError(
          "At least one row must be specified for RHS data record.");
    }

    parse_coefficient(record, state, 0);

    if (!str::all_spaces(record.fields[4])) {
      parse_coefficient(record, state, 1);
    }
  }
};

}  // namespace mps::detail
