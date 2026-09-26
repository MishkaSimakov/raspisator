#pragma once

#include <cassert>
#include <format>

#include "../../../field/FieldTraits.h"
#include "SectionParser.h"

namespace mps::detail {

template <typename Field>
class ColumnsParser final : public SectionParser<Field> {
  enum class MarkerType {
    INTEGER_BEGIN,
    INTEGER_END,
  };

  constexpr static std::string_view marker_string = "'MARKER'";

  bool is_inside_integer_ = false;

  // after begin or end of integer section, column must be changed
  bool must_change_column_ = false;

  void parse_coefficient(const DataRecord& record,
                         MPSParsingState<Field>& state, size_t index) {
    assert(index == 0 || index == 1);

    const auto row = record.fields[2 + 2 * index];

    const auto itr = state.rows_map.find(row);
    if (itr == state.rows_map.end()) {
      throw ParseError(std::format("Unknown row name: '{}'.", row));
    }

    if (str::all_spaces(record.fields[3 + 2 * index])) {
      throw ParseError("Coefficient must be specified.");
    }

    const std::optional<Field> parsed_value = FieldTraits<Field>::from_string(
        str::trim(record.fields[3 + 2 * index]));

    if (!parsed_value) {
      throw ParseError(
          std::format("Failed to parse value in Field {}", 4 + 2 * index));
    }

    auto [_, inserted] =
        state.cols.back().values.emplace(itr->second, *parsed_value);

    if (!inserted) {
      throw ParseError(std::format("Row '{}' is duplicated in column '{}'.",
                                   state.rows[itr->second].name,
                                   state.cols.back().name));
    }
  }

  MarkerType parse_marker_type(std::string_view type) {
    if (type == "'INTORG'") {
      return MarkerType::INTEGER_BEGIN;
    }
    if (type == "'INTEND'") {
      return MarkerType::INTEGER_END;
    }

    throw ParseError(std::format("Unknown marker type: '{}'.", type));
  }

  void parse_marker(const DataRecord& record, MPSParsingState<Field>& state) {
    const auto name = record.fields[1];

    if (state.cols_map.contains(name)) {
      throw ParseError(
          std::format("Marker name must differ from the preceding and "
                      "succeeding column names. Name '{}' is duplicated.",
                      name));
    }

    state.markers.insert(std::string(name));

    const auto type = parse_marker_type(record.fields[3]);

    switch (type) {
      case MarkerType::INTEGER_BEGIN:
        if (is_inside_integer_) {
          throw ParseError(
              "Integer section begin marker while the previous section has not "
              "ended yet.");
        }

        is_inside_integer_ = true;
        must_change_column_ = true;
        break;
      case MarkerType::INTEGER_END:
        if (!is_inside_integer_) {
          throw ParseError(
              "Integer section end marker while not in integer section.");
        }

        is_inside_integer_ = false;
        must_change_column_ = true;
        break;
      default:
        throw ParseError("Unknown marker type.");
    }
  }

  void parse_column(const DataRecord& record, MPSParsingState<Field>& state) {
    const auto name = record.fields[1];

    if (state.markers.contains(name)) {
      throw ParseError(
          std::format("Marker name must differ from the preceding and "
                      "succeeding column names. Name '{}' is duplicated.",
                      name));
    }

    if (must_change_column_ && !state.cols.empty() &&
        state.cols.back().name == name) {
      throw ParseError(
          std::format("Data records for column '{}' exist both inside and "
                      "outside of integer section.",
                      name));
    }

    if (state.cols.empty() || state.cols.back().name != name) {
      // new column name
      const bool inserted = state.add_col(name);
      if (!inserted) {
        throw ParseError(
            std::format("Non-consecutive entries for column '{}'.", name));
      }

      must_change_column_ = false;
    }

    if (is_inside_integer_) {
      state.cols.back().is_integer = true;
    }

    if (str::all_spaces(record.fields[2])) {
      throw ParseError(
          "There should be at least one row specified for column data record.");
    }

    parse_coefficient(record, state, 0);

    if (!str::all_spaces(record.fields[4])) {
      parse_coefficient(record, state, 1);
    }
  }

 public:
  bool has_field_1() const override { return false; }

  void parse(const DataRecord& record, MPSParsingState<Field>& state) override {
    if (record.fields[2] == marker_string) {
      parse_marker(record, state);
    } else {
      parse_column(record, state);
    }
  }

  void teardown() override {
    if (is_inside_integer_) {
      throw ParseError("Integer section is not closed.");
    }
  }
};

}  // namespace mps::detail
