#pragma once

#include <cassert>
#include <format>

#include "SectionParser.h"
#include "linear/FieldTraits.h"

namespace mps {

template <typename Field>
class ColumnsParser final : public SectionParser<Field> {
  enum class MarkerType {
    INTEGER_BEGIN,
    INTEGER_END,
  };

  constexpr static std::string marker_string = "'MARKER'";

  bool is_inside_integer_ = false;

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

    auto [_, inserted] = state.cols.back().values.emplace(
        itr->second,
        FieldTraits<Field>::from_string(record.fields[3 + 2 * index]));

    if (!inserted) {
      throw std::runtime_error(std::format("Row {} is duplicated in column {}.",
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

    throw std::runtime_error(std::format("Unknown marker type: {}.", type));
  }

  void parse_marker(const DataRecord& record) {
    const auto name = record.fields[1];
    const auto type = parse_marker_type(record.fields[3]);

    switch (type) {
      case MarkerType::INTEGER_BEGIN:
        if (is_inside_integer_) {
          throw std::runtime_error(
              "Integer section begin marker while the previous section has not "
              "ended yet.");
        }

        is_inside_integer_ = true;
        break;
      case MarkerType::INTEGER_END:
        if (!is_inside_integer_) {
          throw std::runtime_error(
              "Integer section end marker while not in integer section.");
        }

        is_inside_integer_ = false;
        break;
      default:
        throw std::runtime_error("Unknown marker type.");
    }
  }

 public:
  bool has_field_1() const override { return false; }

  void parse(const DataRecord& record, MPSParsingState<Field>& state) override {
    if (record.fields[2] == marker_string) {
      parse_marker(record);
      return;
    }

    const auto column_name = record.fields[1];

    if (state.cols.empty() || state.cols.back().name != column_name) {
      // new column name
      const bool inserted = state.add_col(column_name);
      if (!inserted) {
        throw std::runtime_error(
            std::format("Non-consecutive entries for column {}.", column_name));
      }
    }

    if (is_inside_integer_) {
      state.cols.back().is_integer = true;
    }

    if (str::all_spaces(record.fields[2])) {
      throw std::runtime_error(
          "There should be at least one row specified for column data record.");
    }

    parse_coefficient(record, state, 0);

    if (!str::all_spaces(record.fields[4])) {
      parse_coefficient(record, state, 1);
    }
  }
};

}  // namespace mps
