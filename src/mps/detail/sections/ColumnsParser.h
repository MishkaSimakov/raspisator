#pragma once

#include <cassert>
#include <format>

#include "SectionParser.h"
#include "linear/FieldTraits.h"

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

  void parse_marker(const DataRecord& record, MPSParsingState<Field>& state) {
    const auto name = record.fields[1];

    if (state.cols_map.contains(name)) {
      throw std::runtime_error(
          std::format("Marker name must differ from the preceding and "
                      "succeeding column names. Name {} is duplicated.",
                      name));
    }

    state.markers.insert(std::string(name));

    const auto type = parse_marker_type(record.fields[3]);

    switch (type) {
      case MarkerType::INTEGER_BEGIN:
        if (is_inside_integer_) {
          throw std::runtime_error(
              "Integer section begin marker while the previous section has not "
              "ended yet.");
        }

        is_inside_integer_ = true;
        must_change_column_ = true;
        break;
      case MarkerType::INTEGER_END:
        if (!is_inside_integer_) {
          throw std::runtime_error(
              "Integer section end marker while not in integer section.");
        }

        is_inside_integer_ = false;
        must_change_column_ = true;
        break;
      default:
        throw std::runtime_error("Unknown marker type.");
    }
  }

  void parse_column(const DataRecord& record, MPSParsingState<Field>& state) {
    const auto name = record.fields[1];

    if (state.markers.contains(name)) {
      throw std::runtime_error(
          std::format("Marker name must differ from the preceding and "
                      "succeeding column names. Name {} is duplicated.",
                      name));
    }

    if (must_change_column_ && !state.cols.empty() &&
        state.cols.back().name == name) {
      throw std::runtime_error(
          std::format("Data records for column {} exist both inside and "
                      "outside of integer section.",
                      name));
    }

    if (state.cols.empty() || state.cols.back().name != name) {
      // new column name
      const bool inserted = state.add_col(name);
      if (!inserted) {
        throw std::runtime_error(
            std::format("Non-consecutive entries for column {}.", name));
      }

      must_change_column_ = false;
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
      throw std::runtime_error("Integer section is not closed.");
    }
  }
};

}  // namespace mps::detail
