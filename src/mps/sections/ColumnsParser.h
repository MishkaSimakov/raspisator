#pragma once

#include <cassert>
#include <format>

#include "SectionParser.h"
#include "linear/FieldTraits.h"

namespace mps {

template <typename Field>
class ColumnsParser final : public SectionParser<Field> {
  void parse_coefficient(const DataRecord& record,
                         MPSParsingState<Field>& state, size_t index) {
    assert(index == 0 || index == 1);

    const auto row = record.fields[2 + 2 * index];

    const auto itr = state.rows_map.find(row);
    if (itr == state.rows_map.end()) {
      throw std::runtime_error(std::format("Unknown row name: {}", row));
    }

    if (str::all_spaces(record.fields[3 + 2 * index])) {
      throw std::runtime_error("Coefficient must be specified.");
    }

    const size_t row_index = itr->second;
    state.cols.back().values.emplace_back(
        row_index,
        FieldTraits<Field>::from_string(record.fields[3 + 2 * index]));
  }

 public:
  bool has_field_1() const override { return false; }

  void parse(const DataRecord& record, MPSParsingState<Field>& state) override {
    const auto column_name = record.fields[1];

    if (state.cols.empty() || state.cols.back().name != column_name) {
      // new column name
      const bool inserted = state.add_col(column_name);
      if (!inserted) {
        throw std::runtime_error(
            "After a matrix element is specified for a column, all other "
            "nonzero elements in that same column should be specified.");
      }
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
