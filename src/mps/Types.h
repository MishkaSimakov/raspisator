#pragma once

#include <deque>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace mps {

enum class Format { FREE, FIXED };

enum class SectionType {
  NAME = 0,
  ROWS,
  COLUMNS,
  RHS,
  RANGES,
  BOUNDS,
  ENDATA,

  // this value must always be the last one
  SECTIONS_COUNT
};
enum class ObjectiveType { MINIMIZE, MAXIMIZE };
enum class RowSense { LESS_THAN, GREATER_THAN, EQUAL, FREE };

template <typename Field>
struct Row {
  RowSense type;
  std::string name;
  Field rhs{0};
  std::optional<Field> range = std::nullopt;

  Row(RowSense type, std::string_view name) : type(type), name(name) {}
};

template <typename Field>
struct Variable {
  std::string name;

  std::vector<std::pair<size_t, Field>> values;

  explicit Variable(std::string_view name) : name(name) {}
};

template <typename Field>
struct MPSParsingState {
  ObjectiveType objective = ObjectiveType::MINIMIZE;

  // Row name is owned by Row class. std::deque never reallocates them.
  // rows_map uses std::string_view to reference row name.
  std::deque<Row<Field>> rows;
  std::unordered_map<std::string_view, size_t> rows_map;

  std::deque<Variable<Field>> cols;
  std::unordered_map<std::string_view, size_t> cols_map;

  bool add_row(RowSense sense, std::string_view name) {
    const size_t index = rows.size();

    const auto& row = rows.emplace_back(sense, name);
    auto [_, inserted] = rows_map.emplace(row.name, index);

    return inserted;
  }

  bool add_col(std::string_view name) {
    const size_t index = cols.size();

    const auto& col = cols.emplace_back(name);
    auto [_, inserted] = cols_map.emplace(col.name, index);

    return inserted;
  }
};

}  // namespace mps
