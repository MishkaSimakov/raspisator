#pragma once

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
  std::vector<std::pair<std::string, Field>> variables;
  Field rhs{0};
  std::optional<Field> range = std::nullopt;

  explicit Row(RowSense type, std::string_view name)
      : type(type), name(name) {}
};

template <typename Field>
struct MPSParsingState {
  ObjectiveType objective_ = ObjectiveType::MINIMIZE;

  std::vector<Row<Field>> rows;
  std::unordered_map<std::string_view, size_t> rows_map;
};

}  // namespace mps
