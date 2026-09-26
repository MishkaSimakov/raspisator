#pragma once

#include <deque>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "problem/Bound.h"

namespace mps::detail {

enum class SectionType {
  NAME = 0,
  ROWS,
  COLUMNS,
  RHS,
  RANGES,
  BOUNDS,
  ENDATA,

  OBJECT,

  // this value must always be the last one
  SECTIONS_COUNT
};

inline std::string section_type_to_string(SectionType type) {
  switch (type) {
    case SectionType::NAME:
      return "NAME";
    case SectionType::ROWS:
      return "ROWS";
    case SectionType::COLUMNS:
      return "COLUMNS";
    case SectionType::RHS:
      return "RHS";
    case SectionType::RANGES:
      return "RANGES";
    case SectionType::BOUNDS:
      return "BOUNDS";
    case SectionType::ENDATA:
      return "ENDATA";
    case SectionType::OBJECT:
      return "OBJECT";
    default:
      return "<UNKNOWN>";
  }
}

enum class ObjectiveType { MINIMIZE, MAXIMIZE };
enum class RowSense { LESS_THAN, GREATER_THAN, EQUAL, FREE };

template <typename Field>
struct Row {
  RowSense type;
  std::string name;
  std::optional<Field> rhs = std::nullopt;
  std::optional<Field> range = std::nullopt;

  Row(RowSense type, std::string_view name) : type(type), name(name) {}
};

template <typename Field>
struct Variable {
  std::string name;
  std::map<size_t, Field> values;

  bool is_integer = false;

  // Bound is (std::nullopt, std::nullopt) by default. lower_specified and
  // upper_specified encode whether lower or upper bound was specified in the
  // MPS file. Default bound behaviour (e.g. default bound is (0, +inf)) is
  // later reconstructed through these values.
  Bound<Field> bound = {std::nullopt, std::nullopt};
  bool lower_specified = false;
  bool upper_specified = false;

  explicit Variable(std::string_view name) : name(name) {}
};

template <typename Field>
struct MPSParsingState {
  std::string problem_name;

  ObjectiveType objective = ObjectiveType::MINIMIZE;

  // Row name is owned by Row class. std::deque never reallocates them.
  // rows_map uses std::string_view to reference row name.
  std::deque<Row<Field>> rows;
  std::unordered_map<std::string_view, size_t> rows_map;

  std::deque<Variable<Field>> cols;
  std::unordered_map<std::string_view, size_t> cols_map;

  // In MPS multiple RHS, RANGES, and BOUNDS vectors may be specified, but only
  // the first one must be selected. These fields capture the name of the first
  // vector name in each section.
  std::optional<std::string> rhs_vector_name = std::nullopt;
  std::optional<std::string> ranges_vector_name = std::nullopt;
  std::optional<std::string> bounds_vector_name = std::nullopt;

  // Markers names must differ from all columns names. This field records all
  // markers names.
  // Ordered set is used to allow for transparent comparator. This allows to
  // call find with std::string_view.
  std::set<std::string, std::less<>> markers;

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

}  // namespace mps::detail
