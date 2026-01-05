#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <sstream>
#include <unordered_map>

#include "ToMatrices.h"
#include "linear/matrix/Matrix.h"
#include "utils/String.h"

template <typename Field>
class MPSReader {
  enum class SectionType {
    NAME,
    ROWS,
    COLUMNS,
    RHS,
    BOUNDS,
    RANGES,
    OBJECT,
    ENDATA
  };

  enum class ObjectiveType { MINIMIZE, MAXIMIZE };

  enum class RowType { LESS_THAN, GREATER_THAN, EQUAL, OBJECTIVE };

  struct Row {
    RowType type;
    std::vector<std::pair<std::string, Field>> variables;
    Field rhs{0};
    std::optional<Field> range = std::nullopt;

    explicit Row(RowType type) : type(type) {}
  };

  struct Bounds {
    std::optional<Field> lower = 0;
    std::optional<Field> upper = std::nullopt;
  };

  static std::array<std::string, 6> get_parts(const std::string& str) {
    constexpr size_t kNumFields = 6;
    constexpr size_t kFieldStartPos[kNumFields] = {1, 4, 14, 24, 39, 49};
    constexpr size_t kFieldLength[kNumFields] = {2, 8, 8, 12, 8, 12};

    std::array<std::string, 6> result;

    for (size_t i = 0; i < kNumFields; ++i) {
      size_t start = kFieldStartPos[i];
      size_t length = kFieldLength[i];

      if (start >= str.size()) {
        break;
      }
      if (start + length > str.size()) {
        length = std::string::npos;
      }

      result[i] = str::rtrim(str.substr(start, length));
    }

    return result;
  }

  static std::optional<SectionType> read_header_card(const std::string& line) {
    std::vector headers = {
        std::pair{"NAME", SectionType::NAME},
        std::pair{"ROWS", SectionType::ROWS},
        std::pair{"COLUMNS", SectionType::COLUMNS},
        std::pair{"RHS", SectionType::RHS},
        std::pair{"BOUNDS", SectionType::BOUNDS},
        std::pair{"RANGES", SectionType::RANGES},
        std::pair{"OBJECT", SectionType::OBJECT},
        std::pair{"ENDATA", SectionType::ENDATA},
    };

    for (const auto& [name, value] : headers) {
      if (line.starts_with(name)) {
        return value;
      }
    }

    throw std::runtime_error("Unknown header card.");
  }

  static RowType decode_row_type(const std::string& type) {
    auto trimmed = str::trim(type);

    switch (trimmed[0]) {
      case 'E':
        return RowType::EQUAL;
      case 'L':
        return RowType::LESS_THAN;
      case 'G':
        return RowType::GREATER_THAN;
      case 'N':
        return RowType::OBJECTIVE;
    }

    throw std::runtime_error("Unknown row type.");
  }

  static Field parse_field(const std::string& str) {
    Field result;
    std::stringstream iss(str);
    iss >> result;

    return result;
  }

  void store_bounds(const std::string& type, Field value, Bounds& bounds) {
    if (type == "LO") {
      bounds.lower = value;
    } else if (type == "UP") {
      bounds.upper = value;
    } else if (type == "FX") {
      bounds.lower = value;
      bounds.upper = value;
    } else if (type == "FR") {
      bounds.lower = std::nullopt;
      bounds.upper = std::nullopt;
    } else if (type == "MI") {
      bounds.lower = std::nullopt;
      bounds.upper = 0;
    } else if (type == "PL") {
      bounds.lower = 0;
      bounds.upper = std::nullopt;
    } else {
      throw std::runtime_error("Unsupported bound type.");
    }
  }

  bool should_skip_line(const std::string& line) {
    if (line.empty() || line[0] == '*') {
      return true;
    }

    if (str::ltrim(line).empty()) {
      return true;
    }

    return false;
  }

  ObjectiveType objective_ = ObjectiveType::MINIMIZE;
  std::unordered_map<std::string, Row> rows_;
  std::unordered_map<std::string, Bounds> bounds_;

 public:
  void read(const std::filesystem::path& filepath) {
    std::ifstream is(filepath);
    std::string line;

    std::optional<SectionType> current_section = std::nullopt;

    while (std::getline(is, line)) {
      if (should_skip_line(line)) {
        continue;
      }

      // header card
      if (line[0] != ' ') {
        current_section = read_header_card(line);

        if (current_section == SectionType::ENDATA) {
          break;
        }

        continue;
      }

      // otherwise proceed with parsing of the current section
      if (!current_section.has_value()) {
        throw std::runtime_error("Section data must be inside of section.");
      }

      auto parts = get_parts(line);

      if (current_section == SectionType::ROWS) {
        rows_.emplace(parts[1], Row(decode_row_type(parts[0])));
      } else if (current_section == SectionType::COLUMNS) {
        std::string variable_name = parts[1];
        bounds_.emplace(variable_name, Bounds{});

        for (size_t i = 2; i < parts.size(); i += 2) {
          std::string row_name = parts[i];

          if (row_name.empty()) {
            break;
          }

          rows_.at(row_name).variables.emplace_back(variable_name,
                                                    parse_field(parts[i + 1]));
        }
      } else if (current_section == SectionType::RHS) {
        // rhs vector name is ignored

        for (size_t i = 2; i < parts.size(); i += 2) {
          std::string row_name = parts[i];
          if (row_name.empty()) {
            continue;
          }

          rows_.at(row_name).rhs = parse_field(parts[i + 1]);
        }
      } else if (current_section == SectionType::BOUNDS) {
        std::string type = parts[0];
        std::string variable_name = parts[2];
        Field value = parse_field(parts[3]);

        store_bounds(type, value, bounds_.at(variable_name));
      } else if (current_section == SectionType::RANGES) {
        for (size_t i = 2; i < parts.size(); i += 2) {
          std::string row_name = parts[i];
          if (row_name.empty()) {
            continue;
          }

          rows_.at(row_name).range = parse_field(parts[i + 1]);
        }
      } else if (current_section == SectionType::OBJECT) {
        // skip for now
        continue;
      } else {
        throw std::runtime_error("Unknown section type.");
      }
    }
  }

  // generates a problem suitable for simplex method:
  // c x -> max, s.t. Ax = b, l <= x <= u
  MILPProblem<Field> get_canonical_representation() {
    MILPProblem<Field> result;
    std::unordered_map<std::string, Variable<Field>> variables;

    for (const Row& row : rows_ | std::views::values) {
      for (const std::string& name : row.variables | std::views::keys) {
        if (!variables.contains(name)) {
          auto lower = bounds_.at(name).lower;
          auto upper = bounds_.at(name).upper;

          auto var = result.new_variable(name, VariableType::REAL, lower ? *lower : -1e7, upper ? *upper : 1e7);
          variables.emplace(name, var);
        }
      }
    }

    auto objective_row = std::ranges::find_if(
        rows_, [](auto p) { return p.second.type == RowType::OBJECTIVE; });

    if (objective_row == rows_.end()) {
      throw std::runtime_error("No objective row found.");
    }

    Expression<Field> objective;
    for (auto [name, coef] : objective_row->second.variables) {
      if (objective_ == ObjectiveType::MINIMIZE) {
        coef = -coef;
      }

      objective += variables.at(name) * coef;
    }

    result.set_objective(objective);

    // process constraints
    for (const Row& row : rows_ | std::views::values) {
      if (row.type == RowType::OBJECTIVE) {
        continue;
      }

      Expression<Field> lhs;
      for (auto [name, coef] : row.variables) {
        lhs += variables.at(name) * coef;
      }

      if (!row.range) {
        if (row.type == RowType::LESS_THAN) {
          result.add_constraint(lhs <= Expression{row.rhs});
        } else if (row.type == RowType::GREATER_THAN) {
          result.add_constraint(lhs >= Expression{row.rhs});
        } else {
          result.add_constraint(lhs == Expression{row.rhs});
        }
      } else {
        Field upper;
        Field lower;

        if (row.type == RowType::LESS_THAN) {
          upper = row.rhs;
          lower = row.rhs - *row.range;
        } else if (row.type == RowType::GREATER_THAN) {
          upper = row.rhs + *row.range;
          lower = row.rhs;
        } else {
          upper = std::max(row.rhs, row.rhs + *row.range);
          lower = std::min(row.rhs, row.rhs + *row.range);
        }

        result.add_constraint(lhs <= Expression{upper});
        result.add_constraint(lhs >= Expression{lower});
      }
    }

    return result;
  }
};
