#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <sstream>
#include <unordered_map>

#include "Matrix.h"

template <typename Field>
class MPSReader {
  enum class SectionType { NAME, ROWS, COLUMNS, RHS, ENDATA };

  enum class ObjectiveType { MINIMIZE, MAXIMIZE };

  enum class RowType { LESS_THAN, GREATER_THAN, EQUAL, OBJECTIVE };

  struct Row {
    RowType type;
    std::vector<std::pair<std::string, Field>> variables;
    Field rhs{0};

    explicit Row(RowType type) : type(type) {}
  };

  static std::vector<std::string> split(const std::string& str) {
    std::vector<std::string> result;

    for (size_t i = 0; i < str.size(); ++i) {
      if (std::isspace(str[i]) != 0) {
        continue;
      }

      if (i == 0 || std::isspace(str[i - 1]) != 0) {
        result.emplace_back(1, str[i]);
      } else {
        result.back().push_back(str[i]);
      }
    }

    return result;
  }

  static std::optional<SectionType> read_section_header(
      const std::string& line) {
    if (line.starts_with("NAME")) {
      return SectionType::NAME;
    }
    if (line.starts_with("ROWS")) {
      return SectionType::ROWS;
    }
    if (line.starts_with("COLUMNS")) {
      return SectionType::COLUMNS;
    }
    if (line.starts_with("RHS")) {
      return SectionType::RHS;
    }
    if (line.starts_with("ENDATA")) {
      return SectionType::ENDATA;
    }

    return std::nullopt;
  }

  static RowType decode_row_type(const std::string& type) {
    switch (type[0]) {
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

  ObjectiveType objective_ = ObjectiveType::MINIMIZE;
  std::unordered_map<std::string, Row> rows_;

 public:
  void read(const std::filesystem::path& filepath) {
    std::ifstream is(filepath);
    std::string line;

    std::optional<SectionType> current_section = std::nullopt;

    while (std::getline(is, line)) {
      auto parts = split(line);

      // skip comments, NAME section and empty lines
      if (parts.empty() || parts[0].starts_with("*") ||
          parts[0].starts_with("NAME")) {
        continue;
      }

      auto section_header = read_section_header(parts[0]);

      // stop if reached end
      if (section_header == SectionType::ENDATA) {
        break;
      }

      // switch section
      if (section_header.has_value()) {
        current_section = section_header;
        continue;
      }

      // otherwise proceed with parsing of the current section
      if (!current_section.has_value()) {
        throw std::runtime_error("Section data must be inside of section.");
      }

      if (current_section == SectionType::ROWS) {
        rows_.emplace(parts[1], Row(decode_row_type(parts[0])));
      } else if (current_section == SectionType::COLUMNS) {
        std::string variable_name = parts[0];

        for (size_t i = 1; i < parts.size(); i += 2) {
          std::string row_name = parts[i];

          // this is too awkward, fix this ASAP, please
          Field coef;
          std::stringstream iss(parts[i + 1]);
          iss >> coef;

          rows_.at(row_name).variables.emplace_back(variable_name, coef);
        }
      } else if (current_section == SectionType::RHS) {
        // rhs vector name is ignored

        for (size_t i = 1; i < parts.size(); i += 2) {
          std::string row_name = parts[i];

          // this is too awkward, fix this ASAP, please
          Field coef;
          std::stringstream iss(parts[i + 1]);
          iss >> coef;

          rows_.at(row_name).rhs = coef;
        }
      }
    }
  }

  // generates a problem suitable for simplex method:
  // c x -> max, s.t. Ax = b, x >= 0
  std::tuple<Matrix<Field>, Matrix<Field>, Matrix<Field>>
  get_canonical_representation() {
    size_t additional_variables_cnt = 0;

    std::unordered_map<std::string, size_t> variables_enumeration;
    for (const Row& row : rows_ | std::views::values) {
      if (row.type == RowType::LESS_THAN || row.type == RowType::GREATER_THAN) {
        ++additional_variables_cnt;
      }

      for (const std::string& var : row.variables | std::views::keys) {
        if (!variables_enumeration.contains(var)) {
          variables_enumeration.emplace(var, variables_enumeration.size());
        }
      }
    }

    size_t total_variables_cnt =
        variables_enumeration.size() + additional_variables_cnt;

    // process objective
    Matrix<Field> c(1, total_variables_cnt);

    auto objective_row = std::ranges::find_if(
        rows_, [](auto p) { return p.second.type == RowType::OBJECTIVE; });

    if (objective_row == rows_.end()) {
      throw std::runtime_error("No objective row found.");
    }

    for (auto [name, coef] : objective_row->second.variables) {
      if (objective_ == ObjectiveType::MINIMIZE) {
        coef = -coef;
      }

      c[0, variables_enumeration[name]] = coef;
    }

    // process constraints
    auto A = Matrix<Field>(rows_.size() - 1, total_variables_cnt, 0);
    auto b = Matrix<Field>(rows_.size() - 1, 1, 0);

    size_t current_additional_variable = variables_enumeration.size();
    size_t i = 0;
    for (const Row& row : rows_ | std::views::values) {
      if (row.type == RowType::OBJECTIVE) {
        continue;
      }

      for (auto [var, coef] : row.variables) {
        A[i, variables_enumeration[var]] = coef;
      }

      if (row.type != RowType::EQUAL) {
        A[i, current_additional_variable] =
            row.type == RowType::LESS_THAN ? 1 : -1;
        ++current_additional_variable;
      }

      b[i, 0] = row.rhs;

      ++i;
    }

    return {A, b, c};
  }
};
