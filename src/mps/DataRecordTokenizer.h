#pragma once

#include <array>
#include <string_view>

#include "Types.h"

namespace mps {

constexpr size_t kFieldsCount = 6;

struct DataRecordTokens {
  std::array<std::string_view, kFieldsCount> fields;
  std::string_view comment;

  bool operator==(const DataRecordTokens&) const = default;
};

// Splits MPS data record into fields. Supports both FREE and FIXED format.
// Return fields values via std::string_view. This means that original string
// must be alive for as long as fields values are used.
// In some sections of MPS file Field 1 may be blank.
class DataRecordTokenizer {
  constexpr static size_t kFieldStartPos[kFieldsCount] = {1, 4, 14, 24, 39, 49};
  constexpr static size_t kFieldLength[kFieldsCount] = {2, 8, 8, 12, 8, 12};

  static DataRecordTokens parse_fixed(std::string_view record,
                                      bool has_field_1) {
    DataRecordTokens result;

    for (size_t i = 0; i < kFieldsCount; ++i) {
      const size_t start = kFieldStartPos[i];
      const size_t length = kFieldLength[i];

      if (start >= record.size()) {
        break;
      }

      if ((i == 2 || i == 4) && record[start] == '$') {
        result.comment = record.substr(start + 1, std::string_view::npos);
        break;
      }

      result.fields[i] = record.substr(start, length);
    }

    if (!has_field_1) {
      for (const char c : result.fields[0]) {
        if (!is_space(c)) {
          throw std::runtime_error(
              "In the current MPS section Field 1 must be empty.");
        }
      }
    }

    return result;
  }

  static DataRecordTokens parse_free(std::string_view record,
                                     bool has_field_1) {
    std::array<std::string_view, kFieldsCount> fields;
    std::string_view comment;

    size_t current_field_index = has_field_1 ? 0 : 1;
    size_t current_field_begin = 0;

    for (size_t i = 0; i < record.size(); ++i) {
      if (i > 0 && !is_space(record[i]) && is_space(record[i - 1])) {
        // field begin

        if (current_field_index >= kFieldsCount) {
          throw std::runtime_error("Too many fields in MPS data record.");
        }

        if ((current_field_index == 2 || current_field_index == 4) &&
            record[i] == '$') {
          comment = record.substr(i + 1, std::string_view::npos);
          break;
        }

        current_field_begin = i;
      }

      if (!is_space(record[i]) &&
          (i + 1 == record.size() || is_space(record[i + 1]))) {
        // field end
        fields[current_field_index] =
            record.substr(current_field_begin, i - current_field_begin + 1);

        ++current_field_index;
      }
    }

    return DataRecordTokens{
        .fields = fields,
        .comment = comment,
    };
  }

  static bool is_space(char symbol) {
    return std::isspace(static_cast<unsigned char>(symbol)) != 0;
  }

 public:
  static DataRecordTokens parse(std::string_view record, Format format,
                                bool has_field_1) {
    if (!is_space(record[0])) {
      throw std::runtime_error("Column 1 in MPS data record must be empty.");
    }

    switch (format) {
      case Format::FREE:
        return parse_free(record, has_field_1);
      case Format::FIXED:
        return parse_fixed(record, has_field_1);
      default:
        throw std::runtime_error("Unknown MPS format.");
    }
  }
};

}  // namespace mps
