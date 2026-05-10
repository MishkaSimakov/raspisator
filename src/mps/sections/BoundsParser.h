#pragma once

#include <optional>
#include <string_view>

#include "SectionParser.h"
#include "linear/model/Bound.h"

namespace mps {

template <typename Field>
class BoundsParser final : public SectionParser<Field> {
  void set_lower(Variable<Field>& variable, std::optional<Field> value) {
    if (variable.lower_specified) {
      throw std::runtime_error(std::format(
          "Lower bound is specified twice for variable: {}", variable.name));
    }

    variable.bound.lower = value;
    variable.lower_specified = true;
  }

  void set_upper(Variable<Field>& variable, std::optional<Field> value) {
    if (variable.upper_specified) {
      throw std::runtime_error(std::format(
          "Upper bound is specified twice for variable: {}", variable.name));
    }

    variable.bound.upper = value;
    variable.upper_specified = true;

    if (!variable.lower_specified && value && *value < 0) {
      variable.bound.lower = std::nullopt;
    }
  }

  void set_integer(Variable<Field>& variable) { variable.is_integer = true; }

  void parse_bound_with_value(Variable<Field>& variable, std::string_view type,
                              Field value) {
    if (type == "LO") {
      set_lower(variable, value);
    } else if (type == "LI") {
      set_lower(variable, value);
      set_integer(variable);
    } else if (type == "UP") {
      set_upper(variable, value);
    } else if (type == "UI") {
      set_upper(variable, value);
      set_integer(variable);
    } else if (type == "FX") {
      set_lower(variable, value);
      set_upper(variable, value);
    } else {
      throw std::runtime_error(std::format("Unknown bound type: {}.", type));
    }
  }

  void parse_bound_without_value(Variable<Field>& variable,
                                 std::string_view type) {
    if (type == "FR") {
      set_lower(variable, std::nullopt);
      set_upper(variable, std::nullopt);
    } else if (type == "MI") {
      set_lower(variable, std::nullopt);
    } else if (type == "PL") {
      set_upper(variable, std::nullopt);
    } else if (type == "BV") {
      set_lower(variable, 0);
      set_upper(variable, 1);
      set_integer(variable);
    } else {
      throw std::runtime_error(std::format("Unknown bound type: {}.", type));
    }
  }

 public:
  bool has_field_1() const override { return true; }

  void parse(const DataRecord& record, MPSParsingState<Field>& state) override {
    const auto vector_name = record.fields[1];

    if (!state.bounds_vector_name) {
      state.bounds_vector_name = vector_name;
    }

    if (*state.bounds_vector_name != vector_name) {
      return;
    }

    auto itr = state.cols_map.find(record.fields[2]);
    if (itr == state.cols_map.end()) {
      throw std::runtime_error(
          std::format("Unknown variable name: {}.", record.fields[2]));
    }

    const auto type = record.fields[0];

    if (type == "LO" || type == "LI" || type == "UP" || type == "UI" ||
        type == "FX") {
      const Field value = FieldTraits<Field>::from_string(record.fields[3]);
      parse_bound_with_value(state.cols[itr->second], type, value);
    } else {
      parse_bound_without_value(state.cols[itr->second], type);
    }
  }
};

}  // namespace mps
