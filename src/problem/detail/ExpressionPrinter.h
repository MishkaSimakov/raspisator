#pragma once

#include <iostream>

#include "problem/Bound.h"

namespace problem::detail {

class ExpressionPrinter {
  std::ostream& os_;

  bool printed_name_{false};
  bool printed_expr_{false};

 public:
  explicit ExpressionPrinter(std::ostream& os) : os_(os) {}

  void name(std::string_view name) {
    if (printed_name_) {
      throw std::logic_error("Name was already printed.");
    }

    os_ << name << " = ";
    printed_name_ = true;
  }

  template <typename Field>
  void print(const Field& value, std::string_view name = std::string_view{}) {
    using std::abs;

    if (printed_expr_) {
      if (value >= 0) {
        os_ << " + " << value;
      } else {
        os_ << " - " << abs(value);
      }
    } else {
      os_ << value;
    }

    if (!name.empty()) {
      os_ << "*" << name;
    }

    printed_expr_ = true;
  }
};

}  // namespace problem::detail
