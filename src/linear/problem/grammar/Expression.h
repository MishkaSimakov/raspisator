#pragma once

#include <iostream>
#include <unordered_map>

#include "Variable.h"
#include "linear/FieldTraits.h"

template <typename Field>
class Expression {
  Field shift_;
  std::unordered_map<std::string, Field> variables_;

  // removes variables with zero coefficient
  void cleanup() {
    for (auto itr = variables_.begin(); itr != variables_.end();) {
      if (FieldTraits<Field>::is_nonzero(itr->second)) {
        ++itr;
      } else {
        itr = variables_.erase(itr);
      }
    }
  }

 public:
  Expression() = default;

  // NOLINTNEXTLINE(google-explicit-constructor)
  Expression(Variable<Field> var) : shift_(0) {
    variables_.emplace(var.get_name(), 1);
  }

  // NOLINTNEXTLINE(google-explicit-constructor)
  Expression(Field value) : shift_(std::move(value)) {}

  const auto& get_variables() const { return variables_; }
  auto& get_variables() { return variables_; }

  // an expression is constant if it doesn't contain any variables
  bool is_constant() const { return variables_.empty(); }

  Expression& operator+=(const Expression& other) {
    shift_ += other.shift_;

    for (const auto& [var, coef] : other.variables_) {
      auto itr = variables_.find(var);

      if (itr == variables_.end()) {
        variables_.emplace(var, coef);
      } else {
        itr->second += coef;
      }
    }

    cleanup();

    return *this;
  }

  Expression& operator-=(const Expression& other) {
    shift_ -= other.shift_;

    for (const auto& [var, coef] : other.variables_) {
      auto itr = variables_.find(var);

      if (itr == variables_.end()) {
        variables_.emplace(var, -coef);
      } else {
        itr->second -= coef;
      }
    }

    cleanup();

    return *this;
  }

  Expression& operator*=(const Field& alpha) {
    for (auto& [_, coef] : variables_) {
      coef *= alpha;
    }
    shift_ *= alpha;

    cleanup();

    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, const Expression& expr) {
    bool has_shift = expr.shift_ != 0;

    for (auto itr = expr.variables_.begin(); itr != expr.variables_.end();
         ++itr) {
      bool is_last = std::next(itr) == expr.variables_.end();

      if (itr->second == -1) {
        os << "-";
      } else if (itr->second != 1) {
        os << itr->second << "*";
      }

      os << itr->first;

      if (!is_last || has_shift) {
        os << " + ";
      }
    }

    if (has_shift || expr.is_constant()) {
      os << expr.shift_;
    }

    return os;
  }

  std::optional<Expression> express(Variable<Field> variable) const {
    auto itr = variables_.find(variable);

    if (itr == variables_.end()) {
      return std::nullopt;
    }

    Expression copy = *this;

    for (auto& coef : copy.variables_ | std::views::values) {
      coef /= itr->second;
    }

    return copy;
  }

  Field get_shift() const { return shift_; }

  friend MILPProblem<Field>;
  friend Constraint<Field>;
};

auto operator+(details::ExpressionLike auto&& left,
               details::ExpressionLike auto&& right) {
  Expression copy = left;
  copy += right;
  return copy;
}

auto operator-(details::ExpressionLike auto&& left,
               details::ExpressionLike auto&& right) {
  Expression copy = left;
  copy -= right;
  return copy;
}

template <typename Field>
auto operator*(const Field& alpha, details::ExpressionLike auto&& expr) {
  Expression copy = expr;
  copy *= alpha;
  return copy;
}

template <typename Field>
auto operator*(details::ExpressionLike auto&& expr, const Field& alpha) {
  Expression copy = expr;
  copy *= alpha;
  return copy;
}

auto operator-(details::ExpressionLike auto&& expr) {
  Expression copy = expr;
  copy *= -1;
  return copy;
}
