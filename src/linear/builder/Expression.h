#pragma once

#include <unordered_map>

#include "Variable.h"

template <typename Field>
class Expression {
  Field shift_;
  std::unordered_map<size_t, Field> variables_;

  // removes variables with zero coefficient
  void cleanup() {
    for (auto itr = variables_.begin(); itr != variables_.end();) {
      if (itr->second == 0) {
        itr = variables_.erase(itr);
      } else {
        ++itr;
      }
    }
  }

 public:
  Expression() = default;

  // NOLINTNEXTLINE(google-explicit-constructor)
  Expression(Variable<Field> var) : shift_(0) {
    variables_.emplace(var.id_, 1);
  }

  // NOLINTNEXTLINE(google-explicit-constructor)
  Expression(Field value) : shift_(std::move(value)) {}

  // an expression is constant if it doesn't contain any variables
  bool is_constant() const { return variables_.empty(); }

  Expression& operator+=(const Expression& other) {
    shift_ += other.shift_;

    for (const auto& [id, coef] : other.variables_) {
      auto itr = variables_.find(id);

      if (itr == variables_.end()) {
        variables_.emplace(id, coef);
      } else {
        itr->second += coef;
      }
    }

    cleanup();

    return *this;
  }

  Expression& operator-=(const Expression& other) {
    shift_ -= other.shift_;

    for (const auto& [id, coef] : other.variables_) {
      auto itr = variables_.find(id);

      if (itr == variables_.end()) {
        variables_.emplace(id, -coef);
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

  friend ProblemBuilder<Field>;
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
  left -= right;
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
