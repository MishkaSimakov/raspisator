#pragma once

#include "fwd.h"

enum class ConstraintType { EQUAL, LESS_OR_EQUAL, GREATER_OR_EQUAL };

template <typename Field>
class Constraint {
  Expression<Field> lhs_;
  Expression<Field> rhs_;
  ConstraintType type_;

 public:
  Constraint(Expression<Field> lhs, Expression<Field> rhs, ConstraintType type)
      : lhs_(std::move(lhs)), rhs_(std::move(rhs)), type_(type) {}

  friend ProblemBuilder<Field>;
};

// TODO: redundant Expression copy here
auto operator==(details::ExpressionLike auto&& left,
                details::ExpressionLike auto&& right) {
  return Constraint(Expression(left), Expression(right), ConstraintType::EQUAL);
}

auto operator<=(details::ExpressionLike auto&& left,
                details::ExpressionLike auto&& right) {
  return Constraint(Expression(left), Expression(right),
                    ConstraintType::LESS_OR_EQUAL);
}

auto operator>=(details::ExpressionLike auto&& left,
                details::ExpressionLike auto&& right) {
  return Constraint(Expression(left), Expression(right),
                    ConstraintType::GREATER_OR_EQUAL);
}
