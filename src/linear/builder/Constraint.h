#pragma once

#include "Expression.h"

enum class ConstraintType { EQUAL, LESS_OR_EQUAL, GREATER_OR_EQUAL };

enum class EvaluationResult { TRUE, FALSE, DONT_KNOW };

template <typename Field>
class Constraint {
  Expression<Field> lhs_;
  Expression<Field> rhs_;
  ConstraintType type_;

 public:
  Constraint(Expression<Field> lhs, Expression<Field> rhs, ConstraintType type)
      : lhs_(std::move(lhs)), rhs_(std::move(rhs)), type_(type) {}

  bool is_constant() const { return lhs_.is_constant() && rhs_.is_constant(); }

  EvaluationResult evaluate() const {
    if (!is_constant()) {
      return EvaluationResult::DONT_KNOW;
    }

    switch (type_) {
      case ConstraintType::EQUAL:
        return lhs_.shift_ == rhs_.shift_ ? EvaluationResult::TRUE
                                          : EvaluationResult::FALSE;
      case ConstraintType::LESS_OR_EQUAL:
        return lhs_.shift_ <= rhs_.shift_ ? EvaluationResult::TRUE
                                          : EvaluationResult::FALSE;
      case ConstraintType::GREATER_OR_EQUAL:
        return lhs_.shift_ >= rhs_.shift_ ? EvaluationResult::TRUE
                                          : EvaluationResult::FALSE;
    }

    std::unreachable();
  }

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
