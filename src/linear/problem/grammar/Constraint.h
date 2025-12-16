#pragma once

#include <format>
#include <iostream>

#include "Expression.h"

enum class ConstraintType { EQUAL_ZERO, LESS_OR_EQUAL_ZERO };

enum class EvaluationResult { TRUE, FALSE, DONT_KNOW };

template <typename Field>
class Constraint {
 public:
  Expression<Field> expr;
  ConstraintType type;

  Constraint(Expression<Field> expr, ConstraintType type)
      : expr(expr), type(type) {}

  bool is_constant() const { return expr.is_constant(); }

  EvaluationResult evaluate() const {
    if (!is_constant()) {
      return EvaluationResult::DONT_KNOW;
    }

    switch (type) {
      case ConstraintType::EQUAL_ZERO:
        return FieldTraits<Field>::is_nonzero(expr.shift_)
                   ? EvaluationResult::FALSE
                   : EvaluationResult::TRUE;
      case ConstraintType::LESS_OR_EQUAL_ZERO:
        return FieldTraits<Field>::is_strictly_positive(expr.shift_)
                   ? EvaluationResult::FALSE
                   : EvaluationResult::TRUE;
    }

    std::unreachable();
  }
};

auto operator==(details::ExpressionLike auto&& left,
                details::ExpressionLike auto&& right) {
  return Constraint(Expression(left) - Expression(right),
                    ConstraintType::EQUAL_ZERO);
}

auto operator<=(details::ExpressionLike auto&& left,
                details::ExpressionLike auto&& right) {
  return Constraint(Expression(left) - Expression(right),
                    ConstraintType::LESS_OR_EQUAL_ZERO);
}

auto operator>=(details::ExpressionLike auto&& left,
                details::ExpressionLike auto&& right) {
  return Constraint(Expression(right) - Expression(left),
                    ConstraintType::LESS_OR_EQUAL_ZERO);
}

template <typename Field>
std::ostream& operator<<(std::ostream& os,
                         const Constraint<Field>& constraint) {
  os << constraint.expr;

  if (constraint.type == ConstraintType::EQUAL_ZERO) {
    os << " == 0";
  } else if (constraint.type == ConstraintType::LESS_OR_EQUAL_ZERO) {
    os << " <= 0";
  }

  return os;
}

template <typename Field>
struct std::formatter<Constraint<Field>, char> : std::formatter<std::string> {
  template <class FmtContext>
  auto format(const Constraint<Field>& value, FmtContext& ctx) const {
    std::ostringstream out;
    out << value;
    return std::ranges::copy(std::move(out).str(), ctx.out()).out;
  }
};
