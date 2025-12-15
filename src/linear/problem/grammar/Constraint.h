#pragma once

#include <format>
#include <iostream>

#include "Expression.h"

enum class ConstraintType { EQUAL, LESS_OR_EQUAL };

enum class EvaluationResult { TRUE, FALSE, DONT_KNOW };

// Constraint is stored in a normalized form. This means that:
// 1. All variables are in left hand side
// 2. All constants are in right hand side
// 3. Constraint sign is either "=" or "<=" (>= is transformed into <=)
template <typename Field>
class Constraint {
 public:
  Expression<Field> lhs;
  Field rhs;
  ConstraintType type;

  Constraint(Expression<Field> lhs, Expression<Field> rhs, ConstraintType type)
      : lhs(lhs - rhs), rhs(0), type(type) {
    this->rhs = -this->lhs.shift_;
    this->lhs.shift_ = 0;
  }

  bool is_constant() const { return lhs.is_constant(); }

  EvaluationResult evaluate() const {
    if (!is_constant()) {
      return EvaluationResult::DONT_KNOW;
    }

    switch (type) {
      case ConstraintType::EQUAL:
        return 0 == rhs ? EvaluationResult::TRUE : EvaluationResult::FALSE;
      case ConstraintType::LESS_OR_EQUAL:
        return 0 <= rhs ? EvaluationResult::TRUE : EvaluationResult::FALSE;
    }

    std::unreachable();
  }
};

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
  return Constraint(Expression(right), Expression(left),
                    ConstraintType::LESS_OR_EQUAL);
}

template <typename Field>
std::ostream& operator<<(std::ostream& os,
                         const Constraint<Field>& constraint) {
  os << constraint.lhs;

  if (constraint.type == ConstraintType::EQUAL) {
    os << " == ";
  } else if (constraint.type == ConstraintType::LESS_OR_EQUAL) {
    os << " <= ";
  }

  os << constraint.rhs;

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
