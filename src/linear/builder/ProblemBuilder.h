#pragma once

#include "Constraint.h"
#include "Expression.h"
#include "Variable.h"

template <typename Field>
class ProblemBuilder {
  std::vector<std::string> variables_names_;

  std::vector<Constraint<Field>> constraints_;
  Expression<Field> objective_;

  void print_expression(std::ostream& os, const Expression<Field>& expr) const {
    for (const auto& [var, coef] : expr.variables_) {
      os << coef << "*" << variables_names_.at(var) << " + ";
    }

    os << expr.shift_;
  }

  void print_constraint(std::ostream& os,
                        const Constraint<Field>& constraint) const {
    print_expression(os, constraint.lhs_);

    if (constraint.type_ == ConstraintType::EQUAL) {
      os << " = ";
    } else if (constraint.type_ == ConstraintType::LESS_OR_EQUAL) {
      os << " <= ";
    } else {
      os << " >= ";
    }

    print_expression(os, constraint.rhs_);
  }

 public:
  // transforms all constraints into equalities (by adding slack variables)
  // moves all variables into lhs, all constants into rhs
  // removes all useless constraints
  // removes linearly dependent rows
  void normalize() {
    size_t slack_index = 0;

    for (Constraint<Field>& constraint : constraints_) {
      if (constraint.type_ == ConstraintType::EQUAL) {
        continue;
      }

      if (constraint.type_ == ConstraintType::GREATER_OR_EQUAL) {
        constraint.lhs_ *= -1;
        constraint.rhs_ *= -1;
      }

      // now it is <= constraint
      auto slack_var = new_variable(fmt::format("slack({})", slack_index));
      ++slack_index;

      constraint.lhs_ += slack_var;
      constraint.type_ = ConstraintType::EQUAL;
    }

    //
    for (Constraint<Field>& constraint : constraints_) {
      constraint.lhs_ -= constraint.rhs_;
      constraint.rhs_ = Expression<Field>(0);
      std::swap(constraint.lhs_.shift_, constraint.rhs_.shift_);
    }
  }

 public:
  explicit ProblemBuilder() = default;

  Variable<Field> new_variable(std::string_view name) {
    size_t index = variables_names_.size();
    variables_names_.emplace_back(name);

    return Variable<Field>{index};
  }

  void add_constraint(Constraint<Field> constraint) {
    constraints_.push_back(std::move(constraint));
  }

  void set_objective(Expression<Field> objective) {
    objective_ = std::move(objective);
  }

  std::tuple<Matrix<Field>, Matrix<Field>, Matrix<Field>> get_matrices() {}

  //
  friend std::ostream& operator<<(std::ostream& os,
                                  const ProblemBuilder& builder) {
    os << "problem builder with " << builder.variables_names_.size()
       << " variables:\n";

    for (const Constraint<Field>& constraint : builder.constraints_) {
      builder.print_constraint(os, constraint);
      os << "\n";
    }

    return os;
  }
};
