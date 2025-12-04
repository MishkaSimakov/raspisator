#pragma once

#include <unordered_set>

#include "Constraint.h"
#include "Expression.h"
#include "Variable.h"
#include "linear/model/MILP.h"

template <typename Field>
class ProblemBuilder {
  struct VariableInfo {
    VariableType type;
    Field lower_bound;
    Field upper_bound;
  };

  std::vector<std::pair<std::string, VariableInfo>> variables_;

  std::vector<Constraint<Field>> constraints_;
  Expression<Field> objective_;

  void print_expression(std::ostream& os, const Expression<Field>& expr) const {
    bool has_shift = expr.shift_ != 0;

    for (auto itr = expr.variables_.begin(); itr != expr.variables_.end();
         ++itr) {
      bool is_last = std::next(itr) == expr.variables_.end();

      if (itr->second == -1) {
        os << "-";
      } else if (itr->second != 1) {
        os << itr->second << "*";
      }

      os << variables_.at(itr->first).first;

      if (!is_last || has_shift) {
        os << " + ";
      }
    }

    if (has_shift || expr.is_constant()) {
      os << expr.shift_;
    }
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

  // transforms all constraints into equalities (by adding slack variables)
  // moves all variables into lhs, all constants into rhs
  // removes all constant constraints
  // throws an exception if any constant constraint is false (e.g. 2 <= 1)
  // removes unused variables
  void normalize() {
    // move all variables to rhs
    for (Constraint<Field>& constraint : constraints_) {
      constraint.lhs_ -= constraint.rhs_;
      constraint.rhs_ = Expression<Field>(0);
      std::swap(constraint.lhs_.shift_, constraint.rhs_.shift_);
      constraint.rhs_.shift_ *= -1;
    }

    // remove constant constraints
    for (auto itr = constraints_.begin(); itr != constraints_.end();) {
      switch (itr->evaluate()) {
        case EvaluationResult::TRUE:
          itr = constraints_.erase(itr);
          break;
        case EvaluationResult::FALSE:
          throw std::runtime_error(
              "Problem is trivially unfeasible. It contains a constant "
              "constraint that evaluates to false.");
          break;
        case EvaluationResult::DONT_KNOW:
          ++itr;
          break;
        default:
          std::unreachable();
      }
    }

    // transform all constraints into equalities
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
      // TODO: slack variable upper bound
      auto slack_var = new_variable(fmt::format("slack({})", slack_index),
                                    VariableType::SLACK, 0, 1'000'000);
      ++slack_index;

      constraint.lhs_ += slack_var;
      constraint.type_ = ConstraintType::EQUAL;
    }

    // remove unused variables
    // TODO: finish this
    // std::unordered_set<size_t> used;
    //
    // for (const Constraint<Field>& constraint : constraints_) {
    //   used.insert_range(constraint.lhs_.variables_ |
    //   std::views::elements<0>); used.insert_range(constraint.rhs_.variables_
    //   | std::views::elements<0>);
    // }
    //
    // for (auto itr = variables_names_.begin(); itr != variables_names_.end();)
    // {
    //   if (itr->)
    // }
  }

  static std::pair<Matrix<Field>, Matrix<Field>> remove_linear_dependent(
      const Matrix<Field>& A, const Matrix<Field>& b) {
    size_t d = A.get_width();
    auto row_basis = linalg::get_row_basis(A);

    Matrix<Field> reduced_A(row_basis.size(), d);
    Matrix<Field> reduced_b(row_basis.size(), 1);

    for (size_t i = 0; i < row_basis.size(); ++i) {
      reduced_A[i, {0, d}] = A[row_basis[i], {0, d}];
      reduced_b[i, 0] = b[row_basis[i], 0];
    }

    return {std::move(reduced_A), std::move(reduced_b)};
  }

 public:
  explicit ProblemBuilder() = default;

  Variable<Field> new_variable(std::string_view name, VariableType type,
                               Field lower_bound, Field upper_bound) {
    size_t index = variables_.size();
    variables_.emplace_back(name, VariableInfo{type, lower_bound, upper_bound});

    return Variable<Field>{index};
  }

  void add_constraint(Constraint<Field> constraint) {
    constraints_.push_back(std::move(constraint));
  }

  void set_objective(Expression<Field> objective) {
    if (objective.shift_ != 0) {
      throw std::invalid_argument("Objective must not contain a shift.");
    }

    objective_ = std::move(objective);
  }

  Field extract_variable(const Matrix<Field>& point,
                         Variable<Field> variable) const {
    return point[variable.id_, 0];
  }

  MILPProblem<Field> get_problem() {
    normalize();

    Matrix<Field> A(constraints_.size(), variables_.size(), 0);
    Matrix<Field> b(constraints_.size(), 1, 0);
    Matrix<Field> c(1, variables_.size(), 0);

    for (const auto& [var, coef] : objective_.variables_) {
      c[0, var] = coef;
    }

    for (size_t i = 0; i < constraints_.size(); ++i) {
      for (const auto& [var, coef] : constraints_[i].lhs_.variables_) {
        A[i, var] = coef;
      }

      b[i, 0] = constraints_[i].rhs_.shift_;
    }

    std::vector<VariableType> types(variables_.size());
    for (size_t i = 0; i < variables_.size(); ++i) {
      types[i] = variables_[i].second.type;
    }

    auto [reduced_A, reduced_b] = remove_linear_dependent(A, b);

    // bounds
    std::vector<Field> lower_bounds(variables_.size());
    std::vector<Field> upper_bounds(variables_.size());

    for (size_t i = 0; i < variables_.size(); ++i) {
      lower_bounds[i] = variables_[i].second.lower_bound;
      upper_bounds[i] = variables_[i].second.upper_bound;
    }

    return {reduced_A, reduced_b, c, types, lower_bounds, upper_bounds};
  }

  //
  friend std::ostream& operator<<(std::ostream& os,
                                  const ProblemBuilder& builder) {
    os << fmt::format("problem builder with {} variables and {} constraints\n",
                      builder.variables_.size(), builder.constraints_.size());

    for (const Constraint<Field>& constraint : builder.constraints_) {
      builder.print_constraint(os, constraint);
      os << "\n";
    }

    return os;
  }
};
