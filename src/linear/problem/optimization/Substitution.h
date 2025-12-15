// TODO
// #pragma once
//
// #include <unordered_map>
//
// #include "BaseOptimizer.h"
//
// template <typename Field>
// class Substitution final : public BaseOptimizer<Field> {
//   std::unordered_map<size_t, Expression<Field>> substitutions_;
//
//   void substitute_in_expression(Expression<Field>& expr,
//                                 const std::string& variable,
//                                 Expression<Field> replacement) {
//     auto itr = expr.get_variables().find(variable);
//
//     if (itr != expr.get_variables().end()) {
//       expr.get_variables().erase(itr);
//       expr += replacement;
//     }
//   }
//
//   void substitute(MILPProblem<Field>& problem, Expression<Field> replacement) {
//     for (auto& constraint : problem.constraints) {
//       substitute_in_expression(constraint.lhs, replacement);
//
//       constraint.rhs -= constraint.lhs.shift_;
//       constraint.lhs.shift_ = 0;
//     }
//
//     substitute_in_expression(problem.objective, variable, replacement);
//   }
//
//  public:
//   Substitution() = default;
//
//   MILPProblem<Field> apply(MILPProblem<Field> problem) override {
//     // TODO: one must choose variable wisely
//     for (auto itr = problem.constraints.begin();
//          itr != problem.constraints.end();) {
//       auto& variables = itr->lhs.get_variables();
//
//       if (variables.size() == 2) {
//         auto substitutee = problem.get_variable(*variables.begin());
//         auto replacement = itr->expr.express(substitutee);
//
//         substitute(problem, replacement);
//       }
//     }
//   }
//
//   Matrix<Field> inverse(const Matrix<Field>& point) override;
// };
