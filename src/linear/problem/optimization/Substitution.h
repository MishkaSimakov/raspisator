#pragma once

#include <cassert>
#include <unordered_map>

#include "BaseOptimizer.h"

// TODO: check bounds
template <typename Field>
class Substitution final : public BaseOptimizer<Field> {
  // stores substitutions
  // (x, y, alpha, beta) means that variable x was replaced with
  // y * alpha + beta
  std::vector<std::tuple<size_t, size_t, Field, Field>> substitutions_;

  void substitute_in_expression(Expression<Field>& expr,
                                const std::string& variable,
                                Expression<Field> replacement) {
    auto itr = expr.get_variables().find(variable);

    if (itr != expr.get_variables().end()) {
      expr -= itr->second * replacement;
    }
  }

  void substitute(MILPProblem<Field>& problem, const std::string& variable,
                  Expression<Field> replacement) {
    for (auto& constraint : problem.constraints) {
      substitute_in_expression(constraint.expr, variable, replacement);
    }

    substitute_in_expression(problem.objective, variable, replacement);
  }

  size_t occurrences_count(const MILPProblem<Field>& problem,
                           const std::string& variable) {
    size_t result = 0;

    for (const auto& constraint : problem.constraints) {
      if (constraint.expr.get_variables().contains(variable)) {
        ++result;
      }
    }

    return result;
  }

  void erase_variable(MILPProblem<Field>& problem,
                      const std::string& variable) {
    for (auto itr = problem.variables.begin(); itr != problem.variables.end();
         ++itr) {
      if (itr->name == variable) {
        problem.variables.erase(itr);
        return;
      }
    }
  }

  bool was_substituted(size_t index) const {
    for (auto x : substitutions_ | std::views::elements<0>) {
      if (x == index) {
        return true;
      }
    }

    return false;
  }

 public:
  Substitution() = default;

  MILPProblem<Field> apply(MILPProblem<Field> problem) override {
    auto enumeration = problem.enumerate_variables();

    for (auto itr = problem.constraints.begin();
         itr != problem.constraints.end();) {
      if (itr->type != ConstraintType::EQUAL_ZERO) {
        ++itr;
        continue;
      }

      auto& variables = itr->expr.get_variables();

      if (variables.size() == 2) {
        std::pair<std::string, Field> a = *variables.begin();
        std::pair<std::string, Field> b = *std::next(variables.begin());

        if (occurrences_count(problem, a.first) >
            occurrences_count(problem, b.first)) {
          std::swap(a, b);
        }

        // replace a with b
        auto substitutee = problem.get_variable(a.first);
        auto replacement = itr->expr.express(substitutee);

        assert(replacement.has_value());

        substitute(problem, a.first, *replacement);

        // store substitution
        substitutions_.emplace_back(enumeration[a.first], enumeration[b.first],
                                    -b.second / a.second,
                                    -itr->expr.get_shift() / a.second);

        erase_variable(problem, a.first);
        itr = problem.constraints.erase(itr);
      } else {
        ++itr;
      }
    }

    return problem;
  }

  Matrix<Field> inverse(const Matrix<Field>& point) override {
    size_t d = point.get_height();

    Matrix<Field> result(d + substitutions_.size(), 1);

    size_t j = 0;
    for (size_t i = 0; i < d + substitutions_.size(); ++i) {
      if (!was_substituted(i)) {
        result[i, 0] = point[j, 0];
        ++j;
      }
    }

    for (auto [x, y, alpha, beta] : substitutions_ | std::views::reverse) {
      result[x, 0] = result[y, 0] * alpha + beta;
    }

    return result;
  }
};
