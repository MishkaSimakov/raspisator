#pragma once

#include <vector>

#include "presolve/Pass.h"
#include "utils/Accumulators.h"

namespace presolve {

template <typename Field>
class TransformToEqualities final : public Pass<Field> {
  size_t initial_var_count_;
  size_t added_slack_count_{0};

 public:
  TransformToEqualities() = default;

  problem::MILP<Field> apply(problem::MILP<Field> problem) override {
    this->register_apply();

    // for each non-zero rhs range a new slack variable is added
    const auto [n, d] = problem.matrix.shape();

    initial_var_count_ = d;

    for (size_t row = 0; row < n; ++row) {
      if (!problem.rhs_bounds[row].is_fixed()) {
        problem.matrix.add_column();
        problem.matrix.push_to_last_column(row, -1);

        const auto slack_bound = problem.rhs_bounds[row];

        problem.var_bounds.push_back(slack_bound);
        problem.implied_var_bounds.push_back(slack_bound);
        problem.rhs_bounds[row] = Bound<Field>{0, 0};

        problem.var_names.push_back(problem.row_name(row) + "_range");

        ++added_slack_count_;
      }
    }

    problem.is_integer.resize(d + added_slack_count_, false);
    problem.implied_is_integer.resize(d + added_slack_count_, false);
    problem.cost.resize(d + added_slack_count_);

    return problem;
  }

  Vector<Field> inverse(Vector<Field> solution) const override {
    if (solution.size() != initial_var_count_ + added_slack_count_) {
      throw std::invalid_argument("Wrong solution size in inverse.");
    }

    solution.resize(initial_var_count_);

    return solution;
  }
};

}  // namespace presolve
