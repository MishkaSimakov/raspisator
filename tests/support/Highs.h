#pragma once

#include <Highs.h>
#include <filesystem>

#include "presolve/Pass.h"
#include "problem/MILP.h"

namespace highs {

struct Solution {
  HighsModelStatus status;
  Vector<double> x;
  double objective;
};

inline Solution inverse_pass(const problem::MILP<double>& problem,
                             const presolve::Pass<double>& pass,
                             Solution solution) {
  if (solution.status != HighsModelStatus::kOptimal) {
    return solution;
  }

  solution.x = pass.inverse(solution.x);
  solution.objective = linalg::dot(problem.cost, solution.x);

  return solution;
}

inline HighsLp read_mps(const std::filesystem::path& path) {
  Highs highs;
  highs.setOptionValue("output_flag", false);

  if (highs.readModel(path) != HighsStatus::kOk) {
    throw std::runtime_error("HiGHS failed to read: " + path.string());
  }

  return highs.getLp();
}

template <typename Field>
HighsLp from_milp(const problem::MILP<Field>& problem) {
  const auto [n, d] = problem.matrix.shape();

  HighsLp lp;

  lp.model_name_ = problem.name;
  lp.objective_name_ = problem.cost_name;

  lp.num_col_ = static_cast<int>(d);
  lp.num_row_ = static_cast<int>(n);

  lp.sense_ = ObjSense::kMaximize;

  // Objective
  lp.col_cost_.resize(d);
  for (size_t i = 0; i < d; ++i) {
    lp.col_cost_[i] = static_cast<double>(problem.cost[i]);
  }

  lp.offset_ = static_cast<double>(problem.cost_offset);

  // Bounds
  lp.col_lower_.resize(d);
  lp.col_upper_.resize(d);
  lp.col_names_.resize(d);

  for (size_t i = 0; i < d; ++i) {
    lp.col_lower_[i] =
        problem.var_bounds[i]
            .lower
            .transform([](Field value) { return static_cast<double>(value); })
            .value_or(-kHighsInf);

    lp.col_upper_[i] =
        problem.var_bounds[i]
            .upper
            .transform([](Field value) { return static_cast<double>(value); })
            .value_or(kHighsInf);

    lp.col_names_[i] = problem.var_names[i];
  }

  // Matrix
  lp.row_lower_.resize(n);
  lp.row_upper_.resize(n);
  lp.row_names_.resize(n);

  for (size_t row = 0; row < n; ++row) {
    lp.row_lower_[row] =
        problem.rhs_bounds[row]
            .lower
            .transform([](Field value) { return static_cast<double>(value); })
            .value_or(-kHighsInf);

    lp.row_upper_[row] =
        problem.rhs_bounds[row]
            .upper
            .transform([](Field value) { return static_cast<double>(value); })
            .value_or(kHighsInf);

    lp.row_names_[row] = problem.row_names[row];
  }

  lp.a_matrix_.start_.resize(d + 1);
  lp.a_matrix_.index_.clear();
  lp.a_matrix_.value_.clear();

  size_t nnz = 0;
  lp.a_matrix_.start_[0] = 0;

  for (size_t col = 0; col < d; ++col) {
    for (auto [row, value] : problem.matrix.get_column(col)) {
      lp.a_matrix_.index_.push_back(static_cast<int>(row));
      lp.a_matrix_.value_.push_back(static_cast<double>(value));
      ++nnz;
    }

    lp.a_matrix_.start_[col + 1] = static_cast<int>(nnz);
  }

  lp.integrality_.resize(d);
  for (size_t j = 0; j < d; ++j) {
    lp.integrality_[j] = problem.is_integer[j] ? HighsVarType::kInteger
                                               : HighsVarType::kContinuous;
  }

  return lp;
}

inline Solution solve(const HighsLp& problem) {
  Highs highs;

  if (highs.setOptionValue("output_flag", false) != HighsStatus::kOk) {
    throw std::runtime_error("Error during highs::setOptionValue call.");
  }

  if (highs.passModel(problem) != HighsStatus::kOk) {
    throw std::runtime_error("Error during highs::passModel call.");
  }

  if (highs.run() != HighsStatus::kOk) {
    throw std::runtime_error("Error during highs::run call.");
  }

  Solution solution;
  solution.status = highs.getModelStatus();

  const auto& highs_sol = highs.getSolution();

  solution.x = Vector<double>(highs_sol.col_value.size());
  for (size_t i = 0; i < highs_sol.col_value.size(); ++i) {
    solution.x[i] = highs_sol.col_value[i];
  }
  solution.objective = highs.getInfo().objective_function_value;

  return solution;
}

}  // namespace highs
