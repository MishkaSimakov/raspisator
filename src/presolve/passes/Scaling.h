#pragma once

#include <cmath>

#include "presolve/Pass.h"

namespace presolve {

template <typename Field>
class Scaling final : public Pass<Field> {
  std::vector<Field> variables_scale_factors_;

  static Field round_scale_factor(double scale_factor) {
    int power = std::round(std::log2(scale_factor));
    return FieldTraits<Field>::exp2(-power);
  }

  static std::vector<Field> get_rows_scale_factor(
      const CSCMatrix<Field>& matrix) {
    std::vector<GeometricMean<double>> scale_factors(matrix.rows());

    for (size_t col = 0; col < matrix.cols(); ++col) {
      for (const auto& [row, value] : matrix.get_column(col)) {
        if (abs(value) > FieldTraits<Field>::tolerance) {
          scale_factors[row].record(static_cast<double>(abs(value)));
        }
      }
    }

    std::vector<Field> result(matrix.rows());

    for (size_t row = 0; row < matrix.rows(); ++row) {
      if (!scale_factors[row].has_value()) {
        result[row] = 1;
      } else {
        result[row] = round_scale_factor(*scale_factors[row]);
      }
    }

    return result;
  }

  static Field get_col_scale_factor(const CSCMatrix<Field>& matrix,
                                    size_t col) {
    using std::abs;

    GeometricMean<double> scale_factor;

    for (const auto& [row, value] : matrix.get_column(col)) {
      if (abs(value) > FieldTraits<Field>::tolerance) {
        scale_factor.record(static_cast<double>(abs(value)));
      }
    }

    return !scale_factor.has_value() ? Field(1)
                                     : round_scale_factor(*scale_factor);
  }

  static Field get_cost_scale_factor(const Vector<Field>& cost) {
    using std::abs;

    GeometricMean<double> scale_factor;

    for (size_t i = 0; i < cost.size(); ++i) {
      if (abs(cost[i]) > FieldTraits<Field>::tolerance) {
        scale_factor.record(static_cast<double>(abs(cost[i])));
      }
    }

    return !scale_factor.has_value() ? Field(1)
                                     : round_scale_factor(*scale_factor);
  }

 public:
  Scaling() = default;

  problem::MILP<Field> apply(problem::MILP<Field> problem) override {
    this->register_apply();

    variables_scale_factors_.resize(problem.matrix.cols());

    // scale rows
    auto scale_factors = get_rows_scale_factor(problem.matrix);
    for (size_t col = 0; col < problem.matrix.cols(); ++col) {
      for (auto& [row, value] : problem.matrix.get_column(col)) {
        value *= scale_factors[row];
      }
    }
    for (size_t row = 0; row < problem.matrix.rows(); ++row) {
      problem.rhs_bounds[row] *= scale_factors[row];
    }

    // scale columns
    for (size_t col = 0; col < problem.matrix.cols(); ++col) {
      if (problem.is_integer[col]) {
        variables_scale_factors_[col] = 1;
        continue;
      }

      Field scale_factor = get_col_scale_factor(problem.matrix, col);

      variables_scale_factors_[col] = scale_factor;

      for (auto& [row, value] : problem.matrix.get_column(col)) {
        value *= scale_factor;
      }

      problem.var_bounds[col] /= scale_factor;
      problem.implied_var_bounds[col] /= scale_factor;
      problem.cost[col] *= scale_factor;

      problem.implied_is_integer[col] = false;
    }

    // scale cost
    problem.cost *= get_cost_scale_factor(problem.cost);

    return problem;
  }

  // Well scaled if returned value < 2
  // https://pure.iiasa.ac.at/id/eprint/4172/7/WP-94-037.pdf
  double get_scaling_quality(const problem::MILP<Field>& problem) const {
    using std::abs;

    Minimum<Field> min;
    Maximum<Field> max;

    for (size_t col = 0; col < problem.matrix.cols(); ++col) {
      for (const auto& [row, value] : problem.matrix.get_column(col)) {
        min.record(abs(value));
        max.record(abs(value));
      }
    }

    return std::log10(static_cast<double>(*max / *min));
  }

  Vector<Field> inverse(Vector<Field> solution) const override {
    for (size_t i = 0; i < solution.size(); ++i) {
      solution[i] *= variables_scale_factors_[i];
    }

    return solution;
  }
};

}  // namespace presolve
