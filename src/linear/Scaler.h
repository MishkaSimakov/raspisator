#pragma once

#include <cmath>

#include "model/MILP.h"
#include "utils/Accumulators.h"

template <typename Field>
class Scaler {
  static_assert(std::is_floating_point_v<Field>,
                "Currently works only for standart types.");

  const MILPProblem<Field>& problem_;

  static double get_coefficient(MatrixSlice<const Field> slice) {
    GeometricAverage<Field> coef;

    auto [n, d] = slice.shape();

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < d; ++j) {
        if (FieldTraits<Field>::is_nonzero(slice[i, j])) {
          coef.record(FieldTraits<Field>::abs(slice[i, j]));
        }
      }
    }

    double power = std::round(std::log2(coef.average()));
    double rounded_coef = std::exp2(-power);

    return rounded_coef;
  }

 public:
  explicit Scaler(const MILPProblem<Field>& problem) : problem_(problem) {}

  // scales only using powers of 2
  MILPProblem<Field> get_scaled() const {
    auto [n, d] = problem_.A.shape();

    auto result = problem_;

    // scale rows
    for (size_t i = 0; i < n; ++i) {
      auto coef = get_coefficient(problem_.A[i, {0, d}]);

      for (size_t j = 0; j < d; ++j) {
        result.A[i, j] *= coef;
      }
      result.b[i, 0] *= coef;
    }

    // scale columns
    for (size_t j = 0; j < d; ++j) {
      auto coef = get_coefficient(problem_.A[{0, n}, j]);

      for (size_t i = 0; i < n; ++i) {
        result.A[i, j] *= coef;
      }
      result.c[0, j] *= coef;

      result.lower_bounds[j] /= coef;
      result.upper_bounds[j] /= coef;
    }

    return result;
  }

  // Well scaled if returned value < 2
  // https://pure.iiasa.ac.at/id/eprint/4172/7/WP-94-037.pdf
  double get_scaling_quality() const {
    Minimum<Field> min;
    Maximum<Field> max;

    auto [n, d] = problem_.A.shape();

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < d; ++j) {
        if (FieldTraits<Field>::is_nonzero(problem_.A[i, j])) {
          Field abs = FieldTraits<Field>::abs(problem_.A[i, j]);

          min.record(abs);
          max.record(abs);
        }
      }
    }

    return std::log10(*max.max() / *min.min());
  }
};
