#pragma once

#include "linalg/Matrix.h"
#include "linalg/Permutation.h"
#include "linalg/lu/EtaFile.h"
#include "linear/FieldTraits.h"
#include "utils/Accumulators.h"

namespace linalg {

template <typename Field>
Field scale_factor(const Matrix<Field>&) {
  return 1;
}

template <>
inline double scale_factor(const Matrix<double>& b) {
  auto [n, d] = b.shape();

  Maximum<double> max;
  Minimum<double> min;

  for (size_t i = 0; i < d; ++i) {
    for (size_t j = 0; j < n; ++j) {
      max.record(FieldTraits<double>::abs(b[j, i]));

      if (FieldTraits<double>::is_nonzero(b[j, i])) {
        min.record(FieldTraits<double>::abs(b[j, i]));
      }
    }
  }

  if (!min.has_value() || !max.has_value()) {
    return 1;
  }

  const double scale_factor = std::sqrt(*max * *min);

  // round to power of 2
  const int power = std::round(std::log2(scale_factor));

  return std::exp2(-power);
}

// solves Ax = b, where PAQ = LU
// ls is eta-file containing L1^-1 ... Ln^-1, where L = L1 ... Ln
// us is eta-file containing U1^-1 ... Un^-1, where U = Un ... U1
template <typename Field>
Vector<Field> solve_linear(Vector<Field> b, const Permutation& P,
                           const Permutation& Q, const EtaFile<Field>& ls,
                           const EtaFile<Field> us) {
  Vector<Field> result = P.apply(std::move(b));

  auto sf = scale_factor(result);
  result /= sf;

  for (auto entry : ls) {
    result = entry.apply(std::move(result));
  }

  for (auto entry : us | std::views::reverse) {
    result = entry.apply(std::move(result));
  }

  result = Q.apply(std::move(result));
  result *= sf;

  return result;
}

// solves A^T x = b, where PAQ = LU
// ls is eta-file containing L1^-1 ... Ln^-1, where L = L1 ... Ln
// us is eta-file containing U1^-1 ... Un^-1, where U = Un ... U1
template <typename Field>
Vector<Field> solve_linear_transposed(Vector<Field> b, const Permutation& P,
                                      const Permutation& Q,
                                      const EtaFile<Field>& ls,
                                      const EtaFile<Field>& us) {
  b = Q.apply_transposed(std::move(b));

  auto sf = scale_factor(b);
  b /= sf;

  for (auto entry : us) {
    b = entry.apply_transposed(std::move(b));
  }

  for (auto entry : ls | std::views::reverse) {
    b = entry.apply_transposed(std::move(b));
  }

  b = P.apply_transposed(std::move(b));
  b *= sf;

  return b;
}

}  // namespace linalg
