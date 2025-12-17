#include <random>

#include "linear/BigInteger.h"
#include "linear/model/LP.h"
#include "linear/simplex/BoundedSimplexMethod.h"

using Field = double;

#include "simplex_core_dump_0.h"

using namespace SimplexDump_0;

Matrix<Field> perturbed_objective(Matrix<Field> c,
                                  const std::vector<Field>& lower,
                                  const std::vector<Field>& upper, size_t n) {
  std::vector<size_t> perturbed;

  for (size_t i = 0; i < n; ++i) {
    if (FieldTraits<Field>::is_nonzero(c[0, i])) {
      continue;
    }

    if (FieldTraits<Field>::is_nonzero(lower[i])) {
      continue;
    }

    if (!FieldTraits<Field>::is_nonzero(lower[i] - upper[i])) {
      continue;
    }

    // remove variables with big upper bounds so that coefficients are not
    // too small
    // if (upper[i] > 500) {
    // continue;
    // }

    c[0, i] -= static_cast<Field>(1) / upper[i];
    perturbed.push_back(i);
  }

  for (size_t i : perturbed) {
    c[0, i] /= static_cast<Field>(2 * perturbed.size());
  }

  return c;
}

std::default_random_engine random_engine;

auto perturb_bounds(const std::vector<Field>& lower,
                    const std::vector<Field>& upper) {
  std::vector<Field> new_lower = lower;
  std::vector<Field> new_upper = upper;

  std::uniform_int_distribution<int> shift_distribution(0, 3);

  for (size_t i = 0; i < new_lower.size(); ++i) {
    new_lower[i] -= Field(shift_distribution(random_engine));
    new_upper[i] += Field(shift_distribution(random_engine));
  }

  return std::pair{new_lower, new_upper};
}

int main() {
  auto basic_variables = linalg::get_row_basis(linalg::transposed(A));

  simplex::Settings<Field> settings = {.max_iterations = 100'000};

  std::ofstream os("iterations.csv");
  os << "count\n";

  for (size_t i = 0; i < c.get_width(); ++i) {
    auto perturbed = perturbed_objective(c, lower, upper, i);

    auto solver =
        simplex::BoundedSimplexMethod(CSCMatrix(A), b, perturbed, settings);
    solver.setup_warm_start(basic_variables);
    auto run_result = solver.dual(lower, upper);

    os << run_result.iterations_count << '\n';

    std::cout << run_result.iterations_count << std::endl;
  }
}
