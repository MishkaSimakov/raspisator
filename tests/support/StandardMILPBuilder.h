#pragma once

#include <random>
#include <vector>

#include "linalg/Matrix.h"
#include "linalg/Random.h"
#include "linalg/Stack.h"
#include "problem/Bound.h"
#include "problem/StandardMILP.h"

// Builder for problem::StandardMILP with guaranteed structural properties.
// All construction is algebraic — no external solver required.
//
// Usage:
//   auto result =
//   StandardMILPBuilder<Rational>{}.rows(3).cols(6).build_feasible(); auto
//   problem  = StandardMILPBuilder<Rational>{}.rows(3).build_infeasible(); auto
//   problem  = StandardMILPBuilder<Rational>{}.rows(3).build_unbounded();
template <typename Field>
class StandardMILPBuilder {
  size_t rows_ = 3;
  size_t cols_ = 6;  // must be > rows_
  int magnitude_ = 10;
  uint64_t seed_ = 0;

 public:
  StandardMILPBuilder& rows(size_t r) {
    rows_ = r;
    return *this;
  }

  StandardMILPBuilder& cols(size_t c) {
    cols_ = c;
    return *this;
  }

  StandardMILPBuilder& magnitude(int m) {
    magnitude_ = m;
    return *this;
  }

  StandardMILPBuilder& seed(uint64_t s) {
    seed_ = s;
    return *this;
  }

  struct FeasibleResult {
    problem::StandardMILP<Field> problem;
    Vector<Field> witness;  // x* used in construction
    std::vector<simplex::VariableState>
        primal_states;  // valid initial primal states
  };

  // Generates a StandardMILP that is guaranteed to be primal feasible.
  //
  // Construction:
  //   - A = [A_basic | A_nonbasic] with A_basic invertible
  //   - x*_basic chosen randomly within bounds, x*_nonbasic = 0 (AT_LOWER)
  //   - rhs = A * x*
  //   - Primal states: first rows_ variables = BASIC, rest = AT_LOWER
  FeasibleResult build_feasible() {
    const size_t n = rows_;
    const size_t d = cols_;
    if (d <= n) {
      throw std::invalid_argument("cols must be > rows for a feasible LP");
    }

    std::default_random_engine rng(seed_);
    std::uniform_int_distribution<int> dist(-magnitude_, magnitude_);
    // Nonzero random integer — avoid 0 so cost is non-degenerate
    std::uniform_int_distribution<int> nonzero_dist(1, magnitude_);

    // Generate invertible n×n basic submatrix
    auto A_basic = linalg::random::dense_invertible<Field>(n, rng, dist);
    // Generate random n×(d-n) nonbasic submatrix
    auto A_nonbasic = linalg::random::dense<Field>(n, d - n, rng, dist);
    auto A = linalg::hstack(A_basic, A_nonbasic);

    // Choose x*_basic randomly
    linalg::Vector<Field> x_basic(n);
    for (size_t i = 0; i < n; ++i) {
      x_basic[i] = Field(dist(rng));
    }

    // x*_nonbasic = 0 (will be AT_LOWER at lower bound 0)
    // rhs = A_basic * x*_basic + A_nonbasic * 0 = A_basic * x*_basic
    linalg::Vector<Field> rhs(n);
    for (size_t row = 0; row < n; ++row) {
      Field val = 0;
      for (size_t col = 0; col < n; ++col) {
        val += A_basic[row, col] * x_basic[col];
      }
      rhs[row] = val;
    }

    // Set up bounds
    std::vector<Bound<Field>> bounds(d);
    for (size_t i = 0; i < n; ++i) {
      // Basic variable bounds contain x*_basic[i]
      bounds[i] = Bound<Field>{x_basic[i] - Field(magnitude_),
                               x_basic[i] + Field(magnitude_)};
    }
    for (size_t j = n; j < d; ++j) {
      // Nonbasic variables: lower bound 0, upper bound magnitude_
      bounds[j] = Bound<Field>{Field(0), Field(magnitude_)};
    }

    // Build witness vector
    linalg::Vector<Field> witness(d);
    for (size_t i = 0; i < n; ++i) {
      witness[i] = x_basic[i];
    }
    for (size_t j = n; j < d; ++j) {
      witness[j] = Field(0);
    }

    // Primal states: first n = BASIC, rest = AT_LOWER
    std::vector<simplex::VariableState> states(d);
    for (size_t i = 0; i < n; ++i) {
      states[i] = simplex::VariableState::BASIC;
    }
    for (size_t j = n; j < d; ++j) {
      states[j] = simplex::VariableState::AT_LOWER;
    }

    // Build cost (random, nonzero for non-degeneracy)
    linalg::Vector<Field> cost(d);
    for (size_t j = 0; j < d; ++j) {
      cost[j] = Field(dist(rng));
    }

    problem::StandardMILP<Field> p;
    p.matrix = linalg::CSCMatrix<Field>(A);
    p.rhs = rhs;
    p.cost = cost;
    p.var_bounds = bounds;
    p.var_names.resize(d);
    p.row_names.resize(n);
    p.is_integer.assign(d, false);

    return {std::move(p), std::move(witness), std::move(states)};
  }

  // Generates a StandardMILP that is guaranteed to be infeasible.
  //
  // Construction:
  //   - Random A (n×d) with finite bounds [0, mag]
  //   - For row 0: compute minimum achievable row value over the box,
  //     then set rhs[0] = that minimum − 1 (strictly outside range)
  //   - For remaining rows: rhs = mid-point of achievable range
  problem::StandardMILP<Field> build_infeasible() {
    const size_t n = rows_;
    const size_t d = cols_;

    std::default_random_engine rng(seed_);
    std::uniform_int_distribution<int> dist(-magnitude_, magnitude_);

    auto A = linalg::random::dense<Field>(n, d, rng, dist);

    std::vector<Bound<Field>> bounds(d);
    for (size_t j = 0; j < d; ++j) {
      bounds[j] = Bound<Field>{Field(0), Field(magnitude_)};
    }

    linalg::Vector<Field> rhs(n);

    // Row 0: set rhs strictly below achievable minimum
    {
      Field row_min = 0;
      for (size_t j = 0; j < d; ++j) {
        Field a = A[0, j];
        if (a > Field(0)) {
          row_min += a * *bounds[j].lower;
        } else {
          row_min += a * *bounds[j].upper;
        }
      }
      rhs[0] = row_min - Field(1);
    }

    // Remaining rows: achievable midpoint
    for (size_t i = 1; i < n; ++i) {
      Field row_min = 0;
      Field row_max = 0;
      for (size_t j = 0; j < d; ++j) {
        Field a = A[i, j];
        if (a > Field(0)) {
          row_min += a * *bounds[j].lower;
          row_max += a * *bounds[j].upper;
        } else {
          row_min += a * *bounds[j].upper;
          row_max += a * *bounds[j].lower;
        }
      }
      rhs[i] = (row_min + row_max) / Field(2);
    }

    linalg::Vector<Field> cost(d);
    for (size_t j = 0; j < d; ++j) {
      cost[j] = Field(dist(rng));
    }

    problem::StandardMILP<Field> p;
    p.matrix = linalg::CSCMatrix<Field>(A);
    p.rhs = rhs;
    p.cost = cost;
    p.var_bounds = bounds;
    p.var_names.resize(d);
    p.row_names.resize(n);
    p.is_integer.assign(d, false);

    return p;
  }

  // Generates a StandardMILP that is guaranteed to be unbounded.
  //
  // Construction:
  //   - d = n+1 (one extra "free" variable)
  //   - A = [A_basic | 0] — last column is zeros, so last variable is
  //     unconstrained by the equalities
  //   - cost = [0,...,0,1] — maximize the unconstrained variable
  //   - bounds: first n vars in [0, mag], last var in [0, +∞)
  //   - rhs = 0 (achieved by x* = 0, basic vars = first n)
  //   - Primal states: first n = BASIC, last = AT_LOWER
  struct UnboundedResult {
    problem::StandardMILP<Field> problem;
    std::vector<simplex::VariableState> primal_states;
  };

  UnboundedResult build_unbounded() {
    const size_t n = rows_;
    const size_t d = n + 1;

    std::default_random_engine rng(seed_);
    std::uniform_int_distribution<int> dist(-magnitude_, magnitude_);

    // A_basic: invertible n×n
    auto A_basic = linalg::random::dense_invertible<Field>(n, rng, dist);

    // Extend to n×(n+1) with a zero column
    Matrix<Field> A(n, d);
    for (size_t row = 0; row < n; ++row) {
      for (size_t col = 0; col < n; ++col) {
        A[row, col] = A_basic[row, col];
      }
      // A[row, n] = 0 (already zeroed)
    }

    // rhs = 0 (x* = 0 with first n as basis gives A_basic * 0 = 0)
    Vector<Field> rhs(n);
    for (size_t i = 0; i < n; ++i) {
      rhs[i] = Field(0);
    }

    // cost: maximize last variable
    Vector<Field> cost(d);
    cost[d - 1] = Field(1);

    // Bounds
    std::vector<Bound<Field>> bounds(d);
    for (size_t j = 0; j < n; ++j) {
      bounds[j] = Bound<Field>{Field(0), Field(magnitude_)};
    }
    // Last variable: lower bound only (unbounded above)
    bounds[d - 1] = Bound<Field>{Field(0), std::nullopt};

    // Primal states: first n = BASIC (at value 0), last = AT_LOWER (at value 0)
    std::vector<simplex::VariableState> states(d);
    for (size_t i = 0; i < n; ++i) {
      states[i] = simplex::VariableState::BASIC;
    }
    states[d - 1] = simplex::VariableState::AT_LOWER;

    problem::StandardMILP<Field> p;
    p.matrix = linalg::CSCMatrix<Field>(A);
    p.rhs = rhs;
    p.cost = cost;
    p.var_bounds = bounds;
    p.var_names.resize(d);
    p.row_names.resize(n);
    p.is_integer.assign(d, false);

    return {std::move(p), std::move(states)};
  }
};
