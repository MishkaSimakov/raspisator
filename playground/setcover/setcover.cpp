#include "Evaluator.h"
#include "Greedy.h"
#include "Reader.h"
#include "Types.h"
#include "linear/problem/MILPProblem.h"
#include "linear/problem/ToMatrices.h"
#include "linear/problem/optimization/TransformToEqualities.h"
#include "linear/simplex/Simplex.h"
#include "utils/Paths.h"

auto to_milp(const setcover::Problem& problem) {
  size_t n = problem.elements_count;
  size_t d = problem.sets.size();

  Matrix<double> A(n, d + n, 0);
  Matrix<double> b(n, 1, 0);
  Matrix<double> c(1, d + n, 0);

  for (size_t j = 0; j < n; ++j) {
    for (size_t i = 0; i < d; ++i) {
      if (problem.sets[i].elements.contains(j)) {
        A[j, i] = 1;
      }
    }

    A[j, d + j] = -1;
  }

  for (size_t i = 0; i < n; ++i) {
    b[i, 0] = 1;
  }

  for (size_t i = 0; i < d; ++i) {
    c[0, i] = -static_cast<double>(problem.sets[i].cost);
  }

  Bounds<double> bounds(d + n);
  for (size_t i = 0; i < d; ++i) {
    bounds[i] = {0, 1};
  }

  for (size_t i = d; i < d + n; ++i) {
    bounds[i] = {0, std::nullopt};
  }

  return std::make_tuple(A, b, c, bounds);
}

// Задача sc_330_0 вызывает значительные вычислительные трудности у текущей
// версии. Скорее всего это происходит из-за ошибок округления в LU или
// симплексе. Интересно, что в конце определитель матрицы улетает в
// бесконечность.
int main() {
  auto problem = setcover::read_problem(paths::resource("setcover/sc_10000_5"));

  auto greedy_solution = setcover::Greedy().solve(problem);
  auto greedy_result = setcover::evaluate(problem, greedy_solution);

  std::cout << greedy_result.score << std::endl;

  auto [A, b, c, bounds] = to_milp(problem);

  auto [n, d] = A.shape();
  std::println("{} x {}", n, d);

  // slightly perturb bounds
  std::default_random_engine engine;
  std::uniform_int_distribution distr(0, 5);
  for (size_t i = 0; i < d; ++i) {
    if (bounds[i].lower) {
      *bounds[i].lower -= static_cast<double>(distr(engine)) * 1e-6;
    }

    if (bounds[i].upper) {
      *bounds[i].upper += static_cast<double>(distr(engine)) * 1e-6;
    }
  }

  // run simplex algorithm
  simplex::Settings<double> settings{.is_strict = true};
  auto simplex = simplex::Simplex<double, simplex::LoggingAccountant<double>>(
      CSCMatrix(A), b, c, settings);

  // construct feasible point
  std::vector states(d, VariableState::AT_UPPER);
  for (size_t i = 0; i < n; ++i) {
    states[d - i - 1] = VariableState::BASIC;
  }

  auto result = simplex.primal(bounds, states);

  std::cout << std::get<FiniteLPSolution<double>>(result.solution).value
            << std::endl;
}
