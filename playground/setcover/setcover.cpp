#include "Reader.h"
#include "Types.h"
#include "linear/problem/MILPProblem.h"
#include "linear/problem/ToMatrices.h"
#include "linear/problem/optimization/TransformToEqualities.h"
#include "linear/simplex/Simplex.h"
#include "utils/Paths.h"

MILPProblem<double> to_milp(const setcover::Problem& problem) {
  MILPProblem<double> result;

  std::vector<Variable<double>> variables;
  for (size_t i = 0; i < problem.sets.size(); ++i) {
    variables.push_back(result.new_variable(std::format("x_{}", i),
                                            VariableType::INTEGER, 0, 1));
  }

  for (size_t j = 0; j < problem.elements_count; ++j) {
    Expression<double> constraint;

    for (size_t i = 0; i < problem.sets.size(); ++i) {
      if (problem.sets[i].elements.contains(j)) {
        constraint += variables[i];
      }
    }

    result.add_constraint(constraint >= Expression<double>{1});
  }

  Expression<double> objective;
  for (size_t i = 0; i < problem.sets.size(); ++i) {
    objective += static_cast<double>(problem.sets[i].cost) * variables[i];
  }
  result.set_objective(-objective);

  return result;
}

// Задача sc_330_0 вызывает значительные вычислительные трудности у текущей
// версии. Скорее всего это происходит из-за ошибок округления в LU или
// симплексе. Интересно, что в конце определитель матрицы улетает в
// бесконечность.
int main() {
  auto problem = setcover::read_problem(paths::resource("setcover/sc_10000_2"));

  auto milp_problem = to_milp(problem);

  auto optimizer = TransformToEqualities<double>();
  auto optimized = optimizer.apply(milp_problem);
  auto matrices = to_matrices(optimized);

  auto [n, d] = matrices.A.shape();
  std::println("{} x {}", n, d);

  // slightly perturb bounds
  std::default_random_engine engine;
  std::uniform_int_distribution distr(0, 5);
  for (size_t i = 0; i < d; ++i) {
    if (matrices.bounds[i].lower) {
      *matrices.bounds[i].lower -= static_cast<double>(distr(engine)) * 1e-5;
    }

    if (matrices.bounds[i].upper) {
      *matrices.bounds[i].upper += static_cast<double>(distr(engine)) * 1e-5;
    }
  }

  // run simplex algorithm
  simplex::Settings<double> settings{.is_strict = true};
  auto simplex = simplex::Simplex<double, simplex::LoggingAccountant<double>>(
      CSCMatrix(matrices.A), matrices.b, matrices.c, settings);

  // construct feasible point
  std::vector states(d, VariableState::AT_UPPER);
  for (size_t i = 0; i < n; ++i) {
    states[d - i - 1] = VariableState::BASIC;
  }

  auto result = simplex.primal(matrices.bounds, states);

  std::cout << std::get<FiniteLPSolution<double>>(result.solution).value
            << std::endl;
}
