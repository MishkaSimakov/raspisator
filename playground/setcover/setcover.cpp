#include <random>

#include "Evaluator.h"
#include "Greedy.h"
#include "Reader.h"
#include "Types.h"
#include "presolve/Presolve.h"
#include "problem/MILP.h"
#include "problem/StandardLP.h"
#include "simplex/Simplex.h"
#include "utils/Paths.h"

auto to_milp(const setcover::Problem& problem) {
  problem::MILP<double> result;

  size_t n = problem.elements_count;
  size_t d = problem.sets.size();

  // constraints
  result.matrix.resize(n, 0);

  for (size_t col = 0; col < d; ++col) {
    result.matrix.add_column();

    for (size_t row = 0; row < n; ++row) {
      if (problem.sets[col].elements.contains(row)) {
        result.matrix.push_to_last_column(row, 1);
      }
    }
  }

  // rhs
  result.rhs_bounds.resize(n);
  result.row_names.resize(n);

  for (size_t i = 0; i < n; ++i) {
    result.rhs_bounds[i] = Bound<double>{1, std::nullopt};
  }

  // cost
  result.cost.resize(d);
  for (size_t i = 0; i < d; ++i) {
    result.cost[i] = -static_cast<double>(problem.sets[i].cost);
  }

  // variables
  result.var_bounds.resize(d);
  result.var_names.resize(d);

  for (size_t i = 0; i < d; ++i) {
    result.var_bounds[i] = {0, 1};
  }

  result.is_integer = std::vector(d, false);

  result.implied_var_bounds = result.var_bounds;
  result.implied_is_integer = result.is_integer;

  result.validate();

  return result;
}

int main() {
  auto problem = setcover::read_problem(paths::resource("setcover/sc_63009_0"));

  auto greedy_solution = setcover::Greedy().solve(problem);
  auto greedy_result = setcover::evaluate(problem, greedy_solution);

  std::cout << "greedy: " << greedy_result.score << std::endl;

  auto milp = to_milp(problem);

  auto optimizer =
      presolve::Chain<double>()
          .add<presolve::TransformToEqualities<double>>()
          .add<presolve::RemoveLinearlyDependentEqualities<double>>()
          .add<presolve::Scaling<double>>();

  problem::StandardLP<double> standard_lp(optimizer.apply(milp));

  auto [n, d] = standard_lp.matrix.shape();
  std::println("{} x {}", n, d);

  // slightly perturb bounds
  // std::default_random_engine engine;
  // std::uniform_int_distribution distr(0, 5);
  // for (size_t i = 0; i < d; ++i) {
  //   if (standard_lp.var_bounds[i].lower) {
  //     *standard_lp.var_bounds[i].lower -=
  //         static_cast<double>(distr(engine)) * 1e-6;
  //   }
  //
  //   if (standard_lp.var_bounds[i].upper) {
  //     *standard_lp.var_bounds[i].upper +=
  //         static_cast<double>(distr(engine)) * 1e-6;
  //   }
  // }

  // run simplex algorithm
  simplex::Simplex<double> simplex;

  simplex.set_problem(standard_lp);
  simplex.set_validate_input(true);
  simplex.set_accountant<simplex::LoggingAccountant<double>>();

  // construct feasible point
  std::vector states(d, simplex::VariableState::AT_UPPER);
  for (size_t i = 0; i < n; ++i) {
    states[d - i - 1] = simplex::VariableState::BASIC;
  }

  auto result = simplex.primal(states);

  auto point = optimizer.inverse(simplex.get_point());
  auto cost = linalg::dot(point, milp.cost);

  std::cout << "simplex: " << cost << std::endl;
}
