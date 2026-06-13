#include <gtest/gtest.h>

#include <variant>

#include "linear/simplex/Simplex.h"
#include "linear/simplex/init/dual/ReducedCost.h"

#include "faker/Catalog.h"
#include "linear/simplex/init/primal/Phase1.h"

TEST(CatalogTests, Dual) {
  const auto problems = faker::catalog<Rational>()
                            .problem_type(faker::ProblemType::StandardLP)
                            .solution_type(faker::SolutionType::FEASIBLE)
                            .all_variables_bounded(true)
                            .know_optimal_objective(true)
                            .all();

  for (auto instance : problems) {
    problem::StandardLP problem(instance.problem);

    auto states = simplex::try_init_dual_by_reduced_cost(problem);

    ASSERT_TRUE(states.has_value());

    simplex::Simplex<Rational> solver;

    solver.set_problem(problem);
    solver.set_validate_input(true);

    auto solution = solver.dual(*states).solution;

    ASSERT_TRUE(std::holds_alternative<FiniteLPSolution<Rational>>(solution));

    ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).value,
              *instance.optimal_objective);

    ASSERT_TRUE(simplex::is_primal_feasible(problem, solver.get_states()));
  }
}

TEST(CatalogTests, UnboundedPrimal) {
  const auto problems = faker::catalog<Rational>()
                            .problem_type(faker::ProblemType::StandardLP)
                            .solution_type(faker::SolutionType::UNBOUNDED)
                            .all();

  for (auto instance : problems) {
    problem::StandardLP problem(instance.problem);

    auto states = simplex::primal_phase1(problem);

    ASSERT_TRUE(states.has_value());

    simplex::Simplex<Rational> solver;

    solver.set_problem(problem);
    solver.set_validate_input(true);

    auto solution = solver.primal(*states).solution;

    ASSERT_TRUE(std::holds_alternative<Unbounded>(solution));
  }
}

TEST(CatalogTests, FeasiblePrimal) {
  const auto problems = faker::catalog<Rational>()
                            .problem_type(faker::ProblemType::StandardLP)
                            .solution_type(faker::SolutionType::FEASIBLE)
                            .know_optimal_objective(true)
                            .all();

  for (auto instance : problems) {
    problem::StandardLP problem(instance.problem);

    auto states = simplex::primal_phase1(problem);

    ASSERT_TRUE(states.has_value());

    simplex::Simplex<Rational> solver;

    solver.set_problem(problem);
    solver.set_validate_input(true);

    auto solution = solver.primal(*states).solution;

    ASSERT_TRUE(std::holds_alternative<FiniteLPSolution<Rational>>(solution));

    ASSERT_EQ(std::get<FiniteLPSolution<Rational>>(solution).value,
              *instance.optimal_objective);

    ASSERT_TRUE(simplex::is_primal_feasible(problem, solver.get_states()));
  }
}
