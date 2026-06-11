#include <gtest/gtest.h>

#include <variant>

#include "linear/simplex/Simplex.h"
#include "linear/simplex/init/dual/ReducedCost.h"

#include "faker/Catalog.h"
#include "linear/simplex/init/primal/Phase1.h"

TEST(CatalogTests, Dual) {
  const auto selector = faker::Tag::KNOWN_OPTIMAL_OBJECTIVE |
                        faker::Tag::ALL_VARIABLES_BOUNDED |
                        faker::Tag::STANDARD_LP;

  const auto problems = faker::catalog<Rational>(selector);

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
  const auto selector = faker::Tag::UNBOUNDED | faker::Tag::STANDARD_LP;

  const auto problems = faker::catalog<Rational>(selector);

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
  const auto selector = faker::Tag::FEASIBLE |
                        faker::Tag::KNOWN_OPTIMAL_OBJECTIVE |
                        faker::Tag::STANDARD_LP;

  const auto problems = faker::catalog<Rational>(selector);

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
