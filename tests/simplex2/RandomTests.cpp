#include <gtest/gtest.h>

#include <variant>

#include "Assertions.h"
#include "linear/BigInteger.h"
#include "linear/model/LP.h"
#include "linear/simplex/Feasibility.h"
#include "linear/simplex/Simplex.h"
#include "linear/simplex/init/dual/ReducedCost.h"
#include "support/StandardMILPBuilder.h"

constexpr size_t kIterations = 500;
constexpr size_t kRows = 5;
constexpr size_t kCols = 10;
constexpr int kMagnitude = 10;

TEST(Simplex2RandomTests, DualFeasible) {
  for (size_t i = 0; i < kIterations; ++i) {
    auto [prob, witness, primal_states] = StandardMILPBuilder<Rational>{}
                                              .rows(kRows)
                                              .cols(kCols)
                                              .magnitude(kMagnitude)
                                              .seed(i)
                                              .build_feasible();

    auto states = simplex::try_init_dual_by_reduced_cost(
        prob.matrix, prob.rhs, prob.cost, prob.var_bounds);

    // try_init_dual_by_reduced_cost is guaranteed to succeed when all
    // variables have both finite bounds
    ASSERT_TRUE(states.has_value()) << "Dual init failed at iteration " << i;

    simplex::Simplex<Rational> solver;
    solver.set_problem(prob);

    auto result = solver.dual(*states);

    ASSERT_TRUE(result.is_feasible())
        << "Dual simplex failed to find solution at iteration " << i;

    auto& sol = std::get<FiniteLPSolution<Rational>>(result.solution);
    ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(prob, sol))
        << "Solution validation failed at iteration " << i;
  }
}

TEST(Simplex2RandomTests, PrimalFeasible) {
  for (size_t i = 0; i < kIterations; ++i) {
    auto [prob, witness, primal_states] = StandardMILPBuilder<Rational>{}
                                              .rows(kRows)
                                              .cols(kCols)
                                              .magnitude(kMagnitude)
                                              .seed(i)
                                              .build_feasible();

    // primal_states are computed by construction — no Phase1 needed
    ASSERT_TRUE(simplex::is_primal_feasible(prob, primal_states))
        << "Builder produced invalid primal states at iteration " << i;

    simplex::Simplex<Rational> solver;
    solver.set_problem(prob);

    auto result = solver.primal(primal_states);

    ASSERT_TRUE(result.is_feasible())
        << "Primal simplex failed to find solution at iteration " << i;

    auto& sol = std::get<FiniteLPSolution<Rational>>(result.solution);
    ASSERT_NO_FATAL_FAILURE(validate_simplex_solution(prob, sol))
        << "Solution validation failed at iteration " << i;
  }
}

TEST(Simplex2RandomTests, Infeasible) {
  for (size_t i = 0; i < kIterations; ++i) {
    auto prob = StandardMILPBuilder<Rational>{}
                    .rows(kRows)
                    .cols(kCols)
                    .magnitude(kMagnitude)
                    .seed(i)
                    .build_infeasible();

    // All variables have both bounds, so dual init is guaranteed to succeed
    auto states = simplex::try_init_dual_by_reduced_cost(
        prob.matrix, prob.rhs, prob.cost, prob.var_bounds);

    ASSERT_TRUE(states.has_value())
        << "Dual init failed for infeasible problem at iteration " << i;

    simplex::Simplex<Rational> solver;
    solver.set_problem(prob);

    auto result = solver.dual(*states);

    EXPECT_TRUE(std::holds_alternative<NoFeasibleElements>(result.solution))
        << "Expected infeasible result at iteration " << i
        << " but got feasible solution";
  }
}

TEST(Simplex2RandomTests, Unbounded) {
  for (size_t i = 0; i < kIterations; ++i) {
    auto [prob, primal_states] = StandardMILPBuilder<Rational>{}
                                     .rows(kRows)
                                     .magnitude(kMagnitude)
                                     .seed(i)
                                     .build_unbounded();

    ASSERT_TRUE(simplex::is_primal_feasible(prob, primal_states))
        << "Builder produced invalid primal states for unbounded problem at "
           "iteration "
        << i;

    simplex::Simplex<Rational> solver;
    solver.set_problem(prob);

    auto result = solver.primal(primal_states);

    EXPECT_TRUE(std::holds_alternative<Unbounded>(result.solution))
        << "Expected unbounded result at iteration " << i;
  }
}
