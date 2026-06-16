#include <fstream>
#include <iostream>
#include <print>
#include <unordered_set>

#include "mps/MPS.h"
#include "presolve/Chain.h"
#include "presolve/passes/RemoveLinearlyDependentEqualities.h"
#include "presolve/passes/TransformToEqualities.h"
#include "simplex/Config.h"
#include "simplex/Simplex.h"
#include "simplex/init/primal/Phase1.h"
#include "utils/Paths.h"

#include "presolve/passes/Scaling.h"
#include "problem/StandardMILP.h"
#include "simplex/pricing/primal/SteepestEdge.h"

using Field = double;

// - WOOD1P fails to find primal feasible in phase1 because row reduction is not
// numerically stable enough
// - QAP15 RemoveLinearlyDependentEqualities takes a lot of time, and phase1
// gets stuck in degenerate iterations
int main() {
  const std::string problem_name = "QAP15";

  auto path = paths::resource(std::format("lp_problems/{}.SIF", problem_name));

  std::ifstream is(path);
  if (!is) {
    std::println("Failed to open problem file.");
    return 0;
  }

  auto problem = mps::read<Field>(is, mps::Format::FIXED);
  std::println("size: {} x {}", problem.matrix.rows(), problem.matrix.cols());

  auto start = std::chrono::steady_clock::now();

  auto optimizer =
      presolve::Chain<Field>()
          .add<presolve::Scaling<Field>>()
          .add<presolve::RemoveLinearlyDependentEqualities<Field>>()
          .add<presolve::TransformToEqualities<Field>>();

  problem::StandardMILP standard_problem(optimizer.apply(problem));

  size_t phase1_iterations;
  auto states = simplex::primal_phase1(
      standard_problem,
      simplex::Config<Field>()
          .set_validate_input(true)
          .set_accountant<simplex::LoggingAccountant<Field>>()
          .set_primal_pricing<simplex::PrimalMostInfeasible<Field>>()
          .set_max_iterations(100'000),
      &phase1_iterations);

  if (!states) {
    std::println("Failed to find primal feasible basis.");
    return 0;
  }

  auto solver = simplex::Simplex<Field>();

  solver.set_validate_input(true);
  solver.set_problem(standard_problem);
  solver.set_accountant<simplex::LoggingAccountant<Field>>();
  solver.set_primal_pricing<simplex::PrimalMostInfeasible<Field>>();
  solver.set_max_iterations(100'000);

  auto result = solver.primal(*states);

  auto end = std::chrono::steady_clock::now();
  std::println("status: {}, time: {}, iterations: {}", to_string(result.status),
               end - start, result.iterations_count + phase1_iterations);

  if (result.status == simplex::Status::OPTIMAL) {
    const auto objective =
        linalg::dot(optimizer.inverse(solver.get_point()), problem.cost);

    std::println("objective: {}", objective);
  }

  return 0;
}
