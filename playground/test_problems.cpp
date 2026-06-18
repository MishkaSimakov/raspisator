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
#include "problem/mutations/RemoveRows.h"
#include "simplex/pricing/primal/SteepestEdge.h"

using Field = double;

int main() {
  const std::string problem_name = "BORE3D";

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
          .add<presolve::TransformToEqualities<Field>>();

  problem::StandardLP standard_problem(optimizer.apply(problem));

  auto phase1 = simplex::primal_phase1(
      standard_problem,
      simplex::Config<Field>()
          .set_validate_input(true)
          .set_accountant<simplex::LoggingAccountant<Field>>()
          .set_primal_pricing<simplex::PrimalMostInfeasible<Field>>()
          .set_max_iterations(100'000));

  if (!phase1) {
    std::println("Failed to find primal feasible basis.");
    return 0;
  }

  standard_problem =
      problem::remove_rows(std::move(standard_problem), phase1->redundant_rows);

  auto solver = simplex::Simplex<Field>();

  solver.set_validate_input(true);
  solver.set_problem(standard_problem);
  solver.set_accountant<simplex::LoggingAccountant<Field>>();
  solver.set_primal_pricing<simplex::PrimalMostInfeasible<Field>>();
  solver.set_max_iterations(100'000);

  auto result = solver.primal(phase1->states);

  auto end = std::chrono::steady_clock::now();
  std::println("status: {}, time: {}, iterations: {}", to_string(result.status),
               end - start, result.iterations_count + phase1->iterations_count);

  if (result.status == simplex::Status::OPTIMAL) {
    const auto objective =
        linalg::dot(optimizer.inverse(solver.get_point()), problem.cost);

    std::println("objective: {}", objective);
  }

  return 0;
}
