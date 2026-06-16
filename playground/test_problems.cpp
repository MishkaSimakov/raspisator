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

// WOOD1P fails to find primal feasible in phase1
int main() {
  std::unordered_set<std::string> problems = {// "SHELL"
                                              // "AFIRO",
                                              // "ADLITTLE",
                                              // "BANDM",
                                              // "BLEND",
                                              // "PILOT"
                                              // "PEROLD",
                                              //"BNL2"
                                              "D6CUBE"};

  auto problems_path = paths::resource("lp_problems");

  std::vector<std::filesystem::directory_entry> all_problems;
  std::ranges::copy(std::filesystem::directory_iterator{problems_path},
                    std::back_inserter(all_problems));

  for (size_t index = 0; index < all_problems.size(); ++index) {
    auto path = all_problems[index].path();

    if (path.extension() != ".SIF") {
      continue;
    }

    path.replace_extension("");

    const auto problem_name = path.filename().string();

    // if (!problems.contains(problem_name)) {
    // continue;
    // }

    std::println("{:02}/{} {}", index + 1, all_problems.size(), problem_name);

    std::ifstream is(all_problems[index]);
    if (!is) {
      std::println("  Failed to open problem file.");
    }

    auto problem = mps::read<Field>(is, mps::Format::FIXED);
    std::println("  size: {} x {}", problem.matrix.rows(),
                 problem.matrix.cols());

    auto start = std::chrono::steady_clock::now();

    auto optimizer =
        presolve::Chain<Field>()
            .add<presolve::TransformToEqualities<Field>>()
            .add<presolve::RemoveLinearlyDependentEqualities<Field>>()
            .add<presolve::Scaling<Field>>();

    problem::StandardMILP standard_problem(optimizer.apply(problem));

    size_t phase1_iterations;
    auto states = simplex::primal_phase1(
        standard_problem,
        simplex::Config<Field>()
            .set_validate_input(true)
            // .set_accountant<simplex::LoggingAccountant<Field>>()
            .set_primal_pricing<simplex::PrimalMostInfeasible<Field>>()
            .set_max_iterations(100'000),
        &phase1_iterations);

    if (!states) {
      std::println("  Failed to find primal feasible basis.");
      continue;
    }

    auto solver = simplex::Simplex<Field>();

    solver.set_validate_input(true);
    solver.set_problem(standard_problem);
    // solver.set_accountant<simplex::LoggingAccountant<Field>>();
    solver.set_primal_pricing<simplex::PrimalMostInfeasible<Field>>();
    solver.set_max_iterations(100'000);

    auto result = solver.primal(*states);

    auto end = std::chrono::steady_clock::now();
    std::println("  status: {}, time: {}, iterations: {}",
                 to_string(result.status), end - start,
                 result.iterations_count + phase1_iterations);

    if (result.status == simplex::Status::OPTIMAL) {
      const auto objective =
          linalg::dot(optimizer.inverse(solver.get_point()), problem.cost);

      std::println("  objective: {}", objective);
    }
  }

  return 0;
}
