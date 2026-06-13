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

using Field = double;

int main() {
  std::unordered_set<std::string> problems = {// "SHELL"
                                              // "AFIRO",
                                              // "ADLITTLE", "BANDM", "BLEND",
                                              // "PILOT"
                                              "PEROLD"};

  auto problems_path = paths::resource("lp_problems");
  for (auto entry : std::filesystem::directory_iterator{problems_path}) {
    if (entry.path().extension() != ".SIF") {
      continue;
    }

    auto path = entry.path();
    path.replace_extension("");

    auto problem_name = path.filename().string();

    // if (!problems.contains(problem_name)) {
    // continue;
    // }

    std::ifstream is(entry);
    if (!is) {
      throw std::runtime_error("Failed to open problem file.");
    }

    auto problem = mps::read<Field>(is, mps::Format::FIXED);
    std::println("{}: {} x {}", problem_name, problem.matrix.rows(),
                 problem.matrix.cols());

    auto optimizer =
        presolve::Chain<Field>()
            .add<presolve::TransformToEqualities<Field>>()
            .add<presolve::RemoveLinearlyDependentEqualities<Field>>()
            .add<presolve::Scaling<Field>>();

    problem::StandardMILP standard_problem(optimizer.apply(problem));

    auto states = simplex::primal_phase1(
        standard_problem,
        simplex::Config<Field>()
            .set_validate_input(true)
            .set_accountant<simplex::LoggingAccountant<Field>>());

    if (!states) {
      std::println("  Failed to find primal feasible basis.");
      continue;
    }

    std::println("  Found primal feasible basis, starting solving.");

    auto solver = simplex::Simplex<Field>(
        simplex::Config<Field>()
            .set_accountant<simplex::LoggingAccountant<Field>>());

    solver.set_validate_input(true);
    solver.set_problem(standard_problem);

    auto solution = solver.primal(*states);

    std::visit(
        Overload{
            [&optimizer,
             &problem](const simplex::FiniteLPSolution<Field>& solution) {
              auto objective =
                  linalg::dot(optimizer.inverse(solution.point), problem.cost);

              std::println("  finite solution: {}", objective);
            },
            [](const simplex::NoFeasibleElements&) {
              std::println("  no feasible elements");
            },
            [](const simplex::ReachedIterationsLimit<Field>&) {
              std::println("  reached iterations limit");
            },
            [](const simplex::Unbounded&) { std::println("  unbounded"); },
        },
        solution.solution);
  }

  return 0;
}
