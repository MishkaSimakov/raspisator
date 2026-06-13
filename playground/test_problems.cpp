#include <fstream>
#include <iostream>
#include <print>
#include <unordered_set>

#include "linear/simplex/Config.h"
#include "linear/simplex/Simplex.h"
#include "linear/simplex/init/primal/Phase1.h"
#include "mps/MPS.h"
#include "presolve/passes/RemoveLinearlyDependentEqualities.h"
#include "presolve/passes/TransformToEqualities.h"
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

    if (!problems.contains(problem_name)) {
      continue;
    }

    std::ifstream is(entry);
    if (!is) {
      throw std::runtime_error("Failed to open problem file.");
    }

    auto problem = mps::read<Field>(is, mps::Format::FIXED);
    std::println("{}: {} x {}", problem_name, problem.matrix.rows(),
                 problem.matrix.cols());

    problem = presolve::TransformToEqualities<Field>().apply(problem);
    problem =
        presolve::RemoveLinearlyDependentEqualities<Field>().apply(problem);
    problem = presolve::Scaling<Field>().apply(problem);

    problem::StandardMILP standard_problem(problem);

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

    std::visit(Overload{
                   [](const FiniteLPSolution<Field>& solution) {
                     std::println("  finite solution: {}", solution.value);
                   },
                   [](const NoFeasibleElements&) {
                     std::println("  no feasible elements");
                   },
                   [](const ReachedIterationsLimit<Field>&) {
                     std::println("  reached iterations limit");
                   },
                   [](const Unbounded&) { std::println("  unbounded"); },
               },
               solution.solution);
  }

  return 0;
}
