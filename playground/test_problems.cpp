#include <fstream>
#include <iostream>
#include <print>
#include <unordered_set>

#include "mps/MPS.h"
#include "presolve/passes/RemoveLinearlyDependentEqualities.h"
#include "presolve/passes/TransformToEqualities.h"
#include "utils/Paths.h"

using Field = double;

int main() {
  std::unordered_set<std::string> problems = {
      // "SHELL"
      // "AFIRO", "ADLITTLE", "BANDM",
      // "BLEND", "PILOT"
  };

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

    std::println("{}", problem_name);
    auto problem = mps::read<Field>(is, mps::Format::FIXED);

    problem = presolve::TransformToEqualities<Field>().apply(problem);
    problem =
        presolve::RemoveLinearlyDependentEqualities<Field>().apply(problem);

    std::println("{}: {} x {}", problem_name, problem.matrix.rows(),
                 problem.matrix.cols());

    // simplex::Settings<Field> settings{.is_strict = true};
    // auto solver = simplex::Simplex<Field, simplex::LoggingAccountant<Field>>(
    //     CSCMatrix(matrices.A), matrices.b, matrices.c, settings);
    //
    // auto states = solver.get_primal_feasible(matrices.bounds);
    //
    // if (!states) {
    //   std::println("  Failed to find primal feasible basis.");
    //   continue;
    // }
    //
    // std::println("  Found primal feasible basis, starting solving.");
    // auto solution = solver.primal(matrices.bounds, *states);
    //
    // std::visit(Overload{
    //                [](const FiniteLPSolution<Field>& solution) {
    //                  std::println("  finite solution: {}", solution.value);
    //                },
    //                [](const NoFeasibleElements&) {
    //                  std::println("  no feasible elements");
    //                },
    //                [](const ReachedIterationsLimit<Field>&) {
    //                  std::println("  reached iterations limit");
    //                },
    //                [](const Unbounded&) { std::println("  unbounded"); },
    //            },
    //            solution.solution);
  }

  return 0;
}
