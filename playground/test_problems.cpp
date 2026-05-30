#include <fstream>
#include <iostream>
#include <print>
#include <unordered_set>

#include "linear/simplex/Config.h"
#include "linear/simplex/Simplex.h"
#include "mps/MPS.h"
#include "presolve/passes/RemoveLinearlyDependentEqualities.h"
#include "presolve/passes/TransformToEqualities.h"
#include "utils/Paths.h"

#include "linear/simplex/pricing/primal/MostInfeasible.h"

using Field = double;

int main() {
  std::unordered_set<std::string> problems = {
      "SHELL"
      "AFIRO",
      "ADLITTLE", "BANDM", "BLEND", "PILOT"};

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

    std::println("{}", problem_name);
    auto problem = mps::read<Field>(is, mps::Format::FIXED);

    problem = presolve::TransformToEqualities<Field>().apply(problem);
    problem =
        presolve::RemoveLinearlyDependentEqualities<Field>().apply(problem);

    std::println("{}: {} x {}", problem_name, problem.matrix.rows(),
                 problem.matrix.cols());

    auto A = problem.matrix;
    auto b = Matrix<Field>(problem.matrix.rows(), 1);
    auto c = Matrix<Field>(1, problem.matrix.cols());

    for (size_t i = 0; i < problem.matrix.rows(); ++i) {
      b[i, 0] = *problem.rhs_bounds[i].lower;
    }
    for (size_t i = 0; i < problem.matrix.cols(); ++i) {
      c[0, i] = problem.cost[i];
    }

    auto bounds = Bounds<Field>(problem.var_bounds);

    simplex::Config<Field> settings{.is_strict = true};
    auto solver = simplex::Simplex<Field, simplex::LoggingAccountant<Field>>(
        A, b, c,
        {
            .is_strict = true,
            .primal_pricing =
                std::make_unique<simplex::PrimalMostInfeasible<Field>>(),
        });

    auto states = solver.get_primal_feasible(bounds);

    if (!states) {
      std::println("  Failed to find primal feasible basis.");
      continue;
    }

    std::println("  Found primal feasible basis, starting solving.");
    auto solution = solver.primal(bounds, *states);

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
