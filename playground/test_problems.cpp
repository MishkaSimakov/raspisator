#include <iostream>
#include <print>
#include <unordered_set>

#include "linear/problem/MPS.h"
#include "linear/problem/optimization/FullOptimizer.h"
#include "linear/problem/optimization/RemoveLinearlyDependentConstraints.h"
#include "linear/problem/optimization/TransformToEqualities.h"
#include "linear/simplex/Simplex.h"
#include "utils/Paths.h"
#include "utils/ShadowFloat.h"
#include "utils/Variant.h"

using Field = double;

int main() {
  std::unordered_set<std::string> problems = {
      // "AFIRO",
      // "ADLITTLE",
      // "BANDM", "BLEND",
      "PILOT"
  };

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

    auto reader = MPSReader<Field>(MPSFieldsMode::FIXED_WIDTH);
    reader.read(entry);

    auto problem = reader.get_canonical_representation();

    problem = Scaling<Field>().apply(problem);
    problem = TransformToEqualities<Field>().apply(problem);
    problem = RemoveLinearlyDependentConstraints<Field>().apply(problem);

    auto matrices = to_matrices(problem);

    std::println("{}: {} x {}", problem_name, matrices.A.get_height(),
                 matrices.A.get_width());

    simplex::Settings<Field> settings{.is_strict = true};
    auto solver = simplex::Simplex<Field, simplex::LoggingAccountant<Field>>(
        CSCMatrix(matrices.A), matrices.b, matrices.c, settings);

    auto states = solver.try_get_primal_feasible(matrices.bounds);

    if (!states) {
      std::println("  Failed to find primal feasible basis.");
      continue;
    }

    auto solution = solver.primal(matrices.bounds, *states);

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
