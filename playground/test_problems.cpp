#include <iostream>
#include <print>

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
  auto problems_path = paths::resource("lp_problems");
  for (auto entry : std::filesystem::directory_iterator{problems_path}) {
    if (entry.path().extension() != ".SIF") {
      continue;
    }

    // if (entry.path().filename().string() != "BORE3D.SIF") {
    // continue;
    // }

    auto reader = MPSReader<Field>(MPSFieldsMode::FIXED_WIDTH);
    reader.read(entry);

    auto problem = reader.get_canonical_representation();

    auto optimized = TransformToEqualities<Field>().apply(problem);
    optimized = RemoveLinearlyDependentConstraints<Field>().apply(optimized);
    auto matrices = to_matrices(optimized);

    std::println("{}: {} x {}", entry.path().filename().string(),
                 matrices.A.get_height(), matrices.A.get_width());

    auto solver =
        simplex::Simplex(CSCMatrix(matrices.A), matrices.b, matrices.c);

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
