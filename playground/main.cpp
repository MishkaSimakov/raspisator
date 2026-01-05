#include <chrono>
#include <print>

#include "linear/matrix/NPY.h"
#include "linear/problem/MPS.h"
#include "linear/scheduling/BlomersHeuristic.h"
#include "problems/Blomer.h"
#include "problems/Dwarf.h"
#include "utils/Drawing.h"
#include "utils/ShadowFloat.h"

using Field = double;

int main() {
  const std::filesystem::path problems{"resources/lp_problems"};

  const std::string small_files[] = {
    "BOEING2"
  };

  for (const auto& dir_entry : small_files) {
    // if (dir_entry.path().extension() != ".SIF") {
    //   continue;
    // }

    std::string name = problems.string() + "/" + dir_entry + ".SIF";

    std::println("Solving {}", name);

    // read problem
    MPSReader<Field> reader;
    reader.read(name);
    auto problem = reader.get_canonical_representation();

    std::cout << problem << std::endl;

    problem = RemoveLinearlyDependentConstraints<Field>().apply(problem);
    problem = TransformToEqualities<Field>().apply(problem);

    auto matrices = to_matrices(problem);

    // run simplex method
    auto solver = simplex::BoundedSimplexMethod(CSCMatrix(matrices.A),
                                                matrices.b, matrices.c);

    auto run_result = solver.dual(matrices.lower, matrices.upper);

    std::println("  iterations: {}", run_result.iterations_count);

    if (run_result.is_feasible()) {
      std::println(
          "  solution: {}",
          std::get<FiniteLPSolution<Field>>(run_result.solution).value);
    } else {
      std::println("  infeasible");
    }
  }
}

// size_t H = 12;
// Field first = 100;
// Field second = 100;
//
// auto problem = small_blomer_problem<Field>(first, second);
//
// auto begin = std::chrono::steady_clock::now();
//
// std::println("precise small_blomer_{}_{} ({})", first, second, H);
//
// TimeGrid time_grid;
//
// time_grid.emplace(0, [](size_t time) { return true; });
// time_grid.emplace(1, [](size_t time) { return true; });
// time_grid.emplace(2, [](size_t time) { return true; });
// time_grid.emplace(3, [](size_t time) { return true; });
//
// auto solution = apply_time_grid_model(problem, H, time_grid);
//
// if (!solution) {
//   std::println("No solution");
// }
//
// if (!solution->check()) {
//   throw std::runtime_error("Invalid solution!");
// }
//
// std::println("    periods: {}", solution->get_total_time());
//
// auto end = std::chrono::steady_clock::now();
// std::println(
//     " Overall: {}",
//     std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin));
