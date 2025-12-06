#include <chrono>
#include <iostream>
#include <print>

#include "encoding/UniformTimeDiscretization.h"
#include "linear/bb/PseudoCost.h"
#include "linear/bb/Settings.h"
#include "linear/bb/TreeStoringAccountant.h"
#include "linear/builder/ProblemBuilder.h"
#include "model/Solution.h"
#include "problems/Blomer.h"
#include "problems/Dwarf.h"
#include "utils/Drawing.h"

using Field = double;

int main() {
  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();

  size_t H = 10;

  auto problem = small_blomer_problem<Field>(100, 50);

  std::cout << to_graphviz(problem) << std::endl;

  auto encoding = to_uniform_time_milp(problem, H);

  std::cout << encoding.builder << std::endl;

  // solve MILP problem
  auto milp_problem = encoding.builder.get_problem();

  auto settings = BranchAndBoundSettings<Field>{.max_nodes = 1'000};
  auto solver = PseudoCostBranchAndBound<Field, TreeStoringAccountant<Field>>(
      milp_problem, settings);
  auto solution = solver.solve();

  std::cout << solver.get_accountant().to_graphviz() << std::endl;

  if (std::holds_alternative<NoFiniteSolution>(solution)) {
    std::println("No finite solution.");
  } else if (std::holds_alternative<ReachedNodesLimit>(solution)) {
    std::println("Reached nodes limit.");
  } else {
    auto finite_solution = std::get<FiniteMILPSolution<Field>>(solution);

    std::println("Finish production in {} time units.\n",
                 -finite_solution.value);

    auto point = finite_solution.point;

    for (const auto& unit : problem.get_units()) {
      std::println("schedule for unit {}:", unit.get_id());

      for (size_t t = 0; t < H; ++t) {
        for (const auto* task : unit.get_tasks() | std::views::keys) {
          Field x = encoding.builder.extract_variable(
              point, encoding.starts.at({&unit, task, t}));
          Field Q = encoding.builder.extract_variable(
              point, encoding.quantities.at({&unit, task, t}));

          std::println("x({}, {}, {}) = {}", unit.get_id(), task->get_id(), t,
                       x);
          std::println("Q({}, {}, {}) = {}", unit.get_id(), task->get_id(), t,
                       Q);
        }
      }
    }
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::println(
      "Overall: {}",
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin));
}

// Overall: 4478023834ns
