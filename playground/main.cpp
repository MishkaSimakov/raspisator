#include <chrono>
#include <fstream>
#include <iostream>
#include <print>

#include "encoding/UniformTimeDiscretization.h"
#include "linear/BigInteger.h"
#include "linear/Scaler.h"
#include "linear/bb/PseudoCost.h"
#include "linear/bb/Settings.h"
#include "linear/bb/TreeStoringAccountant.h"
#include "linear/builder/ProblemBuilder.h"
#include "model/Solution.h"
#include "problems/Blomer.h"
#include "utils/Drawing.h"

using Field = double;

int main() {
  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();

  size_t H = 12;

  auto problem = small_blomer_problem<Field>(100, 100);

  std::cout << to_graphviz(problem) << std::endl;

  auto encoding = to_uniform_time_milp(problem, H);

  std::cout << encoding.builder << std::endl;

  // solve MILP problem
  auto milp_problem = encoding.builder.get_problem();

  auto scaled_problem = Scaler(milp_problem).get_scaled();
  std::cout << "scaling quality: "
            << Scaler(scaled_problem).get_scaling_quality() << std::endl;

  auto settings = BranchAndBoundSettings<Field>{
      .max_nodes = 10'000,
      .perturbation = PerturbationMode::DISABLED,
  };
  auto solver = PseudoCostBranchAndBound<Field, TreeStoringAccountant<Field>>(
      scaled_problem, settings);
  auto solution = solver.solve();

  // std::cout << solver.get_accountant().to_graphviz() << std::endl;

  // {
  //   std::ofstream os("iterations_data.csv");
  //   os << solver.get_accountant().iterations_to_csv();
  // }
  //
  // {
  //   std::ofstream os("strong_branching_data.csv");
  //   os << solver.get_accountant().strong_branching_iterations_to_csv();
  // }

  {
    std::ofstream os("tree.dot");
    os << solver.get_accountant().to_graphviz() << std::endl;
  }

  // print solution
  if (std::holds_alternative<NoFiniteSolution>(solution)) {
    std::println("No finite solution.");
  } else if (std::holds_alternative<ReachedNodesLimit>(solution)) {
    std::println("Reached nodes limit.");
  } else {
    auto finite_solution = std::get<FiniteMILPSolution<Field>>(solution);

    std::println("Finish production in {} time units.\n",
                 -finite_solution.value);

    auto point = finite_solution.point;

    Solution checker(&problem);

    for (const auto& unit : problem.get_units()) {
      std::println("schedule for unit {}:", unit.get_id());
      std::println("(unit, task, time)");

      for (size_t t = 0; t < H; ++t) {
        for (const auto* task : unit.get_tasks() | std::views::keys) {
          Field x = encoding.builder.extract_variable(
              point, encoding.starts.at({&unit, task, t}));
          Field Q = encoding.builder.extract_variable(
              point, encoding.quantities.at({&unit, task, t}));

          if (FieldTraits<Field>::is_strictly_positive(x)) {
            checker.add_instance(
                TaskInstance{unit.get_id(), task->get_id(), Q, t});
          }

          std::println("x({}, {}, {}) = {}", unit.get_id(), task->get_id(), t,
                       x);
          std::println("Q({}, {}, {}) = {}", unit.get_id(), task->get_id(), t,
                       Q);
        }
      }
    }

    if (!checker.check()) {
      throw std::runtime_error("Invalid solution!");
    }
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::println(
      "Overall: {}",
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin));
}

// Overall: 4478023834ns
