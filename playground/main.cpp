#include <chrono>
#include <print>

#include "linear/scheduling/BlomersHeuristic.h"
#include "problems/Blomer.h"
#include "problems/Dwarf.h"
#include "utils/Drawing.h"

using Field = double;

int main() {
  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();

  size_t H = 20;

  auto problem = small_blomer_problem<Field>(100, 100);


  TimeGrid time_grid;

  time_grid.emplace(0, [](size_t time) { return true; });
  time_grid.emplace(1, [](size_t time) { return time % 2 == 0; });
  time_grid.emplace(2, [](size_t time) { return time % 2 == 0; });
  time_grid.emplace(3, [](size_t time) { return true; });

  // time_grid.emplace(0, [](size_t time) { return time % 2 == 0; });
  // time_grid.emplace(1, [](size_t time) { return time % 4 == 0; });
  // time_grid.emplace(2, [](size_t time) { return time % 2 == 0; });
  // time_grid.emplace(3, [](size_t time) { return time % 4 == 0; });
  // time_grid.emplace(4, [](size_t time) { return time % 6 == 0; });
  // time_grid.emplace(5, [](size_t time) { return time % 6 == 0; });
  // time_grid.emplace(6, [](size_t time) { return time % 6 == 0; });
  // time_grid.emplace(7, [](size_t time) { return time % 6 == 0; });
  // time_grid.emplace(8, [](size_t time) { return (time + 1) % 6 == 0; });

  auto solution = blomer_heuristic_model(problem, H, time_grid);

  if (!solution) {
    std::println("No solution");
    return 0;
  }

  if (!solution->check()) {
    throw std::runtime_error("Invalid solution!");
  }

  solution->to_graphviz(std::cout);

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::println(
      "Overall: {}",
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin));
}

// #include <chrono>
// #include <print>
//
// #include "linear/scheduling/BlomersHeuristic.h"
// #include "problems/Blomer.h"
// #include "problems/Dwarf.h"
//
// using Field = double;
//
// int main() {
//   auto problem = dwarf_problem_smaller<Field>();
//
//   Solution<Field> solution(&problem);
//   solution.add_instance(TaskInstance<Field>{0, 0, 50, 0});
//   solution.add_instance(TaskInstance<Field>{0, 0, 50, 2});
//
//   if (!solution.check()) {
//     throw std::runtime_error("wrong solution!");
//   }
//
//   auto refined = apply_left_shift_model(problem, 4, solution);
//
//   std::cout << refined.has_value() << std::endl;
//
//   refined->to_graphviz(std::cout);
// }
