#include <chrono>
#include <print>

#include "linear/scheduling/BlomersHeuristic.h"
#include "problems/Blomer.h"
#include "problems/Dwarf.h"
#include "utils/Drawing.h"

using Field = double;

int main() {
  std::vector<std::pair<Field, Field>> targets = {
      {200, 200}, {50, 50}, {0, 100}, {100, 0}, {100, 100},
  };

  for (auto [first, second] : targets) {
    size_t H = 20;

    auto problem = small_blomer_problem<Field>(first, second);

    {
      auto begin = std::chrono::steady_clock::now();

      std::println("heuristic small_blomer_{}_{} ({})", first, second, H);

      TimeGrid time_grid;

      time_grid.emplace(0, [](size_t time) { return true; });
      time_grid.emplace(1, [](size_t time) { return time % 2 == 0; });
      time_grid.emplace(2, [](size_t time) { return time % 2 == 0; });
      time_grid.emplace(3, [](size_t time) { return true; });

      auto solution = blomer_heuristic_model(problem, H, time_grid);

      if (!solution) {
        std::println("No solution");
      }

      if (!solution->check()) {
        throw std::runtime_error("Invalid solution!");
      }
      //
      // solution->to_graphviz(std::cout);

      auto end = std::chrono::steady_clock::now();
      std::println(
          " Overall: {}",
          std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin));
    }

    {
      auto begin = std::chrono::steady_clock::now();

      std::println("precise small_blomer_{}_{} ({})", first, second, H);

      TimeGrid time_grid;

      time_grid.emplace(0, [](size_t time) { return true; });
      time_grid.emplace(1, [](size_t time) { return true; });
      time_grid.emplace(2, [](size_t time) { return true; });
      time_grid.emplace(3, [](size_t time) { return true; });

      auto solution = apply_time_grid_model(problem, H, time_grid);

      if (!solution) {
        std::println("No solution");
      }

      if (!solution->check()) {
        throw std::runtime_error("Invalid solution!");
      }

      std::println("    periods: {}", solution->get_total_time());

      auto end = std::chrono::steady_clock::now();
      std::println(
          " Overall: {}",
          std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin));
    }
  }
}
