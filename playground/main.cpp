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
  size_t H = 5;
  Field first = 50;
  Field second = 0;

  auto problem = small_blomer_problem<Field>(first, second);

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
