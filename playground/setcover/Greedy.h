#pragma once

#include <vector>

#include "Types.h"

namespace setcover {

class Greedy {
  static size_t argmax_set(const std::vector<CoveringSet>& sets) {
    size_t max_index = 0;
    double max_value = 0;

    for (size_t i = 0; i < sets.size(); ++i) {
      double value = static_cast<double>(sets[i].elements.size()) /
                     static_cast<double>(sets[i].cost);

      if (max_value < value) {
        max_value = value;
        max_index = i;
      }
    }

    return max_index;
  }

 public:
  Greedy() = default;

  Solution solve(const Problem& problem) {
    size_t sets_count = problem.sets.size();

    std::vector<size_t> result;

    // копируем все множества, так как далее из них будут убираться уже покрытые
    // элементы
    auto sets = problem.sets;

    while (true) {
      size_t chosen_set = argmax_set(sets);

      if (sets[chosen_set].elements.empty()) {
        break;
      }

      result.push_back(chosen_set);

      for (size_t i = 0; i < sets_count; ++i) {
        if (i == chosen_set) {
          continue;
        }

        for (size_t element : sets[chosen_set].elements) {
          sets[i].elements.erase(element);
        }
      }
      sets[chosen_set].elements.clear();
    }

    return Solution{std::move(result)};
  }
};

}  // namespace setcover
