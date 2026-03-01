#pragma once

#include <unordered_set>

namespace setcover {

struct CoveringSet {
  size_t cost;
  std::unordered_set<size_t> elements;
};

struct Problem {
  std::vector<CoveringSet> sets;

  size_t elements_count;
};

struct Solution {
  std::vector<size_t> chosen_sets;
};

struct EvaluationResult {
  size_t score;
  bool is_valid;
};

}  // namespace setcover
