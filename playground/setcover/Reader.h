#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>

#include "Types.h"

namespace setcover {

Problem read_problem(const std::filesystem::path& path) {
  std::ifstream is(path);

  if (!is) {
    throw std::runtime_error("Failed to open file.");
  }

  size_t elements_count;
  size_t sets_count;
  is >> elements_count >> sets_count;

  std::vector<CoveringSet> sets(sets_count);
  std::string buffer;

  // skip current line
  std::getline(is, buffer);

  for (size_t i = 0; i < sets_count; ++i) {
    std::getline(is, buffer);
    std::stringstream line(buffer);

    line >> sets[i].cost;

    assert(sets[i].cost > 0 && "cost is expected to be strictly positive");

    while (!line.eof()) {
      size_t element;
      line >> element;

      sets[i].elements.insert(element);
    }
  }

  return Problem{
      .sets = std::move(sets),
      .elements_count = elements_count,
  };
}

}  // namespace setcover
