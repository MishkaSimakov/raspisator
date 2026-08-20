#pragma once

#include <algorithm>
#include <cassert>
#include <random>

namespace rnd {

template <typename Gen>
  requires std::uniform_random_bit_generator<Gen>
size_t index(size_t size, Gen& generator) {
  assert(size != 0);

  return std::uniform_int_distribution<size_t>(0, size - 1)(generator);
}

// Uniformly generates a set of `count` random indices from 0 to (`size` - 1)
// Result is a sorted vector of indices in descending order.
template <typename Gen>
  requires std::uniform_random_bit_generator<Gen>
std::vector<size_t> unique_indices(size_t size, size_t count, Gen& generator) {
  assert(size != 0 && count <= size);

  std::vector<size_t> result(count);
  for (size_t i = 0; i < count; ++i) {
    result[i] = rnd::index(size - i, generator);

    for (size_t j = 0; j < i; ++j) {
      if (result[i - j] >= result[i - j - 1]) {
        ++result[i - j];
        std::swap(result[i - j], result[i - j - 1]);
      } else {
        break;
      }
    }
  }

  return result;
}

// Returns true with given probability
template <typename Gen>
  requires std::uniform_random_bit_generator<Gen>
bool bernoulli(double probability, Gen& generator) {
  assert(0 <= probability && probability <= 1);

  return std::bernoulli_distribution(probability)(generator);
}

}  // namespace rnd
