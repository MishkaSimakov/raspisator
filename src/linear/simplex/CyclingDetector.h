#pragma once

#include <print>
#include <unordered_map>
#include <vector>

#include "linear/FieldTraits.h"
#include "utils/Hashers.h"
#include "linear/model/LP.h"

namespace simplex {

enum class CyclingState { NORMAL, HAS_CYCLING };

// Simplex method cycling detection system.
// It stores iteration on which a basis was last visited and objective value
// on this iteration.
// It is assumed that during simplex method objective value is not increasing.
template <typename Field>
class CyclingDetector {
  std::unordered_map<size_t, std::pair<size_t, Field>> visited_bases;

  static size_t hash_states(const std::vector<VariableState>& states) {
    StreamHasher hasher;

    for (auto state : states) {
      hasher << static_cast<size_t>(state);
    }

    return hasher.get_hash();
  }

 public:
  CyclingState record(size_t iteration,
                      const std::vector<VariableState>& states,
                      Field objective) {
    const auto hash = hash_states(states);
    const auto [itr, was_emplaced] =
        visited_bases.emplace(hash, std::pair{iteration, objective});

    if (!was_emplaced) {
      if (FieldTraits<Field>::is_nonzero(objective - itr->second.second)) {
        // cache collision
        visited_bases.erase(itr);
      } else {
        std::println(
            "Cycling! iteration delta: {}, iteration: {}, last visited on "
            "iteration: {}",
            iteration - itr->second.first, iteration, itr->second.first);

        itr->second = std::pair{iteration, objective};

        return CyclingState::HAS_CYCLING;
      }
    }

    return CyclingState::NORMAL;
  }

  void clear() { visited_bases.clear(); }
};

}  // namespace simplex
