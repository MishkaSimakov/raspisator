#pragma once

#include <memory>
#include <ranges>
#include <vector>

#include "Pass.h"

namespace presolve {

template <typename Field>
class Chain final : Pass<Field> {
  std::vector<std::unique_ptr<Pass<Field>>> passes_;

 public:
  Chain() = default;

  Chain(const Chain&) = delete;
  Chain& operator=(const Chain&) = delete;

  Chain(Chain&&) = default;
  Chain& operator=(Chain&&) = default;

  template <typename T, typename... Args>
  Chain& add(Args&&... args) & {
    passes_.push_back(std::make_unique<T>(std::forward<Args>(args)...));
    return *this;
  }

  template <typename T, typename... Args>
  Chain&& add(Args&&... args) && {
    passes_.push_back(std::make_unique<T>(std::forward<Args>(args)...));
    return std::move(*this);
  }

  problem::MILP<Field> apply(problem::MILP<Field> problem) override {
    for (auto& pass : passes_) {
      problem = pass->apply(std::move(problem));
    }

    return problem;
  }

  Vector<Field> inverse(Vector<Field> solution) const override {
    for (auto& pass : passes_ | std::views::reverse) {
      solution = pass->inverse(std::move(solution));
    }

    return solution;
  }
};

}  // namespace presolve
