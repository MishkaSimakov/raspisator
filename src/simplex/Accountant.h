#pragma once

#include "Move.h"
#include "StateView.h"

namespace simplex {

template <typename Field>
class Accountant {
 public:
  virtual void iteration(StateView<Field> simplex) {}

  virtual void suspicious_pivot(StateView<Field> simplex,
                                ChangeBasisMove<Field> move) {}

  // @culprit is the index of violating variable in simplex.basic_vars array
  virtual void violate_primal_bounds(StateView<Field> simplex, size_t culprit) {
  }

  virtual ~Accountant() = default;
};

//

template <typename Field>
class LoggingAccountant final : public Accountant<Field> {
  using Clock = std::chrono::high_resolution_clock;

  Clock::time_point last_time_;
  size_t iterations_since_last_time_;

 public:
  LoggingAccountant()
      : last_time_(Clock::now()), iterations_since_last_time_(0) {}

  void iteration(StateView<Field> simplex) override {
    ++iterations_since_last_time_;
    auto curr_time = Clock::now();

    if (curr_time - last_time_ > std::chrono::seconds{1}) {
      double speed =
          static_cast<double>(iterations_since_last_time_) /
          std::chrono::duration<double>(curr_time - last_time_).count();

      std::println("  [{}] {:.1f} itr/s, objective: {}", simplex.iteration,
                   speed, simplex.objective);

      last_time_ = curr_time;
      iterations_since_last_time_ = 0;
    }
  }

  void suspicious_pivot(StateView<Field> simplex,
                        ChangeBasisMove<Field> move) override {
    std::println("  [{}] suspicious pivot", simplex.iteration);
  }

  void violate_primal_bounds(StateView<Field> simplex,
                             size_t culprit) override {
    std::println("  [{}] primal violation: variable {} with value {} not in {}",
                 simplex.problem.var_name(simplex.basic_vars[culprit]),
                 simplex.basic_point[culprit],
                 simplex.problem.var_bounds[simplex.basic_vars[culprit]]);
  }
};

}  // namespace simplex
