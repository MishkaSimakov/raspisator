#pragma once
#include "Types.h"

namespace simplex {

template <typename Field>
class EmptyAccountant {
 public:
  void iteration(IterationState<Field> state) {}
};

template <typename Field>
class LoggingAccountant {
  using Clock = std::chrono::high_resolution_clock;

  Clock::time_point last_time_;
  size_t iterations_since_last_time_;

 public:
  LoggingAccountant()
      : last_time_(Clock::now()), iterations_since_last_time_(0) {}

  void iteration(IterationState<Field> state) {
    ++iterations_since_last_time_;
    auto curr_time = Clock::now();

    if (curr_time - last_time_ > std::chrono::seconds{1}) {
      double speed =
          static_cast<double>(iterations_since_last_time_) /
          std::chrono::duration<double>(curr_time - last_time_).count();

      std::println("{:.1f} itr/s, objective: {}, size: {}, elapsed: {}", speed,
                   state.objective, state.lupa.size(), curr_time - last_time_);

      last_time_ = curr_time;
      iterations_since_last_time_ = 0;
    }
  }
};

}  // namespace simplex
