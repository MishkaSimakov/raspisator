#pragma once

#include <vector>

#include "State.h"

class Task {
  size_t id_;

  // TODO: different proportions for outputs
  std::vector<const State*> inputs_;
  std::vector<const State*> outputs_;

  void set_id(size_t id) { id_ = id; }

 public:
  size_t get_id() const { return id_; }

  template <typename... Args>
    requires(std::same_as<std::remove_const_t<Args>, State*> && ...)
  void add_inputs(Args... states) {
    (inputs_.push_back(states), ...);
  }

  template <typename... Args>
    requires(std::same_as<std::remove_const_t<Args>, State*> && ...)
  void add_outputs(Args... states) {
    (outputs_.push_back(states), ...);
  }

  std::vector<const State*> get_inputs() const { return inputs_; }
  std::vector<const State*> get_outputs() const { return outputs_; }

  friend class STN;
};
