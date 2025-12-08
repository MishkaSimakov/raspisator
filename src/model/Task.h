#pragma once

#include <vector>

#include "State.h"

template <typename Field>
class Task {
  size_t id_{0};

  // TODO: variable proportions for outputs
  std::vector<std::pair<const State<Field>*, Field>> inputs_;
  std::vector<std::pair<const State<Field>*, Field>> outputs_;

  void set_id(size_t id) { id_ = id; }

 public:
  size_t get_id() const { return id_; }

  void add_input(const State<Field>* state, Field fraction) {
    inputs_.emplace_back(state, std::move(fraction));
  }

  void add_output(const State<Field>* state, Field fraction) {
    outputs_.emplace_back(state, std::move(fraction));
  }

  const auto& get_inputs() const { return inputs_; }
  const auto& get_outputs() const { return outputs_; }

  friend class STN<Field>;
};
