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

  Field input_fraction_of(const State<Field>& state) const {
    for (auto [input, fraction] : inputs_) {
      if (input->get_id() == state.get_id()) {
        return fraction;
      }
    }

    return 0;
  }

  Field output_fraction_of(const State<Field>& state) const {
    for (auto [output, fraction] : outputs_) {
      if (output->get_id() == state.get_id()) {
        return fraction;
      }
    }

    return 0;
  }

  friend class STN<Field>;
};
