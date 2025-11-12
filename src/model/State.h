#pragma once

#include <variant>

enum class StateType { INPUT, OUTPUT, NON_STORABLE, NORMAL };

class BaseState {
protected:
  size_t id_{};

  friend class STN;
  friend class State;
};

class InputState : public BaseState {};

class OutputState : public BaseState {};

class NonStorableState : public BaseState {};

class NormalState : public BaseState {
  size_t initial_stock_;
  size_t min_level_;
  size_t max_level_;

 public:
  NormalState(size_t initial_stock, size_t min_level, size_t max_level)
      : initial_stock_(initial_stock),
        min_level_(min_level),
        max_level_(max_level) {}
};

using StateVariant =
    std::variant<InputState, OutputState, NonStorableState, NormalState>;

class State : public StateVariant {
 public:
  using StateVariant::variant;

  size_t get_id() const {
    return std::visit([](const auto& state) { return state.id_; }, *this);
  }
};
