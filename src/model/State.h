#pragma once

#include <variant>

enum class StateType { INPUT, OUTPUT, NON_STORABLE, NORMAL };

class BaseState {
 protected:
  size_t id_{0};

  friend class STN;
  friend class State;
};

class InputState : public BaseState {};

class OutputState : public BaseState {};

class NonStorableState : public BaseState {};

class NormalState : public BaseState {
 public:
  size_t initial_stock;
  size_t min_level;
  size_t max_level;

  NormalState(size_t initial_stock, size_t min_level, size_t max_level)
      : initial_stock(initial_stock),
        min_level(min_level),
        max_level(max_level) {}
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
