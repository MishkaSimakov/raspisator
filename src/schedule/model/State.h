#pragma once

#include <variant>

template <typename Field>
class STN;

template <typename Field>
class State;

template <typename Field>
class BaseState {
 protected:
  size_t id_{0};

  friend class STN<Field>;
  friend class State<Field>;
};

template <typename Field>
class InputState : public BaseState<Field> {
 public:
  Field initial_stock;

  explicit InputState(Field initial_stock) : initial_stock(initial_stock) {}
};

template <typename Field>
class OutputState : public BaseState<Field> {
 public:
  Field initial_stock;
  Field target;

  OutputState(Field initial_stock, Field target)
      : initial_stock(initial_stock), target(target) {}
};

template <typename Field>
class NonStorableState : public BaseState<Field> {};

template <typename Field>
class NormalState : public BaseState<Field> {
 public:
  Field initial_stock;
  Field min_level;
  Field max_level;

  NormalState(Field initial_stock, Field min_level, Field max_level)
      : initial_stock(initial_stock),
        min_level(min_level),
        max_level(max_level) {}
};

template <typename Field>
using StateVariant = std::variant<InputState<Field>, OutputState<Field>,
                                  NonStorableState<Field>, NormalState<Field>>;

template <typename Field>
class State : public StateVariant<Field> {
 public:
  using std::variant<InputState<Field>, OutputState<Field>,
                     NonStorableState<Field>, NormalState<Field>>::variant;

  size_t get_id() const {
    return std::visit([](const auto& state) { return state.id_; }, *this);
  }
};
