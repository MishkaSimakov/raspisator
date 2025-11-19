#pragma once
#include <deque>

#include "State.h"
#include "Unit.h"

// State-Task Network
template <typename Field>
class STN {
  std::deque<Unit<Field>> units_;
  std::deque<Task<Field>> tasks_;
  std::deque<State> states_;

  enum class ConnectionType { PRODUCER, CONSUMER };

  auto get_connections_of(const State& state, ConnectionType type) const {
    std::vector<std::pair<const Task<Field>*, Field>> result;

    for (const Task<Field>& task : tasks_) {
      const auto& connections = type == ConnectionType::CONSUMER
                                    ? task.get_inputs()
                                    : task.get_outputs();

      for (const auto& [output, fraction] : connections) {
        if (output == &state) {
          result.emplace_back(&task, fraction);
          break;
        }
      }
    }

    return result;
  }

 public:
  Task<Field>* add(Task<Field> task) {
    task.set_id(tasks_.size() + 1);
    return &tasks_.emplace_back(task);
  }

  Unit<Field>* add(Unit<Field> unit) {
    unit.set_id(units_.size() + 1);
    return &units_.emplace_back(unit);
  }

  State* add(State state) {
    std::visit([id = states_.size() + 1](auto& state) { state.id_ = id; },
               state);
    return &states_.emplace_back(state);
  }

  const std::deque<State>& get_states() const { return states_; }
  const std::deque<Task<Field>>& get_tasks() const { return tasks_; }
  const std::deque<Unit<Field>>& get_units() const { return units_; }

  std::vector<std::pair<const Task<Field>*, Field>> get_producers_of(
      const State& state) const {
    return get_connections_of(state, ConnectionType::PRODUCER);
  }

  std::vector<std::pair<const Task<Field>*, Field>> get_consumers_of(
      const State& state) const {
    return get_connections_of(state, ConnectionType::CONSUMER);
  }

  //
  auto get_task_units(const Task<Field>& task) const {
    std::vector<std::pair<const Unit<Field>*, TaskOnUnitProperties<Field>>>
        result;

    for (const Unit<Field>& unit : units_) {
      for (const auto& [unit_task, props] : unit.get_tasks()) {
        if (unit_task == &task) {
          result.emplace_back(&unit, props);
        }
      }
    }

    return result;
  }
};
