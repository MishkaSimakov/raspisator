#pragma once
#include <deque>

#include "State.h"
#include "Unit.h"

// State-Task Network
class STN {
  std::deque<Unit> units_;
  std::deque<Task> tasks_;
  std::deque<State> states_;

 public:
  Task* add(Task task) {
    task.set_id(tasks_.size() + 1);
    return &tasks_.emplace_back(task);
  }

  Unit* add(Unit unit) {
    unit.set_id(units_.size() + 1);
    return &units_.emplace_back(unit);
  }

  State* add(State state) {
    std::visit([id = states_.size() + 1](auto& state) { state.id_ = id; },
               state);
    return &states_.emplace_back(state);
  }

  const std::deque<State>& get_states() const { return states_; }

  const std::deque<Task>& get_tasks() const { return tasks_; }

  //
  std::vector<std::pair<const Unit*, TaskOnUnitProperties>> get_task_units(
      const Task& task) const {
    std::vector<std::pair<const Unit*, TaskOnUnitProperties>> result;

    for (const Unit& unit : units_) {
      for (const auto& [unit_task, props] : unit.get_tasks()) {
        if (unit_task == &task) {
          result.emplace_back(&unit, props);
        }
      }
    }

    return result;
  }
};
