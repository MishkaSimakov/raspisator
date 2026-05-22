#pragma once

#include <algorithm>
#include <vector>

#include "Task.h"

template <typename Field>
struct TaskOnUnitProperties {
  size_t batch_processing_time;
  Field batch_min_size;
  Field batch_max_size;
};

template <typename Field>
class Unit {
  size_t id_{0};

  std::vector<std::pair<const Task<Field>*, TaskOnUnitProperties<Field>>>
      tasks_;

  void set_id(size_t id) { id_ = id; }

 public:
  size_t get_id() const { return id_; }

  const auto& get_tasks() const { return tasks_; }

  const TaskOnUnitProperties<Field>* get_properties(
      const Task<Field>* task) const {
    for (const auto& p : tasks_) {
      if (p.first == task) {
        return &p.second;
      }
    }

    return nullptr;
  }

  void attach_task(const Task<Field>* task,
                   const TaskOnUnitProperties<Field>& properties) {
    tasks_.emplace_back(task, properties);
  }

  template <typename F>
  friend class STN;
};
