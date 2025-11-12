#pragma once

#include <vector>

#include "Task.h"

struct TaskOnUnitProperties {
  size_t batch_processing_time;
  size_t batch_min_size;
  size_t batch_max_size;
};

class Unit {
  size_t id_{};

  std::vector<std::pair<const Task*, TaskOnUnitProperties>> tasks_;

  void set_id(size_t id) { id_ = id; }

 public:
  size_t get_id() const { return id_; }

  const std::vector<std::pair<const Task*, TaskOnUnitProperties>>& get_tasks()
      const {
    return tasks_;
  }

  void attach_task(const Task* task, const TaskOnUnitProperties& properties) {
    tasks_.emplace_back(task, properties);
  }

  friend class STN;
};
