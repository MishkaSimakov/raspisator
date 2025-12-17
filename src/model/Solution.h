#pragma once

#include <algorithm>
#include <ranges>

#include "STN.h"
#include "linear/FieldTraits.h"
#include "utils/Variant.h"

template <typename Field>
struct Event {
  size_t unit_id;
  size_t task_id;
  Field batch_size;
  size_t time;
  bool is_start;
};

template <typename Field>
struct TaskInstance {
  size_t unit_id;
  size_t task_id;
  Field batch_size;
  size_t start_time;
};

template <typename Field>
class Solution {
  using Event = Event<Field>;
  using TaskInstance = TaskInstance<Field>;
  using Task = Task<Field>;
  using Unit = Unit<Field>;
  using TaskOnUnitProperties = TaskOnUnitProperties<Field>;

  const STN<Field>* stn;
  std::vector<TaskInstance> task_instances{};
  std::vector<Field> state_filled{};
  std::vector<bool> unit_free{};
  size_t current_time = 0;

  static bool compare_events(const Event& a, const Event& b) {
    if (a.time == b.time) {
      return a.is_start < b.is_start;
    }
    return a.time < b.time;
  }

  size_t get_instance_end_time(const TaskInstance& instance) const {
    const Unit* unit = stn->get_unit_by_id(instance.unit_id);
    const Task* task = stn->get_task_by_id(instance.task_id);
    const TaskOnUnitProperties* prop = unit->get_properties(task);

    return instance.start_time + prop->batch_processing_time;
  }

  bool check_task_instance_valid(TaskInstance ti) {
    // this checks that task_instance is valid in this stn
    if (ti.unit_id >= stn->get_units().size()) {
      std::println(
          "Unit index out of range in TaskInstance(unit_id: {}, task_id: {}, "
          "start_time: {}). Unit count: {}",
          ti.unit_id, ti.task_id, ti.start_time, stn->get_units().size());
      return false;
    }

    if (ti.task_id >= stn->get_tasks().size()) {
      std::println(
          "Task index out of range in TaskInstance(unit_id: {}, task_id: {}, "
          "start_time: {}). Task count {}",
          ti.unit_id, ti.task_id, ti.start_time, stn->get_tasks().size());
      return false;
    }

    auto unit_tasks = stn->get_unit_by_id(ti.unit_id)->get_tasks();
    const Task* task = stn->get_task_by_id(ti.task_id);
    for (const auto& p : unit_tasks) {
      if (p.first == task) {
        return true;
      }
    }
    std::println(
        "Task {} cannot be processed on unit {} in TaskInstance(unit_id: {}, "
        "task_id: {}, start_time: {})",
        ti.task_id, ti.unit_id, ti.unit_id, ti.task_id, ti.start_time);
    return false;
  }

  std::optional<std::vector<Event>> generate_events() {
    // this generates sorted TimeInstance vector for further work
    std::vector<Event> out;

    for (TaskInstance ti : task_instances) {
      out.emplace_back(ti.unit_id, ti.task_id, ti.batch_size, ti.start_time,
                       true);

      const Unit* unit = stn->get_unit_by_id(ti.unit_id);
      const Task* task = stn->get_task_by_id(ti.task_id);
      const TaskOnUnitProperties* prop = unit->get_properties(task);
      if (!prop) {
        std::println(
            "Task {} cannot be processed on unit {} in TaskInstance(unit_id: "
            "{}, task_id: {}, start_time: {})",
            ti.task_id, ti.unit_id, ti.unit_id, ti.task_id, ti.start_time);
        return std::nullopt;
      }

      out.emplace_back(ti.unit_id, ti.task_id, ti.batch_size,
                       ti.start_time + prop->batch_processing_time, false);
    }

    std::sort(out.begin(), out.end(), compare_events);
    return std::optional(out);
  }

  bool check_state_bounds(Event e) {
    // checks that states are not overflowed

    // fmt::println("event time: {}, current_time: {}", e.time, current_time);
    // if e happens at current time, overflowed states might be relaxed in e, so
    // no need to worry
    if (e.time == current_time) {
      return true;
    }

    for (const State<Field>& s : stn->get_states()) {
      Field filled = state_filled[s.get_id()];
      bool result = std::visit(
          Overload{[&filled](NonStorableState<Field>) {
                     return !FieldTraits<Field>::is_nonzero(filled);
                   },
                   [&filled](NormalState<Field> state) {
                     // fmt::println("min: {}, max: {}, filled: {}",
                     // state.min_level, state.max_level, filled);
                     return !FieldTraits<Field>::is_strictly_positive(
                                state.min_level - filled) &&
                            !FieldTraits<Field>::is_strictly_positive(
                                filled - state.max_level);
                   },
                   [&filled](const auto&) { return true; }},
          s);

      if (!result) {
        std::println("State {} out of bounds on time {}", s.get_id(), e.time);
        return false;
      }
    }

    return true;
  }

  void increase_state(const State<Field>* s, const Field& value) {
    // fmt::println("current_value: {}, add: {}", state_filled[s->get_id()],
    // value);
    state_filled[s->get_id()] += value;
  }

  bool apply_event(Event e) {
    const Unit* unit = stn->get_unit_by_id(e.unit_id);
    const Task* task = stn->get_task_by_id(e.task_id);
    const TaskOnUnitProperties* prop = unit->get_properties(task);
    // fmt::println("e batch_size: {}, prop_min: {}, prop_max: {}",
    // e.batch_size, prop->batch_min_size, prop->batch_max_size);
    if (FieldTraits<Field>::is_strictly_positive(e.batch_size -
                                                 prop->batch_max_size) ||
        FieldTraits<Field>::is_strictly_negative(e.batch_size <
                                                 prop->batch_min_size)) {
      std::println(
          "event batch size {} out of bounds for task {} on time {}. Bounds: "
          "{}, {}",
          e.batch_size, task->get_id(), e.time, prop->batch_min_size,
          prop->batch_max_size);
      return false;
    }
    if (e.is_start) {
      if (!unit_free[e.unit_id]) {
        std::println("Attempt to access busy unit {} by task {} on time {}",
                     e.unit_id, e.task_id, e.time);
        return false;
      }
      unit_free[e.unit_id] = false;

      auto inputs = task->get_inputs();
      for (auto [s, proportion] : inputs) {
        Field needed_amount = e.batch_size * proportion;
        increase_state(s, -needed_amount);
      }

    } else {
      if (unit_free[e.unit_id]) {
        std::println(
            "Attempt to free already free unit {} by task {} on time {}",
            e.unit_id, e.task_id, e.time);
        return false;
      }
      unit_free[e.unit_id] = true;

      auto outputs = task->get_outputs();
      for (auto [s, proportion] : outputs) {
        Field out_amount = e.batch_size * proportion;
        increase_state(s, out_amount);
      }
    }

    current_time = e.time;

    return true;
  }

  bool check_output_states() {
    for (const State<Field>& s : stn->get_states()) {
      Field filled = state_filled[s.get_id()];
      bool result = std::visit(
          Overload{[&filled](OutputState<Field> state) {
                     return !FieldTraits<Field>::is_strictly_negative(
                         filled - state.target);
                   },
                   [&filled](const auto&) { return true; }},
          s);
      if (!result) {
        std::println("target not acquired for output state {}", s.get_id());
        return false;
      }
    }
    return true;
  }

 public:
  explicit Solution(const STN<Field>* stn)
      : stn(stn),
        state_filled(stn->get_states().size()),
        unit_free(stn->get_units().size(), true) {
    for (auto s : stn->get_states()) {
      state_filled[s.get_id()] = std::visit(
          Overload{[](const auto& state) { return Field(state.initial_stock); },
                   [](NonStorableState<Field>) -> Field { return Field(0); }},
          s);
    }
  }

  void add_instance(TaskInstance instance) {
    task_instances.push_back(instance);
  }

  std::vector<TaskInstance> get_unit_tasks(const Unit& unit) const {
    std::vector<TaskInstance> result;

    for (const auto& instance : task_instances) {
      if (instance.unit_id == unit.get_id()) {
        result.push_back(instance);
      }
    }

    std::ranges::sort(result, {}, [](const TaskInstance& instance) {
      return instance.start_time;
    });

    return result;
  }

  bool check() {
    for (TaskInstance ti : task_instances) {
      if (!check_task_instance_valid(ti)) {
        return false;
      }
    }

    auto events_opt = generate_events();
    if (events_opt == std::nullopt) {
      return false;
    }
    std::vector<Event> events = events_opt.value();

    for (Event e : events) {
      if (!check_state_bounds(e)) {
        return false;
      }
      if (!apply_event(e)) {
        return false;
      }
    }

    if (!check_state_bounds({0, 0, 0, events.back().time + 1, 0})) {
      return false;
    }
    return check_output_states();
  }

  size_t get_total_time() const {
    return std::ranges::max(
        task_instances |
        std::views::transform([this](const TaskInstance& instance) {
          return get_instance_end_time(instance);
        }));
  }

  void to_graphviz(std::ostream& os) const {
    os << "digraph {\n";
    os << "graph [ pad=\"0.5\", nodesep=\"0.5\", ranksep=\"2\" ]\n";
    os << "node  [ shape=plain]\n";

    os << "Foo [label=<\n";
    os << "<table border=\"0\" cellborder=\"1\" cellspacing=\"0\">\n";

    size_t max_time = get_total_time();

    os << "<tr><td></td>";
    for (size_t i = 0; i < max_time; ++i) {
      os << std::format("<td width=\"75\">{}</td>", i);
    }
    os << "</tr>\n";

    for (const auto& unit : stn->get_units()) {
      os << "<tr>";

      os << std::format("<td>unit_{}</td>", unit.get_id());
      std::vector<TaskInstance> unit_task_instances;

      for (const auto& ti : task_instances) {
        if (ti.unit_id == unit.get_id()) {
          unit_task_instances.push_back(ti);
        }
      }

      std::ranges::sort(
          unit_task_instances, {},
          [](const TaskInstance& instance) { return instance.start_time; });

      auto itr = unit_task_instances.begin();
      for (size_t t = 0; t < max_time;) {
        if (itr != unit_task_instances.end() && itr->start_time == t) {
          size_t duration = get_instance_end_time(*itr) - itr->start_time;
          os << std::format(
              "<td colspan=\"{}\" bgcolor=\"grey\">{} ({:.1f})</td>", duration,
              itr->task_id, itr->batch_size);

          t += duration;
          ++itr;
        } else {
          os << "<td></td>";

          ++t;
        }
      }

      os << "</tr>\n";
    }

    os << "</table>>];\n";
    os << "}\n";
  }
};
