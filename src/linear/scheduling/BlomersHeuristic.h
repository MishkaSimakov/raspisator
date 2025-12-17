#pragma once

#include <format>
#include <functional>
#include <ranges>

#include "linear/bb/FullStrongBranching.h"
#include "linear/bb/Settings.h"
#include "linear/model/MILP.h"
#include "linear/problem/MILPProblem.h"
#include "linear/problem/ToMatrices.h"
#include "linear/problem/optimization/FullOptimizer.h"
#include "model/STN.h"
#include "model/Solution.h"
#include "utils/Hashers.h"
#include "utils/Variant.h"

using TimeGrid = std::unordered_map<size_t, std::function<bool(size_t)>>;

// MILP encoding for chemical batch processing scheduling problem
// This encoding from F. Blomer & H.-O. Gunther's paper uses uniform time
// discretization
template <typename Field>
struct TimeGridModel {
  MILPProblem<Field> builder;

  Variable<Field> makespan;
  std::unordered_map<std::tuple<const Unit<Field>*, const Task<Field>*, size_t>,
                     Variable<Field>>
      quantities;
  std::unordered_map<std::tuple<const Unit<Field>*, const Task<Field>*, size_t>,
                     Variable<Field>>
      starts;
  std::unordered_map<std::pair<const State<Field>*, size_t>, Variable<Field>>
      stocks;

  TimeGridModel(const STN<Field>& problem, size_t max_periods,
                const TimeGrid& time_grid);
};

template <typename Field>
TimeGridModel<Field>::TimeGridModel(const STN<Field>& problem,
                                    size_t max_periods,
                                    const TimeGrid& time_grid)
    : makespan(builder.new_variable("MS", VariableType::INTEGER, 0,
                                    max_periods - 1)) {
  const Field kGlobalMaxStock = 2000;

  // decision variables
  // makespan is already defined

  for (const auto& unit : problem.get_units()) {
    for (const auto& [task, info] : unit.get_tasks()) {
      for (size_t t = 0; t < max_periods; ++t) {
        quantities.emplace(
            std::tuple{&unit, task, t},
            builder.new_variable(
                std::format("Q({}, {}, {})", unit.get_id(), task->get_id(), t),
                VariableType::REAL, 0, info.batch_max_size));

        bool can_start = time_grid.at(unit.get_id())(t);

        // optimizers should do the rest of the work
        starts.emplace(
            std::tuple{&unit, task, t},
            builder.new_variable(
                std::format("x({}, {}, {})", unit.get_id(), task->get_id(), t),
                VariableType::INTEGER, 0, can_start ? 1 : 0));
      }
    }
  }

  for (const State<Field>& state : problem.get_states()) {
    if (!std::holds_alternative<NonStorableState<Field>>(state)) {
      Field max_stock =
          std::visit(Overload{
                         [](const NormalState<Field>& state) -> Field {
                           return state.max_level;
                         },
                         [kGlobalMaxStock](const auto&) -> Field {
                           return kGlobalMaxStock;
                         },
                     },
                     state);

      for (size_t t = 0; t < max_periods; ++t) {
        stocks.emplace(
            std::pair{&state, t},
            builder.new_variable(std::format("p({}, {})", state.get_id(), t),
                                 VariableType::REAL, 0, max_stock));
      }
    }
  }

  // makespan
  for (const auto& unit : problem.get_units()) {
    for (const auto& [task, props] : unit.get_tasks()) {
      for (size_t t = 0; t < max_periods; ++t) {
        Field finish_time(t + props.batch_processing_time);
        builder.add_constraint(makespan >=
                               finish_time * starts.at({&unit, task, t}));
      }
    }
  }

  // batch size limits
  for (const auto& unit : problem.get_units()) {
    for (const auto& [task, props] : unit.get_tasks()) {
      for (size_t t = 0; t < max_periods; ++t) {
        builder.add_constraint(Field(props.batch_min_size) *
                                   starts.at({&unit, task, t}) <=
                               quantities.at({&unit, task, t}));

        builder.add_constraint(quantities.at({&unit, task, t}) <=
                               Field(props.batch_max_size) *
                                   starts.at({&unit, task, t}));
      }
    }
  }

  // stock balance
  for (const State<Field>& state : problem.get_states()) {
    if (std::holds_alternative<NonStorableState<Field>>(state)) {
      continue;
    }

    for (size_t t = 0; t < max_periods; ++t) {
      auto new_stock =
          t != 0 ? Expression(stocks.at({&state, t - 1}))
                 : std::visit(Overload{[](const auto& state) {
                                         return Field(state.initial_stock);
                                       },
                                       [](NonStorableState<Field>) -> Field {
                                         std::unreachable();
                                       }},
                              state);

      for (const auto& [task, fraction] : problem.get_producers_of(state)) {
        Expression<Field> task_production(0);

        for (const auto& [unit, props] : problem.get_task_units(*task)) {
          if (t >= props.batch_processing_time) {
            task_production +=
                quantities.at({unit, task, t - props.batch_processing_time});
          }
        }

        new_stock += task_production * fraction;
      }

      for (const auto& [task, fraction] : problem.get_consumers_of(state)) {
        Expression<Field> task_consumption(0);

        for (const auto& [unit, props] : problem.get_task_units(*task)) {
          task_consumption += quantities.at({unit, task, t});
        }

        new_stock -= task_consumption * fraction;
      }

      // TODO: external demand + external supply

      builder.add_constraint(new_stock == stocks.at({&state, t}));
    }
  }

  // stock limits are incorporated into bounds

  // production of non-storable goods
  for (const State<Field>& state : problem.get_states()) {
    if (!std::holds_alternative<NonStorableState<Field>>(state)) {
      continue;
    }

    for (size_t t = 0; t < max_periods; ++t) {
      Expression<Field> production(0);

      for (const auto& [task, fraction] : problem.get_producers_of(state)) {
        Expression<Field> task_production(0);

        for (const auto& [unit, props] : problem.get_task_units(*task)) {
          if (t >= props.batch_processing_time) {
            task_production +=
                quantities.at({unit, task, t - props.batch_processing_time});
          }
        }

        production += task_production * fraction;
      }

      Expression<Field> consumption(0);

      for (const auto& [task, fraction] : problem.get_consumers_of(state)) {
        Expression<Field> task_consumption(0);

        for (const auto& [unit, props] : problem.get_task_units(*task)) {
          task_consumption += quantities.at({unit, task, t});
        }

        consumption += task_consumption * fraction;
      }

      builder.add_constraint(consumption == production);
    }
  }

  // assigning batches to production units
  for (const auto& unit : problem.get_units()) {
    for (size_t t = 0; t < max_periods; ++t) {
      Expression<Field> running_batches_cnt(0);

      for (const auto& [task, props] : unit.get_tasks()) {
        for (size_t t2 = std::max(t, props.batch_processing_time) -
                         props.batch_processing_time;
             t2 < t; ++t2) {
          running_batches_cnt += starts.at({&unit, task, t2});
        }
      }

      builder.add_constraint(running_batches_cnt <= Expression<Field>(1));
    }
  }

  // variables domain are incorporated into bounds

  // desired amounts
  for (const State<Field>& state : problem.get_states()) {
    if (!std::holds_alternative<OutputState<Field>>(state)) {
      continue;
    }

    Expression<Field> desired_amount(
        std::get<OutputState<Field>>(state).target);

    builder.add_constraint(stocks.at({&state, max_periods - 1}) >=
                           desired_amount);
  }

  // objective
  builder.set_objective(-makespan);
}

// clang-format off
template <typename Field>
struct LeftShiftModel {
  MILPProblem<Field> builder;

  // decision variables
  Variable<Field> makespan;
  std::unordered_map<std::pair<size_t, const Unit<Field>*>, Variable<Field>> starts;
  std::unordered_map<std::pair<size_t, const Unit<Field>*>, Variable<Field>> finishes;
  std::unordered_map<std::pair<const State<Field>*, size_t>, Variable<Field>> stocks;

  // auxiliary flags
  std::unordered_map<std::tuple<size_t, const Unit<Field>*, size_t>, Variable<Field>> past_starts;
  std::unordered_map<std::tuple<size_t, const Unit<Field>*, size_t>, Variable<Field>> past_finishes;

  LeftShiftModel(const STN<Field>& problem, size_t max_periods,
                 const Solution<Field>& solution);
};

template <typename Field>
LeftShiftModel<Field>::LeftShiftModel(const STN<Field>& problem,
                                      size_t max_periods,
                                      const Solution<Field>& solution)
    : makespan(builder.new_variable("MS", VariableType::INTEGER, 0, max_periods - 1)) {
  const Field kGlobalMaxStock = 2000;

  // decision variables
  // makespan is already defined

  for (const auto& unit : problem.get_units()) {
    size_t batches_count = solution.get_unit_tasks(unit).size();

    for (size_t i = 0; i < batches_count; ++i) {
      starts.emplace(std::pair{i, &unit}, builder.new_variable(std::format("S({}, {})", i, unit.get_id()), VariableType::REAL, 0, max_periods - 1));
      finishes.emplace(std::pair{i, &unit}, builder.new_variable(std::format("F({}, {})", i, unit.get_id()), VariableType::REAL, 0, max_periods - 1));

      for (size_t t = 0; t < max_periods; ++t) {
        past_starts.emplace(std::tuple{i, &unit, t}, builder.new_variable(std::format("s({}, {}, {})", i, unit.get_id(), t), VariableType::INTEGER, 0, 1));
        past_finishes.emplace(std::tuple{i, &unit, t}, builder.new_variable(std::format("f({}, {}, {})", i, unit.get_id(), t), VariableType::INTEGER, 0, 1));
      }
    }
  }

  for (const State<Field>& state : problem.get_states()) {
    Field max_stock = std::visit(
        Overload{
            [](const NormalState<Field>& state) -> Field {
              return state.max_level;
            },
            [](NonStorableState<Field>) -> Field { return 0; },
            [kGlobalMaxStock](const auto&) -> Field { return kGlobalMaxStock; },
        },
        state);

    for (size_t t = 0; t < max_periods; ++t) {
      stocks.emplace(std::pair{&state, t}, builder.new_variable(std::format("p({}, {})", state.get_id(), t), VariableType::REAL, 0, max_stock));
    }
  }

  // objective
  builder.set_objective(-makespan);

  // makespan
  for (const auto& unit : problem.get_units()) {
    size_t batches_count = solution.get_unit_tasks(unit).size();

    if (batches_count > 0) {
      builder.add_constraint(makespan >= finishes.at({batches_count - 1, &unit}));
    }
  }

  // start and finish times
  for (const auto& unit : problem.get_units()) {
    size_t batches_count = solution.get_unit_tasks(unit).size();

    for (size_t i = 0; i < batches_count; ++i) {
      if (i >= 1) {
        builder.add_constraint(starts.at({i, &unit}) >= finishes.at({i - 1, &unit}));
      }

      size_t task_id = solution.get_unit_tasks(unit)[i].task_id;
      const Task<Field>* task = problem.get_task_by_id(task_id);
      size_t processing_time = unit.get_properties(task)->batch_processing_time;

      builder.add_constraint(finishes.at({i, &unit}) == starts.at({i, &unit}) + Expression<Field>(processing_time));
    }
  }

  // Heaviside function: start & finish events
  for (const auto& unit : problem.get_units()) {
    size_t batches_count = solution.get_unit_tasks(unit).size();

    for (size_t i = 0; i < batches_count; ++i) {
      for (size_t t = 0; t < max_periods; ++t) {
        builder.add_constraint(past_starts.at({i, &unit, t}) <= (Expression<Field>(t) - starts.at({i, &unit})) / Field(max_periods) + Expression<Field>(1));
        builder.add_constraint(past_starts.at({i, &unit, t}) >= (Expression<Field>(t + 1) - starts.at({i, &unit})) / Field(max_periods));

        builder.add_constraint(past_finishes.at({i, &unit, t}) <= (Expression<Field>(t) - finishes.at({i, &unit})) / Field(max_periods) + Expression<Field>(1));
        builder.add_constraint(past_finishes.at({i, &unit, t}) >= (Expression<Field>(t + 1) - finishes.at({i, &unit})) / Field(max_periods));
      }
    }
  }

  // stock balance
  for (const auto& state : problem.get_states()) {
    Field initial_stock = std::visit(Overload{
        [](const auto& state) ->Field { return state.initial_stock; },
        [](NonStorableState<Field>) -> Field { return 0; },
    }, state);

    for (size_t t = 0; t < max_periods; ++t) {
      Expression<Field> income;

      for (const auto& unit : problem.get_units()) {
        for (size_t i = 0; i < solution.get_unit_tasks(unit).size(); ++i) {
          TaskInstance<Field> instance = solution.get_unit_tasks(unit)[i];
          const Task<Field>* task = problem.get_task_by_id(instance.task_id);

          Field input_size = task->input_fraction_of(state) * instance.batch_size;
          Field output_size = task->output_fraction_of(state) * instance.batch_size;

          income += output_size * past_finishes.at({i, &unit, t}) - input_size * past_starts.at({i, &unit, t});
        }
      }

      auto new_stock = Expression<Field>(initial_stock) + income;
      builder.add_constraint(stocks.at({&state, t}) == new_stock);
    }
  }

  // Stock limits are incorporated into variables bounds
}
// clang-format on

template <typename Field>
std::optional<Solution<Field>> apply_time_grid_model(
    const STN<Field>& problem, size_t max_periods, const TimeGrid& time_grid) {
  TimeGridModel<Field> model(problem, max_periods, time_grid);

  std::println("    before: constraints {}, variables {}, average gap: {}",
               model.builder.constraints.size(), model.builder.variables.size(),
               model.builder.average_boundary_gap());

  auto optimizer = FullOptimizer<Field>();
  auto optimized_problem = optimizer.apply(model.builder);

  std::println("    after: constraints {}, variables {}, average gap: {}",
               optimized_problem.constraints.size(),
               optimized_problem.variables.size(),
               optimized_problem.average_boundary_gap());

  // solve MILP problem
  auto matrices = to_matrices(optimized_problem);

  auto settings = BranchAndBoundSettings<Field>{
      .max_nodes = 50'000,
      .strong_branching_max_iterations_factor = 100,
      .strong_branching_min_iterations_limit = 10'000,
  };

  auto solver = FullStrongBranchingBranchAndBound<Field>(
      matrices.A, matrices.b, matrices.c, matrices.lower, matrices.upper,
      matrices.variables, settings);
  auto run_result = solver.solve();

  if (!std::holds_alternative<FiniteMILPSolution<Field>>(run_result.solution)) {
    return std::nullopt;
  }

  std::println("    nodes count: {}, average simplex iterations: {}",
               run_result.nodes_count, run_result.average_simplex_iterations);

  auto point = optimizer.inverse(
      std::get<FiniteMILPSolution<Field>>(run_result.solution).point);

  Solution solution(&problem);

  for (const auto& unit : problem.get_units()) {
    for (size_t t = 0; t < max_periods; ++t) {
      for (const auto* task : unit.get_tasks() | std::views::keys) {
        Field x = model.builder.extract_variable(
            model.starts.at({&unit, task, t}), point);
        Field Q = model.builder.extract_variable(
            model.quantities.at({&unit, task, t}), point);

        if (FieldTraits<Field>::is_strictly_positive(x) &&
            FieldTraits<Field>::is_strictly_positive(Q)) {
          solution.add_instance(
              TaskInstance{unit.get_id(), task->get_id(), Q, t});
        }
      }
    }
  }

  return solution;
}

template <typename Field>
std::optional<Solution<Field>> apply_left_shift_model(
    const STN<Field>& problem, size_t max_periods,
    const Solution<Field>& solution) {
  LeftShiftModel<Field> model(problem, max_periods, solution);

  std::println("    before: constraints {}, variables {}, average gap: {}",
               model.builder.constraints.size(), model.builder.variables.size(),
               model.builder.average_boundary_gap());

  auto optimizer = FullOptimizer<Field>();
  auto optimized_problem = optimizer.apply(model.builder);

  std::println("    after: constraints {}, variables {}, average gap: {}",
               optimized_problem.constraints.size(),
               optimized_problem.variables.size(),
               optimized_problem.average_boundary_gap());

  // solve MILP problem
  auto matrices = to_matrices(optimized_problem);

  auto settings = BranchAndBoundSettings<Field>{
      .max_nodes = 50'000,
      .strong_branching_max_iterations_factor = 100,
      .strong_branching_min_iterations_limit = 10'000,
  };

  auto solver = FullStrongBranchingBranchAndBound<Field>(
      matrices.A, matrices.b, matrices.c, matrices.lower, matrices.upper,
      matrices.variables, settings);
  auto run_result = solver.solve();

  if (!std::holds_alternative<FiniteMILPSolution<Field>>(run_result.solution)) {
    return std::nullopt;
  }

  std::println("    nodes count: {}, average simplex iterations: {}",
               run_result.nodes_count, run_result.average_simplex_iterations);

  auto point = optimizer.inverse(
      std::get<FiniteMILPSolution<Field>>(run_result.solution).point);

  Solution refined_solution(&problem);

  for (const auto& unit : problem.get_units()) {
    const auto& batches = solution.get_unit_tasks(unit);

    for (size_t i = 0; i < batches.size(); ++i) {
      size_t start_time = 0;

      for (size_t t = 0; t < max_periods; ++t) {
        Field value = model.builder.extract_variable(
            model.past_starts.at({i, &unit, t}), point);

        if (FieldTraits<Field>::is_strictly_positive(value)) {
          start_time = t;
          break;
        }
      }

      TaskInstance<Field> instance = batches[i];
      instance.start_time = start_time;

      refined_solution.add_instance(instance);
    }
  }

  return refined_solution;
}

template <typename Field>
std::optional<Solution<Field>> blomer_heuristic_model(
    const STN<Field>& problem, size_t max_periods, const TimeGrid& time_grid) {
  std::println("  solving rough:");
  auto rough_solution = apply_time_grid_model(problem, max_periods, time_grid);

  if (!rough_solution) {
    return std::nullopt;
  }

  if (!rough_solution->check()) {
    throw std::runtime_error("Rough solution is invalid.");
  }

  std::println("    periods: {}", rough_solution->get_total_time());
  // rough_solution->to_graphviz(std::cout);

  std::println("  solving left shift:");
  auto refined_solution =
      apply_left_shift_model(problem, max_periods, *rough_solution);

  std::println("    periods: {}", refined_solution->get_total_time());

  return refined_solution;
}
