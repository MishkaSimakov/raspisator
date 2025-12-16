#pragma once

#include <format>
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

template <typename Field>
struct ProblemEncoding {
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
};

// MILP encoding for chemical batch processing scheduling problem
// This encoding from F. Blomer & H.-O. Gunther's paper uses uniform time
// discretization
template <typename Field>
ProblemEncoding<Field> to_uniform_time_milp(const STN<Field>& problem,
                                            size_t max_periods) {
  const Field kGlobalMaxStock = 2000;

  MILPProblem<Field> builder;

  // decision variables
  auto makespan =
      builder.new_variable("MS", VariableType::INTEGER, 0, max_periods - 1);
  std::unordered_map<std::tuple<const Unit<Field>*, const Task<Field>*, size_t>,
                     Variable<Field>>
      quantities;
  std::unordered_map<std::tuple<const Unit<Field>*, const Task<Field>*, size_t>,
                     Variable<Field>>
      starts;
  std::unordered_map<std::pair<const State<Field>*, size_t>, Variable<Field>>
      stocks;

  for (const auto& unit : problem.get_units()) {
    for (const auto& [task, info] : unit.get_tasks()) {
      for (size_t t = 0; t < max_periods; ++t) {
        quantities.emplace(
            std::tuple{&unit, task, t},
            builder.new_variable(
                std::format("Q({}, {}, {})", unit.get_id(), task->get_id(), t),
                VariableType::REAL, 0, info.batch_max_size));
        starts.emplace(
            std::tuple{&unit, task, t},
            builder.new_variable(
                std::format("x({}, {}, {})", unit.get_id(), task->get_id(), t),
                VariableType::INTEGER, 0, 1));
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

  return ProblemEncoding<Field>{builder, makespan, std::move(quantities),
                                std::move(starts), std::move(stocks)};
}

template <typename Field>
std::optional<Solution<Field>> get_schedule(const STN<Field>& problem,
                                            size_t max_periods) {
  auto encoding = to_uniform_time_milp(problem, max_periods);

  std::println("before constraints {}, variables {}",
             encoding.builder.constraints.size(),
             encoding.builder.variables.size());

  auto optimizer = FullOptimizer<Field>(false);
  auto optimized_problem = optimizer.apply(encoding.builder);

  std::println("after  constraints {}, variables {}",
               optimized_problem.constraints.size(),
               optimized_problem.variables.size());

  // solve MILP problem
  auto matrices = to_matrices(optimized_problem);

  auto settings = BranchAndBoundSettings<Field>{
      .max_nodes = 50'000,
      .strong_branching_max_iterations_factor = std::nullopt};

  auto solver = FullStrongBranchingBranchAndBound<Field>(
      matrices.A, matrices.b, matrices.c, matrices.lower, matrices.upper,
      matrices.variables, settings);
  auto milp_solution = solver.solve();

  if (!std::holds_alternative<FiniteMILPSolution<Field>>(milp_solution)) {
    return std::nullopt;
  }

  auto point = optimizer.inverse(
      std::get<FiniteMILPSolution<Field>>(milp_solution).point);

  Solution solution(&problem);

  for (const auto& unit : problem.get_units()) {
    for (size_t t = 0; t < max_periods; ++t) {
      for (const auto* task : unit.get_tasks() | std::views::keys) {
        Field x = encoding.builder.extract_variable(
            encoding.starts.at({&unit, task, t}), point);
        Field Q = encoding.builder.extract_variable(
            encoding.quantities.at({&unit, task, t}), point);

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
