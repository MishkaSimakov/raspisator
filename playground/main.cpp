#include <iostream>
#include <print>

#include "linear/BigInteger.h"
#include "linear/BranchAndBound.h"
#include "linear/MPS.h"
#include "linear/SimplexMethod.h"
#include "linear/builder/ProblemBuilder.h"
#include "problems/Dwarf.h"
#include "utils/Hashers.h"

using Field = Rational;

int main() {
  auto problem = generate_problem<Field>();

  size_t H = 3;  // number of periods

  ProblemBuilder<Field> builder;

  // decision variables
  auto makespan = builder.new_variable("MS", VariableType::INTEGER);
  std::unordered_map<std::tuple<const Unit<Field>*, const Task<Field>*, size_t>,
                     Variable<Field>>
      quantities;
  std::unordered_map<std::tuple<const Unit<Field>*, const Task<Field>*, size_t>,
                     Variable<Field>>
      starts;
  std::unordered_map<std::pair<const State*, size_t>, Variable<Field>> stocks;

  for (const auto& unit : problem.get_units()) {
    for (const auto* task : unit.get_tasks() | std::views::keys) {
      for (size_t t = 0; t < H; ++t) {
        quantities.emplace(
            std::tuple{&unit, task, t},
            builder.new_variable(
                std::format("Q({}, {}, {})", unit.get_id(), task->get_id(), t),
                VariableType::INTEGER));
        starts.emplace(
            std::tuple{&unit, task, t},
            builder.new_variable(
                std::format("x({}, {}, {})", unit.get_id(), task->get_id(), t),
                VariableType::INTEGER));
      }
    }
  }

  for (const State& state : problem.get_states()) {
    if (!std::holds_alternative<NonStorableState>(state)) {
      for (size_t t = 0; t < H; ++t) {
        stocks.emplace(
            std::pair{&state, t},
            builder.new_variable(fmt::format("p({}, {})", state.get_id(), t),
                                 VariableType::INTEGER));
      }
    }
  }

  // makespan
  for (const auto& unit : problem.get_units()) {
    for (const auto& [task, props] : unit.get_tasks()) {
      for (size_t t = 0; t < H; ++t) {
        Field finish_time(t + props.batch_processing_time - 1);
        builder.add_constraint(makespan >=
                               finish_time * starts.at({&unit, task, t}));
      }
    }
  }

  // batch size limits
  for (const auto& unit : problem.get_units()) {
    for (const auto& [task, props] : unit.get_tasks()) {
      for (size_t t = 0; t < H; ++t) {
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
  for (const State& state : problem.get_states()) {
    if (std::holds_alternative<NonStorableState>(state)) {
      continue;
    }

    for (size_t t = 0; t < H; ++t) {
      auto new_stock =
          t != 0 ? Expression(stocks.at({&state, t - 1}))
                 : std::visit(Overload{[](const auto& state) {
                                         return Field(state.initial_stock);
                                       },
                                       [](NonStorableState) -> Field {
                                         std::unreachable();
                                       }},
                              state);

      for (const auto* task :
           problem.get_producers_of(state) | std::views::elements<0>) {
        for (const auto& [unit, props] : problem.get_task_units(*task)) {
          if (t >= props.batch_processing_time) {
            new_stock +=
                quantities.at({unit, task, t - props.batch_processing_time});
          }
        }
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

  // stock limits
  for (size_t t = 0; t < H; ++t) {
    for (const State& state : problem.get_states()) {
      if (!std::holds_alternative<NormalState>(state)) {
        continue;
      }

      Expression<Field> max_stock(std::get<NormalState>(state).max_level);
      builder.add_constraint(stocks.at({&state, t}) <= max_stock);
    }
  }

  // production of non-storable goods
  for (const State& state : problem.get_states()) {
    if (!std::holds_alternative<NonStorableState>(state)) {
      continue;
    }

    for (size_t t = 0; t < H; ++t) {
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
    for (size_t t = 0; t < H; ++t) {
      Expression<Field> running_batches_cnt(0);

      for (const auto& [task, props] : unit.get_tasks()) {
        if (props.batch_processing_time >= t) {
          for (size_t t2 = t - props.batch_processing_time; t2 < t; ++t2) {
            running_batches_cnt += starts.at({&unit, task, t2});
          }
        }
      }

      builder.add_constraint(running_batches_cnt <= Expression<Field>(1));
    }
  }

  // variables domain
  for (const auto& unit : problem.get_units()) {
    for (const auto* task : unit.get_tasks() | std::views::keys) {
      for (size_t t = 0; t < H; ++t) {
        builder.add_constraint(starts.at({&unit, task, t}) <=
                               Expression<Field>(1));
      }
    }
  }

  // desired amounts
  for (const State& state : problem.get_states()) {
    if (!std::holds_alternative<OutputState>(state)) {
      continue;
    }

    Expression<Field> desired_amount(
        std::get<OutputState>(state).desired_amount);

    builder.add_constraint(stocks.at({&state, H - 1}) >= desired_amount);
  }

  // objective
  builder.set_objective(-makespan);

  std::cout << builder << std::endl;

  // solve MILP problem
  auto milp_problem = builder.get_problem();

  auto solution =
      BranchAndBound<Field, SimplexMethod<Field>>(milp_problem).solve();

  if (std::holds_alternative<NoFiniteSolution>(solution)) {
    std::println("No finite solution.");
  } else {
    auto finite_solution = std::get<FiniteMILPSolution<Field>>(solution);

    std::println("Finish production in {} time units.\n",
                 -finite_solution.value);

    auto point = finite_solution.point;

    for (const auto& unit : problem.get_units()) {
      std::println("schedule for unit {}:", unit.get_id());

      for (size_t t = 0; t < H; ++t) {
        for (const auto* task : unit.get_tasks() | std::views::keys) {
          Field x =
              builder.extract_variable(point, starts.at({&unit, task, t}));
          Field Q =
              builder.extract_variable(point, starts.at({&unit, task, t}));

          std::println("x({}, {}, {}) = {}", unit.get_id(), task->get_id(), t,
                       x);
          std::println("Q({}, {}, {}) = {}", unit.get_id(), task->get_id(), t,
                       Q);
        }
      }
    }
  }
}
