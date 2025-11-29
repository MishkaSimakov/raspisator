#include <iostream>
#include <print>
#include <random>

#include "linear/BigInteger.h"
#include "linear/MPS.h"
#include "linear/SimplexMethod.h"
#include "linear/bb/BranchAndBound.h"
#include "linear/bb/Drawer.h"
#include "linear/builder/ProblemBuilder.h"
#include "linear/sparse/CSCMatrix.h"
#include "linear/sparse/LU.h"
#include "problems/Dwarf.h"
#include "utils/Drawing.h"
#include "utils/Hashers.h"

using Field = double;

template <typename Field>
CSCMatrix<Field> sparse_matrix_1(size_t N, size_t density_multiplier,
                                 size_t seed = 0) {
  auto A = Matrix<Field>::unity(N);

  std::mt19937 engine(seed);
  std::uniform_int_distribution<size_t> uniform_dist(0, N - 1);

  // add O(N) non-zeros in "random" places
  for (size_t i = 0; i < density_multiplier * N; ++i) {
    size_t row = uniform_dist(engine);
    size_t col = uniform_dist(engine);

    A[row, col] = i + 1;
  }

  return CSCMatrix(A);
}

double max_difference(const CSCMatrix<double>& left,
                      const CSCMatrix<Rational>& right) {
  auto dense_left = linalg::to_dense(left);
  auto dense_right = linalg::to_dense(right);

  auto [n, d] = dense_right.shape();
  double max_difference = 0;

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < d; ++j) {
      double difference =
          std::abs(static_cast<double>(dense_right[i, j]) - dense_left[i, j]);

      max_difference = std::max(difference, max_difference);
    }
  }

  return max_difference;
}

int main() {
  size_t density_multiplier = 3;

  std::println("size,error");

  for (size_t N = 10; N < 75; ++N) {
    for (size_t i = 0; i < 10; ++i) {
      auto A_rational = sparse_matrix_1<Rational>(N, density_multiplier, i);
      auto A_double = sparse_matrix_1<double>(N, density_multiplier, i);

      auto [rational_L, rational_U, rational_P] = linalg::sparse_lu(A_rational);
      auto [double_L, double_U, double_P] = linalg::sparse_lu(A_double);

      double error = std::max(max_difference(double_L, rational_L),
                              max_difference(double_U, rational_U));

      std::println("{},{}", N, error);
    }
  }

  // auto matrix = sparse_matrix_1<Rational>(5, 1);
  //
  // std::cout << matrix << std::endl;
  //
  // auto [L, U, P] = linalg::sparse_lu(matrix);
  //
  // auto dense_L = linalg::to_dense(L);
  // auto dense_U = linalg::to_dense(U);
  //
  // for (size_t i = 0; i < dense_L.get_height(); ++i) {
  //   dense_L[i, i] = 1;
  // }
  //
  // std::cout << "L:\n" << dense_L << std::endl;
  // std::cout << "U:\n" << dense_U << std::endl;
  //
  // std::cout << "LU:\n" << (dense_L * dense_U) << std::endl;
  //
  // std::cout << "P: ";
  // for (size_t i : P) {
  //   std::cout << i << " ";
  // }
  // std::cout << std::endl;
  //
  // std::cout << "P^-1 LU:\n" << linalg::apply_permutation(CSCMatrix(dense_L * dense_U), P)
  //           << std::endl;

  // std::filesystem::remove_all("matrices");
  // std::filesystem::create_directory("matrices");
  //
  // auto problem = dwarf_problem_normal<Field>();
  //
  // std::cout << to_graphviz(problem) << std::endl;
  //
  // size_t H = 8;  // number of periods
  //
  // ProblemBuilder<Field> builder;
  //
  // // decision variables
  // auto makespan = builder.new_variable("MS", VariableType::INTEGER);
  // std::unordered_map<std::tuple<const Unit<Field>*, const Task<Field>*,
  // size_t>,
  //                    Variable<Field>>
  //     quantities;
  // std::unordered_map<std::tuple<const Unit<Field>*, const Task<Field>*,
  // size_t>,
  //                    Variable<Field>>
  //     starts;
  // std::unordered_map<std::pair<const State*, size_t>, Variable<Field>>
  // stocks;
  //
  // for (const auto& unit : problem.get_units()) {
  //   for (const auto* task : unit.get_tasks() | std::views::keys) {
  //     for (size_t t = 0; t < H; ++t) {
  //       quantities.emplace(
  //           std::tuple{&unit, task, t},
  //           builder.new_variable(
  //               std::format("Q({}, {}, {})", unit.get_id(), task->get_id(),
  //               t), VariableType::INTEGER));
  //       starts.emplace(
  //           std::tuple{&unit, task, t},
  //           builder.new_variable(
  //               std::format("x({}, {}, {})", unit.get_id(), task->get_id(),
  //               t), VariableType::INTEGER));
  //     }
  //   }
  // }
  //
  // for (const State& state : problem.get_states()) {
  //   if (!std::holds_alternative<NonStorableState>(state)) {
  //     for (size_t t = 0; t < H; ++t) {
  //       stocks.emplace(
  //           std::pair{&state, t},
  //           builder.new_variable(fmt::format("p({}, {})", state.get_id(), t),
  //                                VariableType::INTEGER));
  //     }
  //   }
  // }
  //
  // // makespan
  // for (const auto& unit : problem.get_units()) {
  //   for (const auto& [task, props] : unit.get_tasks()) {
  //     for (size_t t = 0; t < H; ++t) {
  //       Field finish_time(t + props.batch_processing_time);
  //       builder.add_constraint(makespan >=
  //                              finish_time * starts.at({&unit, task, t}));
  //     }
  //   }
  // }
  //
  // // batch size limits
  // for (const auto& unit : problem.get_units()) {
  //   for (const auto& [task, props] : unit.get_tasks()) {
  //     for (size_t t = 0; t < H; ++t) {
  //       builder.add_constraint(Field(props.batch_min_size) *
  //                                  starts.at({&unit, task, t}) <=
  //                              quantities.at({&unit, task, t}));
  //
  //       builder.add_constraint(quantities.at({&unit, task, t}) <=
  //                              Field(props.batch_max_size) *
  //                                  starts.at({&unit, task, t}));
  //     }
  //   }
  // }
  //
  // // stock balance
  // for (const State& state : problem.get_states()) {
  //   if (std::holds_alternative<NonStorableState>(state)) {
  //     continue;
  //   }
  //
  //   for (size_t t = 0; t < H; ++t) {
  //     auto new_stock =
  //         t != 0 ? Expression(stocks.at({&state, t - 1}))
  //                : std::visit(Overload{[](const auto& state) {
  //                                        return Field(state.initial_stock);
  //                                      },
  //                                      [](NonStorableState) -> Field {
  //                                        std::unreachable();
  //                                      }},
  //                             state);
  //
  //     for (const auto* task :
  //          problem.get_producers_of(state) | std::views::elements<0>) {
  //       for (const auto& [unit, props] : problem.get_task_units(*task)) {
  //         if (t >= props.batch_processing_time) {
  //           new_stock +=
  //               quantities.at({unit, task, t - props.batch_processing_time});
  //         }
  //       }
  //     }
  //
  //     for (const auto& [task, fraction] : problem.get_consumers_of(state)) {
  //       Expression<Field> task_consumption(0);
  //
  //       for (const auto& [unit, props] : problem.get_task_units(*task)) {
  //         task_consumption += quantities.at({unit, task, t});
  //       }
  //
  //       new_stock -= task_consumption * fraction;
  //     }
  //
  //     // TODO: external demand + external supply
  //
  //     builder.add_constraint(new_stock == stocks.at({&state, t}));
  //   }
  // }
  //
  // // stock limits
  // for (size_t t = 0; t < H; ++t) {
  //   for (const State& state : problem.get_states()) {
  //     if (!std::holds_alternative<NormalState>(state)) {
  //       continue;
  //     }
  //
  //     Expression<Field> max_stock(std::get<NormalState>(state).max_level);
  //     builder.add_constraint(stocks.at({&state, t}) <= max_stock);
  //   }
  // }
  //
  // // production of non-storable goods
  // for (const State& state : problem.get_states()) {
  //   if (!std::holds_alternative<NonStorableState>(state)) {
  //     continue;
  //   }
  //
  //   for (size_t t = 0; t < H; ++t) {
  //     Expression<Field> production(0);
  //
  //     for (const auto& [task, fraction] : problem.get_producers_of(state)) {
  //       Expression<Field> task_production(0);
  //
  //       for (const auto& [unit, props] : problem.get_task_units(*task)) {
  //         if (t >= props.batch_processing_time) {
  //           task_production +=
  //               quantities.at({unit, task, t - props.batch_processing_time});
  //         }
  //       }
  //
  //       production += task_production * fraction;
  //     }
  //
  //     Expression<Field> consumption(0);
  //
  //     for (const auto& [task, fraction] : problem.get_consumers_of(state)) {
  //       Expression<Field> task_consumption(0);
  //
  //       for (const auto& [unit, props] : problem.get_task_units(*task)) {
  //         task_consumption += quantities.at({unit, task, t});
  //       }
  //
  //       consumption += task_consumption * fraction;
  //     }
  //
  //     builder.add_constraint(consumption == production);
  //   }
  // }
  //
  // // assigning batches to production units
  // for (const auto& unit : problem.get_units()) {
  //   for (size_t t = 0; t < H; ++t) {
  //     Expression<Field> running_batches_cnt(0);
  //
  //     for (const auto& [task, props] : unit.get_tasks()) {
  //       if (props.batch_processing_time >= t) {
  //         for (size_t t2 = t - props.batch_processing_time; t2 < t; ++t2) {
  //           running_batches_cnt += starts.at({&unit, task, t2});
  //         }
  //       }
  //     }
  //
  //     builder.add_constraint(running_batches_cnt <= Expression<Field>(1));
  //   }
  // }
  //
  // // variables domain
  // for (const auto& unit : problem.get_units()) {
  //   for (const auto* task : unit.get_tasks() | std::views::keys) {
  //     for (size_t t = 0; t < H; ++t) {
  //       builder.add_constraint(starts.at({&unit, task, t}) <=
  //                              Expression<Field>(1));
  //     }
  //   }
  // }
  //
  // // desired amounts
  // for (const State& state : problem.get_states()) {
  //   if (!std::holds_alternative<OutputState>(state)) {
  //     continue;
  //   }
  //
  //   Expression<Field> desired_amount(std::get<OutputState>(state).target);
  //
  //   builder.add_constraint(stocks.at({&state, H - 1}) >= desired_amount);
  // }
  //
  // // objective
  // builder.set_objective(-makespan);
  //
  // std::cout << builder << std::endl;
  //
  // // solve MILP problem
  // auto milp_problem = builder.get_problem();
  //
  // auto solver = BranchAndBound<Field, SimplexMethod<Field>>(milp_problem);
  // auto solution = solver.solve();
  //
  // std::cout << GraphvizBuilder<Field>().build(solver.get_tree()) <<
  // std::endl;
  //
  // if (std::holds_alternative<NoFiniteSolution>(solution)) {
  //   std::println("No finite solution.");
  // } else {
  //   auto finite_solution = std::get<FiniteMILPSolution<Field>>(solution);
  //
  //   std::println("Finish production in {} time units.\n",
  //                -finite_solution.value);
  //
  //   auto point = finite_solution.point;
  //
  //   for (const auto& unit : problem.get_units()) {
  //     std::println("schedule for unit {}:", unit.get_id());
  //
  //     for (size_t t = 0; t < H; ++t) {
  //       for (const auto* task : unit.get_tasks() | std::views::keys) {
  //         Field x =
  //             builder.extract_variable(point, starts.at({&unit, task, t}));
  //         Field Q =
  //             builder.extract_variable(point, quantities.at({&unit, task,
  //             t}));
  //
  //         std::println("x({}, {}, {}) = {}", unit.get_id(), task->get_id(),
  //         t,
  //                      x);
  //         std::println("Q({}, {}, {}) = {}", unit.get_id(), task->get_id(),
  //         t,
  //                      Q);
  //       }
  //     }
  //   }
  // }
}
