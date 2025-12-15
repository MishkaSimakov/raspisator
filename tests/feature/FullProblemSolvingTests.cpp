#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <vector>

#include "encoding/UniformTimeDiscretization.h"
#include "linear/bb/FullStrongBranching.h"
#include "linear/bb/PseudoCost.h"
#include "linear/bb/Settings.h"
#include "linear/problem/ToMatrices.h"
#include "linear/problem/optimization/FullOptimizer.h"
#include "model/STN.h"
#include "model/Solution.h"
#include "problems/Dwarf.h"

auto test_problems() {
  std::vector<std::pair<std::string, STN<double>>> problems;

  problems.emplace_back("dwarf_smaller", dwarf_problem_smaller<double>());
  problems.emplace_back("dwarf_small", dwarf_problem_small<double>());
  problems.emplace_back("dwarf_normal_100", dwarf_problem_normal<double>(100));
  problems.emplace_back("dwarf_normal_200", dwarf_problem_normal<double>(200));
  problems.emplace_back("dwarf_normal_300", dwarf_problem_normal<double>(300));

  std::vector<std::pair<std::string, std::shared_ptr<STN<double>>>> result;

  for (auto& [name, stn] : problems) {
    result.emplace_back(name, std::make_shared<STN<double>>(std::move(stn)));
  }

  return result;
}

class FullProblemSolvingTests
    : public ::testing::TestWithParam<
          std::pair<std::string, std::shared_ptr<STN<double>>>> {};

TEST_P(FullProblemSolvingTests, SolveThenCheck) {
  auto stn = GetParam().second;

  size_t H = 10;

  auto encoding = to_uniform_time_milp(*stn, H);

  // solve MILP problem
  auto optimized = FullOptimizer<double>().apply(encoding.builder);

  auto matrices = to_matrices(optimized);

  auto settings = BranchAndBoundSettings<double>{
      .max_nodes = 100'000,
      .perturbation = PerturbationMode::DISABLED,
  };
  auto solver = FullStrongBranchingBranchAndBound(matrices.A, matrices.b, matrices.c,
                                         matrices.lower, matrices.upper,
                                         matrices.variables, settings);
  auto solution = solver.solve();

  // check solution
  ASSERT_TRUE(std::holds_alternative<FiniteMILPSolution<double>>(solution));
  auto finite_solution = std::get<FiniteMILPSolution<double>>(solution);

  auto point = finite_solution.point;

  Solution checker(stn.get());

  for (const auto& unit : stn->get_units()) {
    for (size_t t = 0; t < H; ++t) {
      for (const auto* task : unit.get_tasks() | std::views::keys) {
        double x = matrices.extract_variable(
            encoding.starts.at({&unit, task, t}), point);
        double Q = matrices.extract_variable(
            encoding.quantities.at({&unit, task, t}), point);

        if (FieldTraits<double>::is_strictly_positive(x)) {
          checker.add_instance(
              TaskInstance{unit.get_id(), task->get_id(), Q, t});
        }
      }
    }
  }

  ASSERT_TRUE(checker.check());
}
INSTANTIATE_TEST_SUITE_P(
    FullProblemSolvingTests, FullProblemSolvingTests,
    ::testing::ValuesIn(test_problems()),
    [](const testing::TestParamInfo<FullProblemSolvingTests::ParamType>& info) {
      return info.param.first;
    });
