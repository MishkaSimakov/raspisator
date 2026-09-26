#include <gtest/gtest.h>

#include "schedule/model/Solution.h"
#include "support/GMPRational.h"

STN<GMPRational>* generate_sequential_problem() {
  STN<GMPRational>* stn = new STN<GMPRational>{};

  auto* state1 = stn->add(InputState<GMPRational>{100});
  auto* state2 = stn->add(NormalState<GMPRational>{0, 0, 100});
  auto* state3 = stn->add(NonStorableState<GMPRational>{});
  auto* state4 = stn->add(OutputState<GMPRational>{0, 50});

  auto* unit0 = stn->add(Unit<GMPRational>{});

  auto* task0 = stn->add(Task<GMPRational>{});
  task0->add_input(state1, 1);
  task0->add_output(state2, 1);
  unit0->attach_task(task0, {1, 0, 10});

  auto* task1 = stn->add(Task<GMPRational>{});
  task1->add_input(state2, 1);
  task1->add_output(state3, 1);
  unit0->attach_task(task1, {2, 0, 10});

  auto* task2 = stn->add(Task<GMPRational>{});
  task2->add_input(state3, 1);
  task2->add_output(state4, 1);
  unit0->attach_task(task2, {5, 0, 10});

  return stn;
}

class SolutionCheckerTests : public ::testing::Test {
 protected:
  STN<GMPRational>* sequential;

  void SetUp() override { sequential = generate_sequential_problem(); }

  void TearDown() override { delete sequential; }
};

TEST_F(SolutionCheckerTests, AllOk) {
  Solution<GMPRational> solution(sequential);

  solution.add_instance({0, 0, 10, 0});
  solution.add_instance({0, 0, 10, 1});
  solution.add_instance({0, 0, 10, 2});
  solution.add_instance({0, 0, 10, 3});
  solution.add_instance({0, 0, 10, 4});

  solution.add_instance({0, 1, 10, 5});
  solution.add_instance({0, 2, 10, 7});
  solution.add_instance({0, 1, 10, 12});
  solution.add_instance({0, 2, 10, 14});
  solution.add_instance({0, 1, 10, 19});
  solution.add_instance({0, 2, 10, 21});
  solution.add_instance({0, 1, 10, 26});
  solution.add_instance({0, 2, 10, 28});
  solution.add_instance({0, 1, 10, 33});
  solution.add_instance({0, 2, 10, 35});

  ASSERT_TRUE(solution.check());
}

TEST_F(SolutionCheckerTests, InvalidInstance) {
  Solution<GMPRational> solution(sequential);

  solution.add_instance({0, 0, 10, 0});
  solution.add_instance({0, 100500, 10, 1});
  solution.add_instance({0, 0, 10, 2});
  solution.add_instance({0, 0, 10, 3});
  solution.add_instance({0, 0, 10, 4});

  solution.add_instance({0, 1, 10, 5});
  solution.add_instance({0, 2, 10, 7});
  solution.add_instance({0, 1, 10, 12});
  solution.add_instance({0, 2, 10, 14});
  solution.add_instance({0, 1, 10, 19});
  solution.add_instance({0, 2, 10, 21});
  solution.add_instance({0, 1, 10, 26});
  solution.add_instance({0, 2, 10, 28});
  solution.add_instance({0, 1, 10, 33});
  solution.add_instance({0, 2, 10, 35});

  ASSERT_FALSE(solution.check());

  Solution<GMPRational> solution2(sequential);

  solution2.add_instance({0, 0, 10, 0});
  solution2.add_instance({0, 0, 10, 1});
  solution2.add_instance({100, 0, 10, 2});
  solution2.add_instance({0, 0, 10, 3});
  solution2.add_instance({0, 0, 10, 4});

  solution2.add_instance({0, 1, 10, 5});
  solution2.add_instance({0, 2, 10, 7});
  solution2.add_instance({0, 1, 10, 12});
  solution2.add_instance({0, 2, 10, 14});
  solution2.add_instance({0, 1, 10, 19});
  solution2.add_instance({0, 2, 10, 21});
  solution2.add_instance({0, 1, 10, 26});
  solution2.add_instance({0, 2, 10, 28});
  solution2.add_instance({0, 1, 10, 33});
  solution2.add_instance({0, 2, 10, 35});

  ASSERT_FALSE(solution2.check());
}

TEST_F(SolutionCheckerTests, StateOverflow) {
  Solution<GMPRational> solution(sequential);

  solution.add_instance({0, 0, 10, 0});
  solution.add_instance({0, 0, 10, 1});
  solution.add_instance({0, 0, 10, 2});
  solution.add_instance({0, 0, 10, 3});
  solution.add_instance({0, 0, 10, 4});
  solution.add_instance({0, 0, 10, 5});
  solution.add_instance({0, 0, 10, 6});
  solution.add_instance({0, 0, 10, 7});
  solution.add_instance({0, 0, 10, 8});
  solution.add_instance({0, 0, 10, 9});
  solution.add_instance({0, 0, 10, 10});

  solution.add_instance({0, 1, 10, 12});
  solution.add_instance({0, 2, 10, 14});
  solution.add_instance({0, 1, 10, 19});
  solution.add_instance({0, 2, 10, 21});
  solution.add_instance({0, 1, 10, 26});
  solution.add_instance({0, 2, 10, 28});
  solution.add_instance({0, 1, 10, 33});
  solution.add_instance({0, 2, 10, 35});
  solution.add_instance({0, 1, 10, 40});
  solution.add_instance({0, 2, 10, 42});

  ASSERT_FALSE(solution.check());
}

TEST_F(SolutionCheckerTests, BigBatchSize) {
  Solution<GMPRational> solution(sequential);

  solution.add_instance({0, 0, 10, 0});
  solution.add_instance({0, 0, 10, 1});
  solution.add_instance({0, 0, 10, 2});
  solution.add_instance({0, 0, 10, 3});
  solution.add_instance({0, 0, 15, 4});

  solution.add_instance({0, 1, 10, 5});
  solution.add_instance({0, 2, 10, 7});
  solution.add_instance({0, 1, 10, 12});
  solution.add_instance({0, 2, 10, 14});
  solution.add_instance({0, 1, 10, 19});
  solution.add_instance({0, 2, 10, 21});
  solution.add_instance({0, 1, 10, 26});
  solution.add_instance({0, 2, 10, 28});
  solution.add_instance({0, 1, 10, 33});
  solution.add_instance({0, 2, 10, 35});

  ASSERT_FALSE(solution.check());
}

TEST_F(SolutionCheckerTests, BusyUnit) {
  Solution<GMPRational> solution(sequential);

  solution.add_instance({0, 0, 10, 0});
  solution.add_instance({0, 0, 10, 1});
  solution.add_instance({0, 0, 10, 2});
  solution.add_instance({0, 0, 10, 3});
  solution.add_instance({0, 0, 10, 4});

  solution.add_instance({0, 1, 10, 5});
  solution.add_instance({0, 2, 10, 6});
  solution.add_instance({0, 1, 10, 11});
  solution.add_instance({0, 2, 10, 13});
  solution.add_instance({0, 1, 10, 18});
  solution.add_instance({0, 2, 10, 20});
  solution.add_instance({0, 1, 10, 25});
  solution.add_instance({0, 2, 10, 27});
  solution.add_instance({0, 1, 10, 32});
  solution.add_instance({0, 2, 10, 34});

  ASSERT_FALSE(solution.check());
}

TEST_F(SolutionCheckerTests, TargetNotAcquired) {
  Solution<GMPRational> solution(sequential);

  solution.add_instance({0, 0, 10, 0});
  solution.add_instance({0, 0, 10, 1});
  solution.add_instance({0, 0, 10, 2});
  solution.add_instance({0, 0, 10, 3});
  solution.add_instance({0, 0, 10, 4});

  solution.add_instance({0, 1, 10, 5});
  solution.add_instance({0, 2, 10, 7});
  solution.add_instance({0, 1, 10, 12});
  solution.add_instance({0, 2, 10, 14});
  solution.add_instance({0, 1, 10, 19});
  solution.add_instance({0, 2, 10, 21});
  solution.add_instance({0, 1, 10, 26});
  solution.add_instance({0, 2, 10, 28});

  ASSERT_FALSE(solution.check());
}
