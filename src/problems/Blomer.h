#pragma once
#include "model/STN.h"

// figure 2 in Blomer et al., 2000
inline STN generate_problem() {
  STN stn;

  // States
  auto* state1 = stn.add(InputState{});
  auto* state2 = stn.add(InputState{});
  auto* state3 = stn.add(InputState{});
  auto* state4 = stn.add(NormalState{0, 0, 100});
  auto* state5 = stn.add(NormalState{0, 0, 200});
  auto* state6 = stn.add(NonStorableState{});
  auto* state7 = stn.add(NormalState{0, 0, 150});
  auto* state8 = stn.add(NonStorableState{});
  auto* state9 = stn.add(OutputState{});
  auto* state10 = stn.add(OutputState{});

  // Units
  auto* unit1 = stn.add(Unit{});
  auto* unit2 = stn.add(Unit{});
  auto* unit3 = stn.add(Unit{});
  auto* unit4 = stn.add(Unit{});

  // Tasks
  auto* task1 = stn.add(Task{});
  task1->add_inputs(state1);
  task1->add_outputs(state4);
  unit1->attach_task(task1, {1, 0, 100});

  auto* task2 = stn.add(Task{});
  task2->add_inputs(state4, state5);
  task2->add_outputs(state7, state9);
  unit2->attach_task(task2, {2, 0, 50});
  unit3->attach_task(task2, {2, 0, 90});

  auto* task3 = stn.add(Task{});
  task3->add_inputs(state2, state3);
  task3->add_outputs(state5);
  unit2->attach_task(task3, {2, 0, 50});
  unit3->attach_task(task3, {2, 0, 90});

  auto* task4 = stn.add(Task{});
  task4->add_inputs(state3, state7);
  task4->add_outputs(state6);
  unit2->attach_task(task4, {1, 0, 50});
  unit3->attach_task(task4, {1, 0, 90});

  auto* task5 = stn.add(Task{});
  task5->add_inputs(state6);
  task5->add_outputs(state8, state10);
  unit4->attach_task(task5, {1, 0, 200});

  auto* task6 = stn.add(Task{});
  task6->add_inputs(state8);
  task6->add_outputs(state7);
  unit4->attach_task(task6, {2, 0, 200});

  return stn;
}
