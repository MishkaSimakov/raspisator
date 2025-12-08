#pragma once
#include "model/STN.h"

// figure 2 in Blomer et al., 2000
template <typename Field>
STN<Field> small_blomer_problem(Field target1, Field target2) {
  STN<Field> stn;

  // States
  auto* state1 = stn.add(InputState{1000});
  auto* state2 = stn.add(InputState{1000});
  auto* state3 = stn.add(InputState{1000});
  auto* state4 = stn.add(NormalState{0, 0, 100});
  auto* state5 = stn.add(NormalState{0, 0, 200});
  auto* state6 = stn.add(NonStorableState<Field>{});
  auto* state7 = stn.add(NormalState{0, 0, 150});
  auto* state8 = stn.add(NonStorableState<Field>{});
  auto* state9 = stn.add(OutputState<Field>{0, target1});
  auto* state10 = stn.add(OutputState<Field>{0, target2});

  // Units
  auto* unit1 = stn.add(Unit<Field>{});
  auto* unit2 = stn.add(Unit<Field>{});
  auto* unit3 = stn.add(Unit<Field>{});
  auto* unit4 = stn.add(Unit<Field>{});

  // Tasks
  auto* task1 = stn.add(Task<Field>{});
  task1->add_input(state1, 1);
  task1->add_output(state4, 1);
  unit1->attach_task(task1, {1, 0, 100});

  auto* task2 = stn.add(Task<Field>{});
  task2->add_input(state4, Field{2} / 5);
  task2->add_input(state5, Field{3} / 5);
  task2->add_output(state9, Field{2} / 5);
  task2->add_output(state7, Field{3} / 5);
  unit2->attach_task(task2, {2, 0, 50});
  unit3->attach_task(task2, {2, 0, 90});

  auto* task3 = stn.add(Task<Field>{});
  task3->add_input(state2, Field{1} / 2);
  task3->add_input(state3, Field{1} / 2);
  task3->add_output(state5, 1);
  unit2->attach_task(task3, {2, 0, 50});
  unit3->attach_task(task3, {2, 0, 90});

  auto* task4 = stn.add(Task<Field>{});
  task4->add_input(state3, Field{1} / 5);
  task4->add_input(state7, Field{4} / 5);
  task4->add_output(state6, 1);
  unit2->attach_task(task4, {1, 0, 50});
  unit3->attach_task(task4, {1, 0, 90});

  auto* task5 = stn.add(Task<Field>{});
  task5->add_input(state6, 1);
  task5->add_output(state8, Field{1} / 10);
  task5->add_output(state10, Field{9} / 10);
  unit4->attach_task(task5, {1, 0, 200});

  auto* task6 = stn.add(Task<Field>{});
  task6->add_input(state8, 1);
  task6->add_output(state7, 1);
  unit4->attach_task(task6, {2, 0, 200});

  return stn;
}
