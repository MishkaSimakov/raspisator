#pragma once

#include "model/STN.h"

// a miniscule problems straight out of my head
template <typename Field>
STN<Field> dwarf_problem_normal() {
  STN<Field> stn;

  // States
  auto* state1 = stn.add(InputState{200});
  auto* state2 = stn.add(NormalState{0, 0, 50});
  auto* state3 = stn.add(NormalState{0, 0, 50});
  auto* state4 = stn.add(OutputState{0, 100});

  // Units
  auto* unit1 = stn.add(Unit<Field>{});
  auto* unit2 = stn.add(Unit<Field>{});

  // Tasks
  auto* task1 = stn.add(Task<Field>{});
  task1->add_input(state1, 1);
  task1->add_output(state2, 1);
  unit1->attach_task(task1, {1, 5, 30});

  auto* task2 = stn.add(Task<Field>{});
  task2->add_input(state1, 1);
  task2->add_output(state3, 1);
  unit2->attach_task(task2, {1, 5, 30});

  auto* task3 = stn.add(Task<Field>{});
  task3->add_input(state2, Field{1} / 2);
  task3->add_input(state3, Field{1} / 2);
  task3->add_output(state4, 1);
  unit1->attach_task(task3, {1, 50, 200});

  return stn;
}

template <typename Field>
STN<Field> dwarf_problem_small() {
  STN<Field> stn;

  // States
  auto* state1 = stn.add(InputState{100});
  auto* state2 = stn.add(NormalState{0, 0, 100});
  auto* state3 = stn.add(OutputState{0, 100});

  // Units
  auto* unit1 = stn.add(Unit<Field>{});
  auto* unit2 = stn.add(Unit<Field>{});

  // Tasks
  auto* task1 = stn.add(Task<Field>{});
  task1->add_input(state1, 1);
  task1->add_output(state2, 1);
  unit1->attach_task(task1, {1, 10, 200});

  auto* task2 = stn.add(Task<Field>{});
  task2->add_input(state2, 1);
  task2->add_output(state3, 1);
  unit2->attach_task(task2, {1, 10, 200});

  return stn;
}

// even smaller!
// https://vk.com/wall-166124324_5897
template <typename Field>
STN<Field> dwarf_problem_smaller() {
  STN<Field> stn;

  // States
  auto* state1 = stn.add(InputState{100});
  auto* state2 = stn.add(OutputState{0, 100});

  // Units
  auto* unit1 = stn.add(Unit<Field>{});

  // Tasks
  auto* task1 = stn.add(Task<Field>{});
  task1->add_input(state1, 1);
  task1->add_output(state2, 1);
  unit1->attach_task(task1, {1, 10, 200});

  return stn;
}
