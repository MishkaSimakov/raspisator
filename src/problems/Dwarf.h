#pragma once

#include "linear/BigInteger.h"
#include "model/STN.h"

// a miniscule problem straight out of my head
template <typename Field>
inline STN<Field> generate_problem() {
  STN<Field> stn;

  // States
  auto* state1 = stn.add(InputState{});
  auto* state2 = stn.add(NormalState{0, 0, 50});
  auto* state3 = stn.add(NormalState{0, 0, 50});
  auto* state4 = stn.add(OutputState{});

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
  task2->add_input(state2, Rational{1} / 2);
  task2->add_input(state3, Rational{1} / 2);
  task2->add_output(state4, 1);
  unit1->attach_task(task3, {1, 50, 200});

  return stn;
}
