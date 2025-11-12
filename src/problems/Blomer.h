#pragma once
#include "model/STN.h"

// NOLINTBEGIN(readability-magic-numbers)
// figure 2 in Blomer et al., 2000
inline STN generate_problem() {
  STN stn;

  // states
  auto* state1 = stn.add(InputState{});
  auto* state2 = stn.add(InputState{});
  auto* state3 = stn.add(InputState{});
  auto* state4 = stn.add(NormalState{0, 0, 100});
  auto* state5 = stn.add(NormalState{0, 0, 200});
  auto* state6 = stn.add(NonStorableState{});

  // units
  auto* unit1 = stn.add(Unit{});
  auto* unit2 = stn.add(Unit{});
  auto* unit3 = stn.add(Unit{});
  auto* unit4 = stn.add(Unit{});

  auto* task1 = stn.add(Task{});
  task1->add_inputs(state1);
  task1->add_outputs(state4);
  unit1->attach_task(task1, {1, 0, 100});

  auto* task2 = stn.add(Task{});
  task2->add_inputs(state4, state5);
  unit2->attach_task(task2, {2, 0, 50});
  unit3->attach_task(task2, {2, 0, 90});

  return stn;
}

// NOLINTEND(readability-magic-numbers)
