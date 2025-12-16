#pragma once
#include "model/STN.h"

// figure 2 in Blomer et al., 2000
template <typename Field>
STN<Field> small_blomer_problem(Field target1, Field target2) {
  STN<Field> stn;

  // States
  auto* state1 = stn.add(InputState<Field>{1000});
  auto* state2 = stn.add(InputState<Field>{1000});
  auto* state3 = stn.add(InputState<Field>{1000});
  auto* state4 = stn.add(NormalState<Field>{0, 0, 100});
  auto* state5 = stn.add(NormalState<Field>{0, 0, 200});
  auto* state6 = stn.add(NonStorableState<Field>{});
  auto* state7 = stn.add(NormalState<Field>{0, 0, 150});
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

// figure 3 in Blomer et al., 2000
template <typename Field>
STN<Field> big_blomer_problem(Field target1, Field target2, Field target3,
                              Field target4, Field target5) {
  STN<Field> stn;

  // States
  auto* state0 = stn.add(InputState<Field>{1000});
  auto* state1 = stn.add(NormalState<Field>{10, 0, 30});
  auto* state2 = stn.add(NormalState<Field>{10, 0, 30});
  auto* state3 = stn.add(NormalState<Field>{0, 0, 15});
  auto* state4 = stn.add(NormalState<Field>{10, 0, 30});
  auto* state5 = stn.add(NonStorableState<Field>{});
  auto* state6 = stn.add(NormalState<Field>{0, 0, 10});
  auto* state7 = stn.add(NormalState<Field>{0, 0, 10});
  auto* state8 = stn.add(NormalState<Field>{0, 0, 10});
  auto* state9 = stn.add(NonStorableState<Field>{});
  auto* state10 = stn.add(NonStorableState<Field>{});
  auto* state11 = stn.add(NormalState<Field>{0, 0, 10});
  auto* state12 = stn.add(NonStorableState<Field>{});
  auto* state13 = stn.add(NormalState<Field>{0, 0, 10});
  auto* state14 = stn.add(OutputState<Field>{0, target1});
  auto* state15 = stn.add(OutputState<Field>{0, target2});
  auto* state16 = stn.add(OutputState<Field>{0, target3});
  auto* state17 = stn.add(OutputState<Field>{0, target4});
  auto* state18 = stn.add(OutputState<Field>{0, target5});

  // Units
  auto* unit1 = stn.add(Unit<Field>{});
  auto* unit2 = stn.add(Unit<Field>{});
  auto* unit3 = stn.add(Unit<Field>{});
  auto* unit4 = stn.add(Unit<Field>{});
  auto* unit5 = stn.add(Unit<Field>{});
  auto* unit6 = stn.add(Unit<Field>{});
  auto* unit7 = stn.add(Unit<Field>{});
  auto* unit8 = stn.add(Unit<Field>{});
  auto* unit9 = stn.add(Unit<Field>{});

  // Tasks
  auto* task1 = stn.add(Task<Field>{});
  task1->add_input(state0, 1);
  task1->add_output(state1, 1);
  unit1->attach_task(task1, {2, 3, 10});

  auto* task2 = stn.add(Task<Field>{});
  task2->add_input(state1, 1);
  task2->add_output(state2, Field(2) / 5);
  task2->add_output(state3, Field(3) / 5);
  unit2->attach_task(task2, {4, 5, 20});

  auto* task3 = stn.add(Task<Field>{});
  task3->add_input(state3, 1);
  task3->add_output(state1, Field(31) / 100);
  task3->add_output(state4, Field(69) / 100); // this is Blomer, not me!
  unit3->attach_task(task3, {2, 4, 10});

  auto* task4 = stn.add(Task<Field>{});
  task4->add_input(state2, 1);
  task4->add_output(state5, 1);
  unit4->attach_task(task4, {4, 4, 10});

  auto* task5 = stn.add(Task<Field>{});
  task5->add_input(state2, 1);
  task5->add_output(state6, 1);
  unit4->attach_task(task5, {4, 4, 10});

  auto* task6 = stn.add(Task<Field>{});
  task6->add_input(state4, 1);
  task6->add_output(state7, 1);
  unit4->attach_task(task6, {4, 4, 10});

  auto* task7 = stn.add(Task<Field>{});
  task7->add_input(state4, 1);
  task7->add_output(state8, 1);
  unit4->attach_task(task7, {4, 4, 10});

  auto* task8 = stn.add(Task<Field>{});
  task8->add_input(state2, 1);
  task8->add_output(state9, 1);
  unit5->attach_task(task8, {6, 4, 10});

  auto* task9 = stn.add(Task<Field>{});
  task9->add_input(state4, 1);
  task9->add_output(state10, 1);
  unit5->attach_task(task9, {6, 4, 10});

  auto* task10 = stn.add(Task<Field>{});
  task10->add_input(state6, 1);
  task10->add_output(state11, 1);
  unit6->attach_task(task10, {4, 3, 7});
  unit7->attach_task(task10, {5, 3, 7});

  auto* task11 = stn.add(Task<Field>{});
  task11->add_input(state7, 1);
  task11->add_output(state12, 1);
  unit6->attach_task(task11, {5, 3, 7});
  unit7->attach_task(task11, {6, 3, 7});

  auto* task12 = stn.add(Task<Field>{});
  task12->add_input(state8, 1);
  task12->add_output(state13, 1);
  unit6->attach_task(task12, {6, 3, 7});
  unit7->attach_task(task12, {6, 3, 7});

  auto* task13 = stn.add(Task<Field>{});
  task13->add_input(state9, 1);
  task13->add_output(state14, 1);
  unit8->attach_task(task13, {4, 4, 12});
  unit9->attach_task(task13, {6, 4, 12});

  auto* task14 = stn.add(Task<Field>{});
  task14->add_input(state10, 1);
  task14->add_output(state15, 1);
  unit8->attach_task(task14, {4, 4, 12});
  unit9->attach_task(task14, {6, 4, 12});

  auto* task15 = stn.add(Task<Field>{});
  task15->add_input(state5, Field(1) / 2);
  task15->add_input(state11, Field(1) / 2);
  task15->add_output(state16, 1);
  unit8->attach_task(task15, {4, 4, 12});

  auto* task16 = stn.add(Task<Field>{});
  task16->add_input(state12, 1);
  task16->add_output(state17, 1);
  unit8->attach_task(task16, {6, 4, 12});
  unit9->attach_task(task16, {6, 4, 12});

  auto* task17 = stn.add(Task<Field>{});
  task17->add_input(state13, 1);
  task17->add_output(state18, 1);
  unit8->attach_task(task17, {6, 4, 12});
  unit9->attach_task(task17, {6, 4, 12});

  return stn;
}
