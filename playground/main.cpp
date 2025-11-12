#include <iostream>

#include "../src/utils/Drawing.h"
#include "model/STN.h"
#include "problems/Blomer.h"

int main() {
  auto problem = generate_problem();

  std::cout << to_graphviz(problem) << std::endl;

  return 0;
}
