#include "linear/matrix/NPY.h"
#include "linear/problem/MPS.h"
#include "linear/problem/optimization/TransformToEqualities.h"
using Field = double;

int main() {
  std::string filepath;
  std::cin >> filepath;

  MPSReader<Field> reader(MPSFieldsMode::SPACE_SEPARATED);
  reader.read(filepath);

  auto problem = reader.get_canonical_representation();
  auto optimized = TransformToEqualities<Field>().apply(problem);

  std::ofstream output("output.npy");
  linalg::to_npy(output, to_matrices(optimized).A);
}
