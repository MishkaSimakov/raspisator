#include "linear/matrix/NPY.h"
#include "linear/problem/MPS.h"
#include "linear/problem/optimization/TransformToEqualities.h"
using Field = double;

int main(int argc, char* argv[]) {
  if (argc != 3) {
    throw std::invalid_argument("Wrong number of arguments. Expected 2.");
  }

  std::string input_path = argv[1];
  std::string output_path = argv[2];

  MPSReader<Field> reader(MPSFieldsMode::SPACE_SEPARATED);
  reader.read(input_path);

  auto problem = reader.get_canonical_representation();
  auto optimized = TransformToEqualities<Field>().apply(problem);

  std::ofstream output(output_path);
  linalg::to_npy(output, to_matrices(optimized).A);
}
