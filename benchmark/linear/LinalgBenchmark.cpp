#include <benchmark/benchmark.h>

#include <random>

#include "linalg/CSCMatrix.h"
#include "linalg/Matrix.h"
#include "linalg/Random.h"
#include "linalg/Transpose.h"
#include "linalg/Vector.h"

using namespace linalg;

static size_t N = 10'000;

// Must always return the same matrices so that benchmark results are reliable.
auto get_matrices() {
  std::default_random_engine random(0);
  std::uniform_real_distribution<double> value_distribution(1, 10);

  auto matrix = random::sparse(N, N, 10, random, value_distribution);

  Vector cost =
      Matrix<double>::generate(N, 1, [](size_t i, size_t j) { return i + j; });
  Vector pi =
      Matrix<double>::generate(N, 1, [](size_t i, size_t j) { return i + j; });

  return std::tuple{std::move(matrix), std::move(cost), std::move(pi)};
}

static void NewLinalgLibrary(benchmark::State& state) {
  const auto [matrix, cost, pi] = get_matrices();

  for (auto _ : state) {
    Matrix result = cost - transpose(matrix) * pi;

    benchmark::DoNotOptimize(result);
    benchmark::ClobberMemory();
  }
}

static void OldWay(benchmark::State& state) {
  const auto [matrix, cost, pi] = get_matrices();

  for (auto _ : state) {
    auto result = Matrix<double>::zeros(N, 1);

    for (size_t i = 0; i < N; ++i) {
      result[i, 0] = cost[i, 0];

      for (const auto& [row, value] : matrix.get_column(i)) {
        result[i, 0] -= value * pi[row, 0];
      }
    }

    benchmark::DoNotOptimize(result);
    benchmark::ClobberMemory();
  }
}

BENCHMARK(NewLinalgLibrary);
BENCHMARK(OldWay);

BENCHMARK_MAIN();
