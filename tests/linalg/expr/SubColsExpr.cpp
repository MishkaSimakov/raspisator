#include <gtest/gtest.h>

#include "linalg/CSCMatrix.h"
#include "linalg/expr/SubColsExpr.h"

using namespace linalg;

static_assert(ColWiseMatrixRange<
              detail::SubColsExpr<CSCMatrix<int>, std::span<size_t>>>);
