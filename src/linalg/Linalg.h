#pragma once

// Umbrella header that includes main linalg types. Also adds usings for most
// commonly used types.

// Core classes
#include "linalg/CSCMatrix.h"
#include "linalg/Matrix.h"
#include "linalg/Vector.h"

// Methods for matrices
#include "Arithmetics.h"
#include "Print.h"
#include "Transpose.h"

using linalg::Matrix, linalg::Vector, linalg::CSCMatrix;
